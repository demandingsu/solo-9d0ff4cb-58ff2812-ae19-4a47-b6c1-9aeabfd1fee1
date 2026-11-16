"""逐对流式拆分核心:一次只处理一条记录,不整份读入内存。

输入为三份未压缩的 4 行 FASTQ(R1/R2/I1)。任何结构损坏、读段号不一致、
序列与质量长度不符都会抛 FastqFormatError,由上层把运行标记为失败,
保证不发布半套文件。

输出句柄策略:样本数可能有数百个,不能为每个样本常开 R1/R2 两个文件。
所有输出文件先以空文件建好(零读段组也有文件),写记录时通过有界 LRU
句柄池按需以追加方式打开,最多同时占用 ``max_open_outputs`` 个句柄。
同一分组始终按读段到达顺序追加,因此中途被换出再打开也不会打乱顺序。
"""
from __future__ import annotations

import os
from collections import OrderedDict
from itertools import zip_longest

BASES = "ACGT"
_MISSING = object()  # zip_longest 哨兵:某份文件先读完

# 默认最多同时打开 16 个输出文件:3 个输入 + 16 个输出,在 nofile=64
# 的容器里仍留出足够余量给解释器、监听套接字和数据库连接。
# 实际取值由 config.max_open_outputs() 从 DEMUX_MAX_OPEN_OUTPUTS 读取,
# 再经 demultiplex(..., max_open_outputs=...) 传入。
OUTPUT_BUFFER_BYTES = 1 << 16


class FastqFormatError(Exception):
    """FASTQ 结构损坏或三份文件读对不一致。"""


def read_id(header: str) -> str:
    """读段号 = 头部第一个空白之前的部分,如 '@INST:1:FC:1:1101:100:100 1:N:0:0'。"""
    parts = header.split(None, 1)
    return parts[0] if parts else ""


def iter_fastq(path, label):
    """按 4 行一条产出 (header, seq, plus, qual);结构损坏即抛错。"""
    rec = 0
    with open(path, "r", encoding="ascii", newline="") as fh:
        while True:
            try:
                block = [fh.readline() for _ in range(4)]
            except UnicodeDecodeError as exc:
                raise FastqFormatError(
                    f"{label}: 第 {rec + 1} 条附近不是纯文本(需要未压缩 FASTQ): {exc}"
                ) from exc
            if block[0] == "":
                return  # 恰好结束在记录边界
            rec += 1
            if block[1] == "" or block[2] == "" or block[3] == "":
                raise FastqFormatError(f"{label}: 第 {rec} 条记录在文件末尾被截断")
            header, seq, plus, qual = (line.rstrip("\r\n") for line in block)
            if not header.startswith("@"):
                raise FastqFormatError(f"{label}: 第 {rec} 条记录头缺少 '@'")
            if not plus.startswith("+"):
                raise FastqFormatError(f"{label}: 第 {rec} 条记录分隔行缺少 '+'")
            if len(seq) != len(qual):
                raise FastqFormatError(
                    f"{label}: 第 {rec} 条({header})序列长度 {len(seq)} 与质量长度 {len(qual)} 不一致"
                )
            yield header[1:], seq, plus, qual


def build_matcher(samples):
    """预计算精确匹配表与单碱基错配候选表(错配串 -> 候选样本集合)。"""
    exact = {}
    neighbors: dict[str, set[str]] = {}
    for sid, index in samples:
        exact[index] = sid
    for sid, index in samples:
        for pos in range(len(index)):
            for base in BASES:
                if base != index[pos]:
                    mut = index[:pos] + base + index[pos + 1 :]
                    neighbors.setdefault(mut, set()).add(sid)
    return exact, neighbors


def assign_index(index, exact, neighbors):
    """精确匹配优先;否则仅当单碱基错配候选唯一时纠错;其余返回未分配。"""
    sid = exact.get(index)
    if sid is not None:
        return sid, False
    cands = neighbors.get(index)
    if cands is not None and len(cands) == 1:
        return next(iter(cands)), True
    return None, False


def _default_max_open_outputs() -> int:
    """环境变量可覆盖句柄上限;非法值回落到默认值。"""
    from .config import max_open_outputs

    return max_open_outputs()


def _prepare_outputs(out_dir, order):
    """为每个分组预建空的 R1/R2 文件,返回 {组: (r1路径, r2路径)}。

    逐个用 with 打开后立刻关闭,不会累积句柄;这样零读段组也一定有文件。
    """
    paths = {}
    for sid in order:
        r1 = os.path.join(out_dir, f"{sid}_R1.fastq")
        r2 = os.path.join(out_dir, f"{sid}_R2.fastq")
        for path in (r1, r2):
            with open(path, "w", encoding="ascii"):
                pass
        paths[sid] = (r1, r2)
    return paths


class _OutputPool:
    """有界 LRU 输出句柄池:超出上限时关闭最久未用的文件。

    ``write`` 中的打开/写入/换出错误都会直接抛给调用方;调用方在失败路径
    用 ``close_quietly`` 收尾,避免关闭错误掩盖最初的异常。
    """

    def __init__(self, max_open: int):
        self._max_open = max(1, int(max_open))
        self._handles: OrderedDict[str, object] = OrderedDict()
        self.peak = 0

    def write(self, path: str, text: str) -> None:
        fh = self._handles.get(path)
        if fh is None:
            if len(self._handles) >= self._max_open:
                _, stale = self._handles.popitem(last=False)
                stale.close()  # 换出前落盘;写失败在此暴露
            fh = open(
                path, "a", encoding="ascii", buffering=OUTPUT_BUFFER_BYTES
            )
            self._handles[path] = fh
            self.peak = max(self.peak, len(self._handles))
        else:
            self._handles.move_to_end(path)
        fh.write(text)

    def close(self) -> None:
        """关闭全部句柄;若有多个关闭失败,报第一个,其余仍会关闭。"""
        first_error = None
        while self._handles:
            _, fh = self._handles.popitem(last=False)
            try:
                fh.close()
            except Exception as exc:  # noqa: BLE001 - 需要关闭其余句柄后再抛
                if first_error is None:
                    first_error = exc
        if first_error is not None:
            raise first_error

    def close_quietly(self) -> None:
        while self._handles:
            _, fh = self._handles.popitem(last=False)
            try:
                fh.close()
            except Exception:  # noqa: BLE001 - 失败路径只做清理
                pass


def demultiplex(r1_path, r2_path, i1_path, samples, out_dir,
                max_open_outputs=None):
    """逐对拆分。samples 为 [(样本号, 索引序列)],输出写入 out_dir。

    返回 {'total_pairs', 'counts': {组: [对数, 纠错对数]}, 'order': [组...],
    'peak_open_outputs': 峰值句柄数}。
    输出文件按输入顺序追加,保证各组内读段顺序与原始一致。
    """
    exact, neighbors = build_matcher(samples)
    order = [sid for sid, _ in samples] + ["undetermined"]
    counts = {sid: [0, 0] for sid in order}
    paths = _prepare_outputs(out_dir, order)
    if max_open_outputs is None:
        max_open_outputs = _default_max_open_outputs()
    pool = _OutputPool(max_open_outputs)
    readers = {
        "R1": iter_fastq(r1_path, "R1"),
        "R2": iter_fastq(r2_path, "R2"),
        "I1": iter_fastq(i1_path, "I1"),
    }
    total = 0
    try:
        streams = zip_longest(*readers.values(), fillvalue=_MISSING)
        for rec1, rec2, reci in streams:
            if rec1 is _MISSING or rec2 is _MISSING or reci is _MISSING:
                raise FastqFormatError(f"第 {total + 1} 条起,三份文件的记录数不一致")
            total += 1
            h1, s1, p1, q1 = rec1
            h2, s2, p2, q2 = rec2
            hi, si, _pi, _qi = reci
            id1, id2, idi = read_id(h1), read_id(h2), read_id(hi)
            if not id1 or id1 != id2 or id1 != idi:
                raise FastqFormatError(
                    f"第 {total} 条读段号不一致: R1={h1!r} R2={h2!r} I1={hi!r}"
                )
            sid, corrected = assign_index(si.upper(), exact, neighbors)
            key = sid if sid is not None else "undetermined"
            w1, w2 = paths[key]
            pool.write(w1, f"@{h1}\n{s1}\n{p1}\n{q1}\n")
            pool.write(w2, f"@{h2}\n{s2}\n{p2}\n{q2}\n")
            counts[key][0] += 1
            if corrected:
                counts[key][1] += 1
    except BaseException:
        pool.close_quietly()  # 不掩盖原始异常
        raise
    else:
        pool.close()  # 所有缓冲在此落盘;写入失败会抛给上层,阻止发布
    finally:
        for reader in readers.values():
            reader.close()
    return {
        "total_pairs": total,
        "counts": counts,
        "order": order,
        "peak_open_outputs": pool.peak,
    }
