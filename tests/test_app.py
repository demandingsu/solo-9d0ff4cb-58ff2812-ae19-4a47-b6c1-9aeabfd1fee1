"""关键情形自测:拆分正确性、纠错规则、损坏处理、原子发布、恢复、不可改写。"""
import sys
import time
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from app.demux import (  # noqa: E402
    FastqFormatError,
    assign_index,
    build_matcher,
    iter_fastq,
    read_id,
)
from app.samples import SampleSheetError, load_samples  # noqa: E402

SAMPLES_CSV = "sample_id,index\nS1,AAAACCCC\nS2,AAAACCGG\nS3,GGGGTTTT\n"

# S1=AAAACCCC 与 S2=AAAACCGG 相距 2,便于构造歧义错配;S3 远离二者
READS = [
    ("r001", "AAAACCCC"),  # S1 精确
    ("r002", "AAAACCGG"),  # S2 精确
    ("r003", "AAAACCCA"),  # 距 S1 唯一单碱基错配 -> 纠错到 S1
    ("r004", "AAAACCGC"),  # 距 S1、S2 各 1 -> 歧义 -> 未分配
    ("r005", "NNNNNNNN"),  # 无候选 -> 未分配
    ("r006", "GGGGTTTA"),  # 距 S3 唯一单碱基错配 -> 纠错到 S3
    ("r007", "AAAACCCC"),  # S1 精确(检验同组内顺序)
]

SEQ1, SEQ2 = "ACGT" * 10, "TGCA" * 10
QUAL = "I" * 40


def fq(records):
    return "".join(f"@{h}\n{s}\n+\n{q}\n" for h, s, q in records)


def make_inputs(reads=READS):
    r1, r2, i1 = [], [], []
    for n, (rid, idx) in enumerate(reads):
        head = f"INST:1:FC:1:1101:{1000 + n}:{2000 + n}"
        r1.append((f"{head} 1:N:0:0", SEQ1, QUAL))
        r2.append((f"{head} 2:N:0:0", SEQ2, QUAL))
        i1.append((f"{head} 1:N:0:0", idx, "I" * len(idx)))
    return fq(r1).encode(), fq(r2).encode(), fq(i1).encode(), SAMPLES_CSV.encode()


@pytest.fixture()
def client(tmp_path, monkeypatch):
    monkeypatch.setenv("DATA_DIR", str(tmp_path / "data"))
    from fastapi.testclient import TestClient
    from app.main import app

    with TestClient(app) as c:
        c.data_dir = tmp_path / "data"
        yield c


def submit(client, run_number, r1, r2, i1, samples):
    return client.post(
        "/runs",
        data={"run_number": run_number},
        files={
            "r1": ("r1.fastq", r1),
            "r2": ("r2.fastq", r2),
            "i1": ("i1.fastq", i1),
            "samples": ("samples.csv", samples),
        },
        follow_redirects=False,
    )


def submit_good(client, run_number="RUN1"):
    r1, r2, i1, csv = make_inputs()
    return submit(client, run_number, r1, r2, i1, csv)


def run_id_of(resp):
    assert resp.status_code == 303, resp.text
    return int(resp.headers["location"].rsplit("/", 1)[1])


def wait_done(client, run_id, timeout=15):
    deadline = time.time() + timeout
    while time.time() < deadline:
        st = client.get(f"/api/runs/{run_id}").json()
        if st["status"] != "processing":
            return st
        time.sleep(0.05)
    raise TimeoutError("处理超时")


def parse_fastq_ids(text):
    lines = text.splitlines()
    return [lines[i][1:].split()[0] for i in range(0, len(lines), 4)]


def head_of(n):
    """第 n 条 READS 记录对应的读段号(见 make_inputs)。"""
    return f"INST:1:FC:1:1101:{1000 + n}:{2000 + n}"


def counts_of(status):
    return {c["sample_id"]: (c["pairs"], c["corrected"]) for c in status["counts"]}


# ------------------------------------------------------------ 单元级

def test_assign_index_rules():
    samples = [("S1", "AAAACCCC"), ("S2", "AAAACCGG"), ("S3", "GGGGTTTT")]
    exact, neighbors = build_matcher(samples)
    assert assign_index("AAAACCCC", exact, neighbors) == ("S1", False)   # 精确优先
    assert assign_index("AAAACCCA", exact, neighbors) == ("S1", True)    # 唯一错配
    assert assign_index("AAAACCGC", exact, neighbors) == (None, False)   # 歧义
    assert assign_index("NNNNNNNN", exact, neighbors) == (None, False)   # 无候选
    assert assign_index("AAAA", exact, neighbors) == (None, False)       # 长度不符


def test_read_id():
    assert read_id("INST:1:FC:1:1101:100:100 1:N:0:0") == "INST:1:FC:1:1101:100:100"
    assert read_id("PLAINID") == "PLAINID"
    assert read_id("") == ""


def test_iter_fastq_truncated(tmp_path):
    p = tmp_path / "bad.fastq"
    p.write_text("@a\nACGT\n+\n")  # 缺质量行
    with pytest.raises(FastqFormatError, match="截断"):
        list(iter_fastq(str(p), "R1"))


def test_iter_fastq_len_mismatch(tmp_path):
    p = tmp_path / "bad.fastq"
    p.write_text("@a\nACGT\n+\nIII\n")  # 质量短 1
    with pytest.raises(FastqFormatError, match="长度"):
        list(iter_fastq(str(p), "R1"))


def test_samples_validation(tmp_path):
    def load(text):
        p = tmp_path / "s.csv"
        p.write_text(text)
        return load_samples(str(p))

    assert load("sid,idx\nA,ACGT\n") == [("A", "ACGT")]
    assert load("A,acgt\n") == [("A", "ACGT")]  # 小写归一、无表头
    with pytest.raises(SampleSheetError, match="重复"):
        load("A,ACGT\nB,ACGT\n")
    with pytest.raises(SampleSheetError, match="长度不一致"):
        load("A,ACGT\nB,ACGTA\n")
    with pytest.raises(SampleSheetError, match="非法字符"):
        load("A,ACGT\nB,ACGB\n")
    with pytest.raises(SampleSheetError, match="保留字"):
        load("undetermined,ACGT\n")


# ------------------------------------------------------------ 端到端

def test_happy_path(client):
    rid = run_id_of(submit_good(client))
    st = wait_done(client, rid)
    assert st["status"] == "published", st["error"]
    assert st["total_pairs"] == 7
    assert counts_of(st) == {
        "S1": (3, 1), "S2": (1, 0), "S3": (1, 1), "undetermined": (2, 0),
    }

    s1_r1 = client.get(f"/runs/{rid}/download/S1_R1.fastq")
    assert s1_r1.status_code == 200
    assert parse_fastq_ids(s1_r1.text) == [head_of(0), head_of(2), head_of(6)]  # 原顺序
    s1_r2 = client.get(f"/runs/{rid}/download/S1_R2.fastq")
    assert parse_fastq_ids(s1_r2.text) == [head_of(0), head_of(2), head_of(6)]
    assert SEQ2 in s1_r2.text  # R2 内容确实是 R2
    und = client.get(f"/runs/{rid}/download/undetermined_R1.fastq")
    assert parse_fastq_ids(und.text) == [head_of(3), head_of(4)]

    log = client.get(f"/runs/{rid}/download/demux.log")
    assert log.status_code == 200 and "总读段对: 7" in log.text

    page = client.get(f"/runs/{rid}")
    assert page.status_code == 200 and "AAAACCCC" in page.text and "已发布" in page.text


def test_truncated_input_fails_then_retry_succeeds(client):
    r1, r2, i1, csv = make_inputs()
    bad_r2 = r2[: r2.rfind(b"\n", 0, r2.rfind(b"\n"))]  # 末条记录缺质量行
    rid = run_id_of(submit(client, "RUN2", r1, bad_r2, i1, csv))
    st = wait_done(client, rid)
    assert st["status"] == "failed"
    assert "截断" in st["error"]
    # 不发布半套文件
    assert not (client.data_dir / "published" / "RUN2").exists()
    assert client.get(f"/runs/{rid}/download/S1_R1.fastq").status_code == 409

    # 同一运行号重试可完成
    rid2 = run_id_of(submit_good(client, "RUN2"))
    st2 = wait_done(client, rid2)
    assert st2["status"] == "published"
    assert st2["total_pairs"] == 7
    assert client.get(f"/runs/{rid2}/download/S1_R1.fastq").status_code == 200


def test_read_id_mismatch_fails(client):
    r1, r2, i1, csv = make_inputs()
    bad_r2 = r2.replace(b"1101:1001:2001", b"1101:9999:2001", 1)
    rid = run_id_of(submit(client, "RUN3", r1, bad_r2, i1, csv))
    st = wait_done(client, rid)
    assert st["status"] == "failed" and "读段号不一致" in st["error"]
    assert not (client.data_dir / "published" / "RUN3").exists()


def test_seq_qual_len_mismatch_fails(client):
    r1, r2, i1, csv = make_inputs()
    bad_i1 = i1.replace(b"IIIIIIII\n", b"IIIIIII\n", 1)  # 首条 I1 质量少 1
    rid = run_id_of(submit(client, "RUN4", r1, r2, bad_i1, csv))
    st = wait_done(client, rid)
    assert st["status"] == "failed" and "长度" in st["error"]


def test_record_count_mismatch_fails(client):
    r1, r2, i1, csv = make_inputs()
    extra = f"@INST:1:FC:1:1101:9:9 1:N:0:0\n{SEQ1}\n+\n{QUAL}\n".encode()
    rid = run_id_of(submit(client, "RUN5", r1 + extra, r2, i1, csv))
    st = wait_done(client, rid)
    assert st["status"] == "failed" and "记录数不一致" in st["error"]


def test_published_run_immutable(client):
    rid = run_id_of(submit_good(client, "RUN6"))
    assert wait_done(client, rid)["status"] == "published"
    before = client.get(f"/runs/{rid}/download/S1_R1.fastq").text

    resp = submit_good(client, "RUN6")  # 同号再提交
    assert resp.status_code == 409
    assert "不可改写" in resp.json()["detail"]
    after = client.get(f"/runs/{rid}/download/S1_R1.fastq").text
    assert before == after  # 文件未被改写


def test_bad_sample_sheet_rejected(client):
    r1, r2, i1, _ = make_inputs()
    dup = b"A,ACGT\nB,ACGT\n"
    resp = submit(client, "RUN7", r1, r2, i1, dup)
    assert resp.status_code == 400
    assert "重复" in resp.json()["detail"]


def test_interrupted_run_recovered_on_startup(client):
    """模拟进程中断:遗留 processing 记录 + 半成品目录,重启后标失败、清理干净。"""
    from app.main import recover_interrupted
    from app.db import db_conn
    from app.config import db_path, staging_dir, published_dir

    staging = staging_dir() / "deadbeef"
    (staging / "out").mkdir(parents=True)
    (staging / "out" / "S1_R1.fastq").write_text("半成品")
    stray_pub = published_dir() / "RUN8"
    stray_pub.mkdir(parents=True)
    (stray_pub / "S1_R1.fastq").write_text("改名后来不及落库的文件")
    with db_conn(db_path()) as conn:
        conn.execute(
            "INSERT INTO runs(run_number, status, staging_dir, created_at, updated_at) "
            "VALUES ('RUN8', 'processing', ?, 't', 't')",
            (str(staging),),
        )

    recover_interrupted()  # 模拟重启

    with db_conn(db_path()) as conn:
        row = conn.execute("SELECT * FROM runs WHERE run_number='RUN8'").fetchone()
    assert row["status"] == "failed" and "中断" in row["error"]
    assert not staging.exists()
    assert not stray_pub.exists()
    # 恢复后同一运行号可重试并完成
    rid = run_id_of(submit_good(client, "RUN8"))
    assert wait_done(client, rid)["status"] == "published"


def test_index_page(client):
    resp = client.get("/")
    assert resp.status_code == 200
    assert "新建拆样运行" in resp.text


# ------------------------------------------------------------ 句柄有界与交错

def _many_samples(n, length=8):
    """生成 n 个长度唯一、仅含 ACGT 的索引(按四进制编码样本号)。"""
    samples, seen = [], set()
    for i in range(n):
        x, chars = i, []
        for _ in range(length):
            chars.append("ACGT"[x % 4])
            x //= 4
        index = "".join(chars)
        assert index not in seen
        seen.add(index)
        samples.append((f"S{i:04d}", index))
    return samples


def _write_interleaved_inputs(dir_path, samples, pairs):
    """第 r 对循环落到不同样本,使输出句柄反复换入换出。"""
    expected = {}
    with open(dir_path / "r1.fastq", "w") as r1, \
         open(dir_path / "r2.fastq", "w") as r2, \
         open(dir_path / "i1.fastq", "w") as i1:
        for r in range(pairs):
            sid, idx = samples[(r * 7) % len(samples)]
            if r % 13 == 0:  # 少量无候选读段进未分配组
                sid, idx = "undetermined", "T" * len(idx)
            head = f"INST:1:FC:1:1101:{r}:{r}"
            r1.write(f"@{head} 1:N:0:0\nACGTACGTAC\n+\nIIIIIIIIII\n")
            r2.write(f"@{head} 2:N:0:0\nTGCATGCATG\n+\nIIIIIIIIII\n")
            i1.write(f"@{head} 1:N:0:0\n{idx}\n+\n{'I' * len(idx)}\n")
            expected[sid] = expected.get(sid, 0) + 1
    return expected


def test_bounded_output_handles_interleaved(tmp_path):
    """数百样本反复交错时,输出句柄峰值受控、顺序与计数不变、空组也有文件。"""
    from app.demux import demultiplex

    samples = _many_samples(300)
    expected = _write_interleaved_inputs(tmp_path, samples, pairs=1500)
    out = tmp_path / "out"
    out.mkdir()

    stats = demultiplex(
        str(tmp_path / "r1.fastq"), str(tmp_path / "r2.fastq"),
        str(tmp_path / "i1.fastq"), samples, str(out), max_open_outputs=4,
    )

    assert stats["total_pairs"] == 1500
    assert stats["peak_open_outputs"] <= 4  # 与句柄上限严格一致
    assert stats["order"][-1] == "undetermined"
    for sid, (pairs, _corrected) in stats["counts"].items():
        assert pairs == expected.get(sid, 0), sid
        r1 = out / f"{sid}_R1.fastq"
        r2 = out / f"{sid}_R2.fastq"
        assert r1.is_file() and r2.is_file()  # 零读段组也有文件
        ids = parse_fastq_ids(r1.read_text())
        nums = [int(x.rsplit(":", 1)[1]) for x in ids]
        assert len(ids) == pairs
        assert nums == sorted(nums), sid  # 组内保持输入顺序


def test_under_nofile_64(tmp_path):
    """把进程 fd 上限压到 64,数百样本仍能完成拆分(原实现会 EMFILE)。"""
    import resource
    from app.demux import demultiplex

    soft, hard = resource.getrlimit(resource.RLIMIT_NOFILE)
    try:
        limit = 64 if hard in (resource.RLIM_INFINITY, -1) else min(64, hard)
        resource.setrlimit(resource.RLIMIT_NOFILE, (limit, hard))
        samples = _many_samples(400)
        expected = _write_interleaved_inputs(tmp_path, samples, pairs=2000)
        out = tmp_path / "out"
        out.mkdir()
        stats = demultiplex(
            str(tmp_path / "r1.fastq"), str(tmp_path / "r2.fastq"),
            str(tmp_path / "i1.fastq"), samples, str(out),
        )
    finally:
        resource.setrlimit(resource.RLIMIT_NOFILE, (soft, hard))
    assert stats["total_pairs"] == 2000
    assert stats["peak_open_outputs"] <= 16
    assert sum(v[0] for v in stats["counts"].values()) == 2000
    assert all(stats["counts"][sid][0] == pairs for sid, pairs in expected.items())
    assert len(list(out.iterdir())) == 2 * len(stats["order"])


def test_tail_corruption_and_retry_consistency(client):
    """末尾损坏 → 不发布;同号重试后页面计数与下载内容一致。"""
    r1, r2, i1, csv = make_inputs()
    # 破坏末条 I1 的质量行:改成长度不符,前面记录全部合法
    bad_i1 = i1.rstrip(b"\n")
    bad_i1 = bad_i1[: bad_i1.rfind(b"\n")] + b"\nIII\n"
    rid = run_id_of(submit(client, "RUN9", r1, r2, bad_i1, csv))
    st = wait_done(client, rid)
    assert st["status"] == "failed" and "长度" in st["error"]
    assert not (client.data_dir / "published" / "RUN9").exists()
    assert not any((client.data_dir / "staging").iterdir())
    assert client.get(f"/runs/{rid}/download/S1_R1.fastq").status_code == 409

    rid2 = run_id_of(submit_good(client, "RUN9"))
    st2 = wait_done(client, rid2)
    assert st2["status"] == "published"
    page = client.get(f"/runs/{rid2}").text
    for sid, (pairs, corrected) in counts_of(st2).items():
        assert f"<td>{sid}</td>" in page
        r1_text = client.get(f"/runs/{rid2}/download/{sid}_R1.fastq").text
        r2_text = client.get(f"/runs/{rid2}/download/{sid}_R2.fastq").text
        assert len(parse_fastq_ids(r1_text)) == pairs
        assert len(parse_fastq_ids(r2_text)) == pairs
        assert parse_fastq_ids(r1_text) == parse_fastq_ids(r2_text)
        assert f"<td>{pairs}</td>" in page
        assert f"<td>{corrected}</td>" in page


def test_write_failure_does_not_publish(client, monkeypatch):
    """输出写入在末尾失败(模拟磁盘写满):不发布半套文件,重试可成功。"""
    from app.demux import _OutputPool

    real_write = _OutputPool.write
    calls = {"n": 0}

    def failing_write(self, path, text):
        calls["n"] += 1
        if calls["n"] >= 9:  # 前几条正常,模拟处理到尾段时写盘失败
            raise OSError(28, "No space left on device")
        return real_write(self, path, text)

    monkeypatch.setattr(_OutputPool, "write", failing_write)
    rid = run_id_of(submit_good(client, "RUN10"))
    st = wait_done(client, rid)
    assert st["status"] == "failed" and "No space left" in st["error"]
    assert not (client.data_dir / "published" / "RUN10").exists()
    assert not any((client.data_dir / "staging").iterdir())

    monkeypatch.setattr(_OutputPool, "write", real_write)  # 修复后同号重试
    rid2 = run_id_of(submit_good(client, "RUN10"))
    st2 = wait_done(client, rid2)
    assert st2["status"] == "published"
    assert counts_of(st2) == {
        "S1": (3, 1), "S2": (1, 0), "S3": (1, 1), "undetermined": (2, 0),
    }
    for sid, (pairs, _c) in counts_of(st2).items():
        assert len(parse_fastq_ids(
            client.get(f"/runs/{rid2}/download/{sid}_R1.fastq").text)) == pairs


def test_flush_failure_on_close_does_not_publish(client, monkeypatch):
    """末尾缓冲落盘失败(close 时写盘):运行标失败、不发布、暂存清理干净。"""
    from app.demux import _OutputPool

    real_close = _OutputPool.close

    def failing_close(self):
        real_close(self)  # 句柄确实都关了,不掩盖后续异常
        raise OSError(28, "No space left on device")

    monkeypatch.setattr(_OutputPool, "close", failing_close)
    rid = run_id_of(submit_good(client, "RUN13"))
    st = wait_done(client, rid)
    assert st["status"] == "failed" and "No space left" in st["error"]
    assert not (client.data_dir / "published" / "RUN13").exists()
    assert not any((client.data_dir / "staging").iterdir())

    monkeypatch.setattr(_OutputPool, "close", real_close)
    rid2 = run_id_of(submit_good(client, "RUN13"))
    st2 = wait_done(client, rid2)
    assert st2["status"] == "published" and st2["total_pairs"] == 7


def test_zero_read_groups_have_empty_files(client):
    """从未命中的样本发布后仍有空 R1/R2 文件可下载,计数为 0。"""
    r1, r2, i1, _csv = make_inputs()
    csv = b"sample_id,index\nS1,AAAACCCC\nS2,AAAACCGG\nS3,GGGGTTTT\nS4,TTTTGGGG\n"
    rid = run_id_of(submit(client, "RUN11", r1, r2, i1, csv))
    st = wait_done(client, rid)
    assert st["status"] == "published"
    assert counts_of(st)["S4"] == (0, 0)
    for mate in ("R1", "R2"):
        resp = client.get(f"/runs/{rid}/download/S4_{mate}.fastq")
        assert resp.status_code == 200 and resp.content == b""


def test_api_many_samples_under_nofile_64(client):
    """端到端:API 上传 120 个样本、读段交错,在 nofile=64 下发布完整结果。"""
    import resource

    samples = _many_samples(120)
    tmp = Path(client.data_dir.parent) / "inputs"
    tmp.mkdir()
    expected = _write_interleaved_inputs(tmp, samples, pairs=600)
    csv = "".join(f"{sid},{idx}\n" for sid, idx in samples).encode()
    r1 = (tmp / "r1.fastq").read_bytes()
    r2 = (tmp / "r2.fastq").read_bytes()
    i1 = (tmp / "i1.fastq").read_bytes()

    soft, hard = resource.getrlimit(resource.RLIMIT_NOFILE)
    try:
        limit = 64 if hard in (resource.RLIM_INFINITY, -1) else min(64, hard)
        resource.setrlimit(resource.RLIMIT_NOFILE, (limit, hard))
        rid = run_id_of(submit(client, "RUN12", r1, r2, i1, csv))
        st = wait_done(client, rid)
    finally:
        resource.setrlimit(resource.RLIMIT_NOFILE, (soft, hard))

    assert st["status"] == "published", st["error"]
    assert st["total_pairs"] == 600
    got = counts_of(st)
    assert got.get("undetermined") == (expected.get("undetermined", 0), 0)
    assert sum(p for p, _c in got.values()) == 600
    for sid, _idx in samples:
        assert got[sid][0] == expected.get(sid, 0)
        assert len(parse_fastq_ids(
            client.get(f"/runs/{rid}/download/{sid}_R1.fastq").text)) == got[sid][0]
    log = client.get(f"/runs/{rid}/download/demux.log").text
    assert "峰值同时打开输出文件:" in log
