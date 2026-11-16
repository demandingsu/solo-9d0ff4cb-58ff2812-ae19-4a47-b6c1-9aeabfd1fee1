"""路径配置:所有状态(数据库、暂存区、发布区)都在 DATA_DIR 下,保证同一数据卷。"""
import os
from pathlib import Path


def data_dir() -> Path:
    env = os.environ.get("DATA_DIR")
    base = Path(env) if env else Path(__file__).resolve().parent.parent / "data"
    return base.resolve()


def staging_dir() -> Path:
    """处理中的工作区:上传件与半成品输出都在这里,绝不直接对外提供。"""
    return data_dir() / "staging"


def published_dir() -> Path:
    """发布区:只有完整处理成功的运行才会出现在这里。"""
    return data_dir() / "published"


def db_path() -> Path:
    return data_dir() / "app.db"


# 同时打开的拆样输出文件上限。三份输入另占 3 个句柄,默认值在 nofile=64
# 的容器里也有充足余量;可用环境变量按需调整。
DEFAULT_MAX_OPEN_OUTPUTS = 16


def max_open_outputs() -> int:
    """DEMUX_MAX_OPEN_OUTPUTS 覆盖输出句柄上限;缺失或非法时用默认值。"""
    raw = os.environ.get("DEMUX_MAX_OPEN_OUTPUTS", "")
    try:
        value = int(raw)
    except (TypeError, ValueError):
        return DEFAULT_MAX_OPEN_OUTPUTS
    return value if value >= 1 else DEFAULT_MAX_OPEN_OUTPUTS
