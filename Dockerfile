FROM python:3.12-slim

ENV PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1

WORKDIR /app

COPY requirements.txt ./
RUN pip install --no-cache-dir -r requirements.txt

COPY app ./app

RUN useradd --system --uid 10001 demux \
    && mkdir -p /app/data \
    && chown -R demux:demux /app/data
USER demux

# 所有状态(SQLite、暂存区、发布区)都在这一个数据目录里
ENV DATA_DIR=/app/data

# 输出用有界句柄池,数十个 fd 即可处理数百个样本;这里按测试条件收紧到 64,
# 防止镜像默认的高上限掩盖句柄泄露。docker-compose.yml 也设置了同一限额。
# 注意:ulimit 是每进程设置,不能写在 Dockerfile 里,由 compose/运行参数传入。
ENV DEMUX_MAX_OPEN_OUTPUTS=16

EXPOSE 8000
# 单 worker:后台线程在本进程内处理,多 worker 会与中断恢复逻辑冲突
CMD ["uvicorn", "app.main:app", "--host", "0.0.0.0", "--port", "8000"]
