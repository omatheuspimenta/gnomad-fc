# docker build -t aprovar-import:1.0 -t aprovar-import:latest .
FROM python:3.11-slim

# Dependências
RUN apt-get update && \
    apt-get install -y --no-install-recommends openjdk-21-jdk curl && \
    rm -rf /var/lib/apt/lists/*

# Python
RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir hail

COPY export_to_es.py /app/

WORKDIR /data

CMD ["python3.11"]