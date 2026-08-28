FROM python:3.11-slim

# Install uv
COPY --from=ghcr.io/astral-sh/uv:latest /uv /bin/uv

WORKDIR /app

# Create non-root user
RUN useradd -m -u 1000 appuser

# Copy files
COPY pyproject.toml uv.lock README.md ./
COPY api.py ./
COPY app ./app
COPY gunicorn_conf.py ./

# Change ownership BEFORE installing deps
RUN chown -R appuser:appuser /app

# Switch user BEFORE uv sync
USER appuser

# Install dependencies
RUN uv sync --frozen

EXPOSE 8000

CMD ["uv", "run", "gunicorn", "-c", "gunicorn_conf.py", "api:app"]