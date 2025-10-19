FROM python:3.11

WORKDIR /app
COPY requirements.txt .

RUN pip install --no-cache-dir -r requirements.txt

ARG MODEL_URL="https://huggingface.co/TheBloke/openchat-3.5-1210-GGUF/resolve/4b8bafa8b8a69336dca6e306a129c55f1ebbac05/openchat-3.5-1210.Q3_K_S.gguf"
ENV MODEL_PATH="/app/model/openchat-3.5-1210.Q3_K_S.gguf"

RUN mkdir -p /app/model \
    && curl -L "$MODEL_URL" -o "$MODEL_PATH"

COPY . .

EXPOSE 8000

CMD ["uvicorn", "main:app", "--host", "0.0.0.0", "--port", "8000"]
