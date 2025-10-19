FROM python:3.11

COPY requirements.txt .

RUN pip install -r requirements.txt

WORKDIR /app
COPY . /app

EXPOSE 8501

CMD ["snakemake", "--cores", "1"]
