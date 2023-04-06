# syntax=docker/dockerfile:1

FROM python:3.6.1

WORKDIR /walker

COPY requirements.txt .
RUN pip install --upgrade pip
RUN pip install -r requirements.txt
COPY . .

ENTRYPOINT [ "python", "walker.py"]

