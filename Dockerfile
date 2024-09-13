# Use an official Python runtime as a parent image
FROM python:3.12-slim

# Set the working directory in the container
WORKDIR /usr/src/app

# Install system dependencies
RUN apt-get update && apt-get install -y \
    wget \
    git \
    tabix \
    bgzip \
    && rm -rf /var/lib/apt/lists/*

# Install pipenv
RUN pip install pipenv

# Copy the current directory contents into the container at /usr/src/app
COPY . /usr/src/app

# Install Python dependencies
RUN pipenv install --deploy --ignore-pipfile

# Set up seqrepo
RUN mkdir -p /usr/local/share/seqrepo \
    && chown -R $USER /usr/local/share/seqrepo
RUN pipenv run seqrepo init -i auto-acmg
RUN pipenv run seqrepo fetch-load -i auto-acmg -n RefSeq <sequence identifiers>

# Make port 8000 available to the world outside this container
EXPOSE 8000

# Define environment variable
ENV NAME World

# Run uvicorn when the container launches
CMD ["pipenv", "run", "uvicorn", "app.main:app", "--host", "0.0.0.0", "--port", "8000"]
