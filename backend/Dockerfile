# Use an official Python runtime as a parent image
FROM python:3.10-slim

# Set the working directory in the container
WORKDIR /app

# Install system dependencies
# build-essential for compiling some python packages
# git for installing packages from git
RUN apt-get update && apt-get install -y \
    build-essential \
    git \
    default-jre \
    autodock-vina \
    && rm -rf /var/lib/apt/lists/*

# Copy the requirements file into the container
COPY requirements.api.txt .

# Install any needed packages specified in requirements.txt
RUN pip install --no-cache-dir -r requirements.api.txt

# Copy the rest of the application code
COPY . .

# Expose the port that Uvicorn will run on
# Cloud Run expects the container to listen on $PORT (default 8080)
ENV PORT=8080

# Run api.py when the container launches
CMD ["uvicorn", "api:app", "--host", "0.0.0.0", "--port", "8080"]
