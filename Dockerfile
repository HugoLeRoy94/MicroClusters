# Use an official Ubuntu as the base image
FROM ubuntu:22.04

# Set the working directory in the container
WORKDIR /app

# Install necessary packages
RUN apt-get update && apt-get install -y \
    build-essential \
    python3 \
    python3-pip \
    cmake \
    && rm -rf /var/lib/apt/lists/*

# Copy the requirements file and install Python dependencies
COPY requirements.txt .
RUN pip3 install --no-cache-dir -r requirements.txt

# Copy the source code into the container
COPY src/ ./src/

# Build the C++ application
WORKDIR /app/src
RUN make

# Set the command to run your application
CMD ["./build/your_cpp_executable"]
