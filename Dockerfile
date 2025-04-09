# Use an official Python image as base
FROM python:3.9-slim

# Image metadata
LABEL maintainer="Chromatinsight Support <support@example.com>"
LABEL description="Docker image for Chromatinsight_py3, a tool for ChIP-seq data analysis"
LABEL version="1.0"

# Prevent Python from writing .pyc files
ENV PYTHONDONTWRITEBYTECODE=1
# Prevent Python from buffering stdout/stderr
ENV PYTHONUNBUFFERED=1

# Create a non-root directory for the application
WORKDIR /app

# Install dependencies
RUN pip install --no-cache-dir pandas scikit-learn numpy matplotlib

# Copy the source code
COPY . /app/

# Create a non-root user to run the application
RUN groupadd -r chromatinuser && \
    useradd -r -g chromatinuser -d /home/chromatinuser -s /sbin/nologin -c "Chromatinsight User" chromatinuser && \
    chown -R chromatinuser:chromatinuser /app

# Switch to non-root user
USER chromatinuser

# Default command to run a unit test as an example
CMD ["python", "chromatinsight_utest.py"]
