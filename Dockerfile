# Start from the official Python base image
FROM python:3.11-slim

# Set working directory
WORKDIR /.

# Default command to run Python shell
CMD ["python"]
