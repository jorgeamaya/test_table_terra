# Use a pinned, stable Miniconda image instead of "latest"
FROM continuumio/miniconda3:25.3.1-1

# Set a working directory
WORKDIR /app

# Copy the environment specification into the image
COPY . /app

# Create the Conda environment
RUN conda env create -f /app/environment.yml

# Install your package inside the environment
RUN /opt/conda/bin/conda run -n protbindscreen_env pip install /app

# Make this environment the default Python
ENV PATH=/opt/conda/envs/protbindscreen_env/bin:$PATH

# Default shell
CMD ["/bin/bash"]
