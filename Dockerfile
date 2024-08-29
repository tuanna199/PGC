FROM ubuntu:20.04

# Install system dependencies

ARG DEBIAN_FRONTEND=noninteractive
ENV DEBIAN_FRONTEND=${DEBIAN_FRONTEND}

RUN apt-get update && apt-get install -y \
    software-properties-common \
    build-essential \
    libssl-dev \
    libffi-dev \
    libncurses5-dev \
    libncursesw5-dev \
    libreadline-dev \
    libsqlite3-dev \
    libgdbm-dev \
    libc6-dev \
    zlib1g-dev \
    libpython3-dev \
    libffi-dev \
    git \
    wget \
    unzip \
    r-base

# Install Miniconda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh \
    && /bin/bash /tmp/miniconda.sh -b -p /opt/conda \
    && rm /tmp/miniconda.sh
ENV PATH="/opt/conda/bin:$PATH"

# Create a Conda environment and install Python packages
COPY environment.yml .
RUN conda env create -f environment.yml

# Make RUN commands use the new environment:
SHELL ["conda", "run", "-n", "docker_pgc", "/bin/bash", "-c"]

RUN R -e "install.packages(c('argparser', 'r-base', 'dplyr', 'futile.logger', 'ggplot2', 'ggpubr', 'gridExtra', 'qqman', 'readr', 'reshape2', 'scales', 'scrime', 'VennDiagram'), repos='https://cloud.r-project.org/')"

# Clone git repo
RUN git clone https://github.com/DReichLab/EIG.git /app/EIG

# Download, unzip, and remove the ZIP file
RUN wget https://data.broadinstitute.org/igv/projects/downloads/2.18/IGV_Linux_2.18.2_WithJava.zip -O /tmp/IGV.zip \
    && unzip /tmp/IGV.zip -d /app/ \
    && rm /tmp/IGV.zip

# Copy your Snakemake pipeline files
COPY . /app
WORKDIR /app

# Activate the environment in the container entry point
CMD ["conda", "run", "-n", "docker_pgc", "/bin/bash"]
# CMD ["/bin/bash"]

