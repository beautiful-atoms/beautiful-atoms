# Base image for the blender environment. Used for CD/CI test and debug
# Dockerfile adapted from https://github.com/nytimes/rd-blender-docker
# To build the Docker image, create the build from root directory:
# docker build --build-arg BLENDER_VERSION=3.0 . -f Dockerfiles/Dockerfile.base_3.0 -t <image_name>:blender3.0

ARG BLENDER_VERSION=3.1
FROM nytimes/blender:${BLENDER_VERSION}-gpu-ubuntu18.04

LABEL Author="T.Tian <alchem0x2a@gmail.com>,X.Wang <xingwang1991@gmail.com>"
LABEL Title="Blender base environment for testing beautiful-atoms"
LABEL RepoUrl="https://github.com/beautiful-atoms/beautiful-atoms"
LABEL ImageUrl="https://github.com/beautiful-atoms/beautiful-atoms/pkgs/container/beautiful-env"


RUN rm /bin/sh && ln -s /bin/bash /bin/sh


# Enviorment variables
ENV CONDA_DIR="/opt/conda"
ENV PATH="${CONDA_DIR}/bin:${PATH}"


# Install miniconda 4.11.0 (minimal support py 3.10)
ENV MINICONDA_VERSION=4.11.0 \
    MINICONDA_MD5=4e2f31e0b2598634c80daa12e4981647 \
    CONDA_VERSION=4.12.0 \
    PYTHON_VERSION=3.10.2


WORKDIR /tmp
RUN wget --quiet https://repo.continuum.io/miniconda/Miniconda3-py39_${MINICONDA_VERSION}-Linux-x86_64.sh && \
    echo "${MINICONDA_MD5} *Miniconda3-py39_${MINICONDA_VERSION}-Linux-x86_64.sh" | md5sum -c - && \
    /bin/bash Miniconda3-py39_${MINICONDA_VERSION}-Linux-x86_64.sh -f -b -p $CONDA_DIR && \
    rm Miniconda3-py39_${MINICONDA_VERSION}-Linux-x86_64.sh && \
    conda config --system --prepend channels conda-forge && \
    conda config --system --set show_channel_urls true && \
    conda install --quiet --yes \
    -c conda-forge \
    conda=${CONDA_VERSION} \
    python=$PYTHON_VERSION \
    pip &&\
    conda update --all --quiet --yes &&\
    conda clean --all -f -y

ADD ./ /tmp

RUN . ${CONDA_DIR}/bin/activate &&\
    yes | python install.py --generate-env-file ./env.yml ${BLENDER_PATH} &&\
    rm -rf ${BLENDER_PATH}/python &&\
    ln -s ${CONDA_DIR} ${BLENDER_PATH}/python


RUN conda env update -n base --file ./env.yml &&\
    conda clean --all -f -y &&\
    rm -rf /tmp/*


# Set the working directory
ENV B_USER="blender"
ENV B_UID="1000"
ENV B_GID="100"
ENV HOME=/home/${B_USER}


# Add user
RUN useradd -m -s /bin/bash -N -u $B_UID $B_USER && \
    mkdir -p /workdir && \
    chmod g+w /etc/passwd && \
    chown -R ${B_UID}:${B_GID} ${HOME} &&\
    chown -R ${B_UID}:${B_GID} /workdir


USER ${B_USER}
SHELL ["/bin/bash", "-c"]
WORKDIR /workdir
