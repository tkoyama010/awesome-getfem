FROM getfemdoc/getfem:stable
ENV DEBIAN_FRONTEND noninteractive
RUN apt update && apt -y install python3-pip

# Xvfb
RUN apt-get install -yq --no-install-recommends \
    xvfb \
    x11-utils \
    libx11-dev \
    qt5-default \
    mayavi2 \
    && apt-get clean

ENV DISPLAY=:99

# install the notebook package
RUN pip3 install --no-cache --upgrade pip && \
    pip3 install --no-cache jupyterlab

# create user with a home directory
ARG NB_USER
ARG NB_UID
ENV USER ${NB_USER}
ENV HOME /home/${NB_USER}

RUN adduser --disabled-password \
    --gecos "Default user" \
    --uid ${NB_UID} \
    ${NB_USER}
WORKDIR ${HOME}
USER root
COPY . ${HOME}
RUN pip3 install -r requirements.txt
RUN chown -R ${NB_USER} ${HOME}
USER ${USER}

# Install Python packages
RUN pip3 install vtk && \
    pip3 install boto && \
    pip3 install h5py && \
    pip3 install nose && \
    pip3 install ipyevents && \
    pip3 install ipywidgets && \
    pip3 install nibabel && \
    pip3 install numpy && \
    pip3 install pillow && \
    pip3 install pyqt5 && \
    pip3 install scikit-learn && \
    pip3 install scipy && \
    pip3 install xvfbwrapper

# Add an x-server to the entrypoint. This is needed by Mayavi
ENTRYPOINT ["tini", "-g", "--", "xvfb-run"]
