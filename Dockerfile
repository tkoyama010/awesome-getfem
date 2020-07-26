FROM getfemdoc/getfem:stable
ENV DEBIAN_FRONTEND noninteractive
RUN apt update && apt -y install python3-pip

# install the notebook package
RUN pip3 install --no-cache --upgrade pip && \
    pip3 install --no-cache jupyterlab

# install pyvista
RUN apt-get -y install libgl1-mesa-dev
RUN apt-get -y install xvfb
RUN set -x
ENV DISPLAY :99.0
ENV PYVISTA_OFF_SCREEN true
ENV PYVISTA_USE_PANEL true
RUN which Xvfb
RUN Xvfb :99 -screen 0 1024x768x24 > /dev/null 2>&1 &
RUN sleep 3
RUN set +x
RUN exec "$@"
RUN pip3 install --no-cache pyvista

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
