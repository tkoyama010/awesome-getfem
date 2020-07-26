FROM getfemdoc/getfem:stable
ENV DEBIAN_FRONTEND noninteractive
RUN apt update && apt -y install python3-pip

# install the notebook package
RUN apt -y install libgl1-mesa-dev
RUN apt -y install xvfb
RUN apt -y install git
RUN pip3 install --no-cache --upgrade pip && \
    pip3 install --no-cache jupyterlab && \
    pip3 install --no-cache git+git://github.com/pyvista/pyvista@master

# create user with a home directory
ARG NB_USER
ARG NB_UID
ENV USER ${NB_USER}
ENV HOME /home/${NB_USER}
ENV PYVISTA_OFF_SCREEN true
ENV PYVISTA_USE_PANEL true
ENV PYVISTA_PLOT_THEME document
# This is needed for Panel - use with cuation!
ENV PYVISTA_AUTO_CLOSE false

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