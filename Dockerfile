FROM getfemdoc/getfem:stable
ENV DEBIAN_FRONTEND noninteractive
RUN apt update
RUN apt-get -y -f install
RUN apt -y install python3-pip

# install the notebook package
RUN apt -y install git
RUN apt-get -y install libgl1-mesa-dev
RUN pip3 install PyQt5
RUN apt-get -y install xvfb
RUN pip3 install --no-cache --upgrade pip && \
    pip3 install --no-cache jupyterlab && \
    pip3 install --no-cache pyviz && \
    pip3 install --no-cache bokeh && \
    pip3 install --no-cache pyviz_comms && \
    pip3 install --no-cache panel && \
    pip3 install --no-cache lxml && \
    pip3 install --no-cache matplotlib && \
    pip3 install --no-cache pyct && \
    pip3 install --no-cache git+git://github.com/tkoyama010/pyvista@patch-3

# create user with a home directory
ARG NB_USER
ARG NB_UID
ENV USER ${NB_USER}
ENV HOME /home/${NB_USER}
ENV PYVISTA_VIRTUAL_DISPLAY true
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
