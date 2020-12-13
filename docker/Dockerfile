FROM ubuntu
RUN apt-get update
RUN apt-get -y install tzdata
RUN apt-get -y install python3-getfem++
RUN apt-get -y install python3-pip
RUN apt-get -y install libgl1-mesa-dev
RUN apt-get -y install nodejs npm
RUN apt-get -y install xvfb
RUN pip3 install jupyterlab
RUN pip3 install pyvista
RUN pip3 install pyvirtualdisplay
ENV PYVISTA_OFF_SCREEN true
ENV PYVISTA_USE_PANEL true
ENV PYVISTA_PLOT_THEME document
ENV PYVISTA_AUTO_CLOSE false
