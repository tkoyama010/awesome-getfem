FROM getfemdoc/getfem:latest

RUN apt-get update
RUN apt-get install -y sudo
RUN apt-get -y install tzdata
ENV TZ=Asia/Tokyo
RUN apt-get -y install python3-pip
RUN apt-get -y install nodejs npm
RUN apt-get -y install language-pack-ja
RUN sudo update-locale LANG=ja_JP.UTF-8
RUN pip3 install jupyterlab
RUN apt-get -y install python3-getfem++
RUN apt-get -y install xvfb
RUN pip3 install pyvista
RUN pip3 install pyvirtualdisplay
ENV PYVISTA_OFF_SCREEN true
ENV PYVISTA_USE_PANEL true
ENV PYVISTA_PLOT_THEME document
ENV PYVISTA_AUTO_CLOSE false

