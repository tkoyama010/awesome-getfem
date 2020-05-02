FROM getfemdoc/getfem:latest
ENV DEBIAN_FRONTEND noninteractive
RUN apt update && \
    apt -y install python3-pip && \
    apt -y install python3-meshio && \
    apt -y install python3-matplotlib

# install the notebook package
RUN pip3 install --no-cache --upgrade pip && \
    pip3 install --no-cache notebook

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
RUN chown -R ${NB_USER} ${HOME}
USER ${USER}
