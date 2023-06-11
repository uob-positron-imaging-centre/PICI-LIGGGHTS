FROM --platform=linux/amd64 dealii/dealii:v9.4.0-focal

# Set user to what Binder expects
ENV NB_USER dealii
ENV NB_UID 1000
ENV USER ${NB_USER}
ENV HOME /home/${NB_USER}

# Make sure the contents of our repo are in ${HOME}
COPY ./install_liggghts.sh ./requirements.txt ${HOME}
USER root
RUN chown -R ${NB_UID} ${HOME}
USER ${NB_USER}

# Python binaries will be installed here
ENV PATH="${PATH}:${HOME}/.local/bin"

# Install the notebook packages
RUN sudo apt update && sudo apt install python3-pip -y && sudo rm -rf /var/lib/apt/lists/*
RUN python3 -m pip install --no-cache-dir -r requirements.txt

# Jupyter needs the /usr/etc folder
USER root
RUN mkdir /usr/etc

# Install VTK and PICI-LIGGGHTS
RUN bash install_liggghts.sh

# LIGGGHTS Python interface
RUN mkdir /usr/etc/PICI-LIGGGHTS-Python && cp PICI-LIGGGHTS/python/liggghts.py /usr/etc/PICI-LIGGGHTS-Python
ENV PYTHONPATH="${PYTHONPATH}:/usr/etc/PICI-LIGGGHTS-Python"

# Cleanup
RUN rm -rf PICI-LIGGGHTS VTK-9.2.6 install_liggghts.sh requirements.txt

# Set user
USER ${NB_USER}
