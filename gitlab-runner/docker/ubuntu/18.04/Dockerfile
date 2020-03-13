# Use an official Debian runtime as a base image
FROM ubuntu:18.04
MAINTAINER Yann Pouillon "yann.pouillon@materialsevolution.es"

# Set the working directory to /app
WORKDIR /siesta

# Copy the current directory contents into the container at /app
# ADD . /siesta #This is a bad idea, every single change in the configuration will require to compute entire image

# Prepare the environment
ENV DEBIAN_FRONTEND=noninteractive

# Install the required development packages
COPY packages.txt ./
RUN apt update -y && \
    apt install apt-utils -y && \
    apt upgrade -y && \
    apt install locales -y && \
    /bin/bash -c "xargs -a <(awk '/^\s*[^#]/' packages.txt) -r -- apt install -y" && \
    rm -r /var/cache/apt

RUN sed -i -e 's/# en_US.UTF-8 UTF-8/en_US.UTF-8 UTF-8/' /etc/locale.gen && \
    dpkg-reconfigure --frontend=noninteractive locales && \
    update-locale LANG=en_US.UTF-8

ENV LANG en_US.UTF-8
ENV LC_ALL en_US.UTF-8

# Install any needed packages specified in requirements.txt
COPY requirements.txt ./
RUN pip3 install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir --upgrade wheel && \
    pip install --no-cache-dir -r requirements.txt && \
    pip install --no-cache-dir --upgrade git+https://gitlab.com/ase/ase.git@master

# Build the currently available versions of SIESTA
COPY install-siesta-4.0.2.sh ./
COPY arch.make-siesta-4.0.2-gnu-mpi ./
COPY arch.make-siesta-4.0.2-gnu-serial ./
RUN git clone https://gitlab.com/siesta-project/siesta.git /siesta/trunk
RUN /bin/bash /siesta/install-siesta-4.0.2.sh

# Define container name
ENV NAME siesta-runner

# Open a terminal when the container launches
CMD ["/bin/bash", "-login"]
