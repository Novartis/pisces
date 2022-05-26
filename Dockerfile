FROM ubuntu:18.04

COPY . /pisces
WORKDIR /pisces

# Install R
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran40/'
RUN apt update
RUN apt install r-base

# Install PISCES
RUN python setup.py install
RUN pisces dependencies

ENTRYPOINT ["/usr/local/bin/pisces"]