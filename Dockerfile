FROM ubuntu:22.04

ARG DEBIAN_FRONTEND=noninteractive

COPY . /pisces
WORKDIR /pisces

# Install dependencies for R key
RUN apt-get update && apt-get install -y gnupg software-properties-common
# Install R
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu jammy-cran40/'
RUN apt update
RUN apt install -y r-base

# Install PISCES
RUN apt-get install -y python3 python3-pip git && ln -s /usr/bin/python3 /usr/bin/python
RUN pip install -e .
RUN pisces_dependencies

ENTRYPOINT ["/usr/local/bin/pisces"]