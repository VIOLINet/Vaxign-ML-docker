FROM brinkmanlab/psortb_commandline:1.0.2

#############
# Vaxign-ML #
#############

MAINTAINER Edison Ong <edong@umich.edu>

RUN mkdir -p /tmp/results && chmod 777 /tmp/results

WORKDIR /app

COPY . /app

RUN apt-get update
RUN apt-get install -y software-properties-common
RUN add-apt-repository ppa:jonathonf/python-3.6
RUN apt-get update
RUN apt-get install -y python3.6 python3-pip
RUN python3.6 -m pip install --trusted-host numpy==1.14.2 scipy==1.2.1 scikit-learn==0.20.3 xgboost==0.81 biopython==1.72 matplotlib==2.2.2 pandas==0.20.3 pathlib

ENTRYPOINT ["/usr/bin/env"]
CMD ["python3.6", "VaxignML.py"]
