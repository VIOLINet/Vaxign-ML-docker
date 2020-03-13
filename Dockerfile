FROM brinkmanlab/psortb_commandline:1.0.2

#############
# Vaxign-ML #
#############

MAINTAINER Edison Ong <edong@umich.edu>

RUN mkdir -p /tmp/results && chmod 777 /tmp/results

WORKDIR /app

COPY . /app

WORKDIR /app/lib/spaan/SPAAN
RUN gcc -o standard.o standard.c -lm -w
RUN gcc -o filter.o filter.c -lm -w
RUN gcc -o annotate.o annotate.c -lm -w
RUN gcc -o AAcompo/AAcompo.o AAcompo/AAcompo.c -lm -w
RUN gcc -o AAcompo/recognize.o AAcompo/recognize.c -lm -w
RUN gcc -o charge/charge.o charge/charge.c -lm -w
RUN gcc -o charge/recognize.o charge/recognize.c -lm -w
RUN gcc -o hdr/hdr.o hdr/hdr.c -lm -w
RUN gcc -o hdr/recognize.o hdr/recognize.c -lm -w
RUN gcc -o multiplets/multiplets.o multiplets/multiplets.c -lm -w
RUN gcc -o multiplets/recognize.o multiplets/recognize.c -lm -w
RUN gcc -o dipep/dipep.o dipep/dipep.c -lm -w
RUN gcc -o dipep/recognize.o dipep/recognize.c -lm -w
RUN gcc -o finalp1.o finalp1.c -lm -w

WORKDIR /app
RUN apt-get update
RUN apt-get install -y software-properties-common
RUN add-apt-repository ppa:deadsnakes/ppa
RUN apt-get update
RUN apt-get install -y python3.6 python3-pip
RUN python3.6 -m pip install --trusted-host numpy==1.14.2 scipy==1.2.1 scikit-learn==0.20.3 xgboost==0.81 biopython==1.72 matplotlib==2.2.2 pandas==0.20.3 pathlib

ENTRYPOINT ["/usr/bin/env"]
CMD ["python3.6", "VaxignML.py"]
