# golob/ya16sdb:0.2C

FROM      ubuntu:18.04

env DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y \
python3-dev \
python3-pip \
postgresql \
wget \
git \
&& apt-get clean \
&& apt-get purge \
&& rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN pip3 install pip --upgrade
RUN pip3 install \
        awscli>=1.15.14 \
        boto3>=1.7.14 \
        numpy>=1.14.2 \
        bucket_command_wrapper==0.3.1 \
        biopython>=1.68 \
        Cython \
        taxtastic

RUN mkdir -p /fh && mkdir -p /app && mkdir -p /src
RUN mkdir -p /mnt/inputs/file && mkdir -p /mnt/outputs/file && mkdir /scratch
RUN mkdir /logs && mkdir /records

RUN cd /src && \
wget https://github.com/torognes/vsearch/releases/download/v2.13.6/vsearch-2.13.6-linux-x86_64.tar.gz && \
tar xzvf vsearch-2.13.6-linux-x86_64.tar.gz && \
cp /src/vsearch-2.13.6-linux-x86_64/bin/vsearch /usr/local/bin/ && \
rm -r /src/*

RUN cd /src && \
wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.9.0/ncbi-blast-2.9.0+-x64-linux.tar.gz && \
tar xzvf ncbi-blast-2.9.0+-x64-linux.tar.gz && \
cp  /src/ncbi-blast-2.9.0+/bin/makeblastdb /usr/local/bin && \
rm -r /src/*

ADD requirements.txt /src/
RUN pip3 install -r /src/requirements.txt && rm -r /src/*

RUN mkdir -p /db/
ADD data/rdp_16s_type_strains.fasta.bz2 /db/
RUN cd /db/ && bunzip2 rdp_16s_type_strains.fasta.bz2 -c > gzip > rdp_16s_type_strains.fasta.gz

ADD bin/* /usr/local/bin/
RUN chmod +x /usr/local/bin/



WORKDIR /working