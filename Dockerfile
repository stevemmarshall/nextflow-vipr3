FROM conda/miniconda2

MAINTAINER Thanh Le Viet lethanhx2k@gmail.com

RUN apt update &&\
    apt install -y wget &&\
    conda config --add channels defaults &&\
    conda config --add channels conda-forge &&\
    conda config --add channels r &&\
    conda config --add channels bioconda

#Install common bio tools
RUN conda install bwa=0.7.15 \
                  spades=3.10.1\
                  samtools=1.3 \
                  lofreq=2.1.2 \
                  bedtools=2.25.0 \
                  mummer=3.23 \
                  fastp \
                  bbmap=37.17 \
                  matplotlib \
    && conda clean --all

#Install custom scripts
RUN mkdir /scripts && \
    cd /scripts &&\
    wget https://raw.githubusercontent.com/andreas-wilm/simple-contig-joiner/master/simple_contig_joiner.py &&\
    chmod +x simple_contig_joiner.py &&\
    wget https://raw.githubusercontent.com/gis-rpd/pipelines/master/germs/vipr/aux/plot.py &&\
    sed -i 's/python3/python/g' plot.py &&\
    chmod +x plot.py


ENV PATH /scripts:$PATH

# Add user biodocker with password biodocker
# https://github.com/BioContainers/containers/blob/master/biodocker/Dockerfile
RUN groupadd fuse \
    && useradd --create-home --shell /bin/bash --user-group --uid 1000 --groups sudo,fuse biodocker && \
    echo `echo "biodocker\nbiodocker\n" | passwd biodocker`

# Change user
USER biodocker

CMD ["/bin/bash"]
