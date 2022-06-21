FROM r-base:4.0.2

LABEL maintainer="LAMOUCHE Jean-Baptiste <jbaptiste.lamouche@gmail.com>" \
    Software="DPNI" \
    Version="0.1" \
    Description="DPNI"

###########################
## Miniconda & Snakemake  #
###########################

RUN apt-get update -y
RUN cd /tmp
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-py39_4.12.0-Linux-x86_64.sh && sh Miniconda3-py39_4.12.0-Linux-x86_64.sh -b -p /usr/local/lib/miniconda3
ENV PATH="/usr/local/lib/miniconda3/bin:$PATH"
RUN conda install -y -c bioconda -c conda-forge snakemake samtools
COPY . /app
RUN conda env create -f /app/environment.yml
RUN echo "source activate dpni" >> ~/.bashrc
RUN chmod -R +x /app
ENV PATH="/app/bin:$PATH"
WORKDIR /app

ENTRYPOINT ["/bin/bash"]