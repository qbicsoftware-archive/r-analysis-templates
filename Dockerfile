FROM nfcore/base
LABEL authors="gisela.gabernet@qbic.uni-tuebingen.de" \
      description="Docker image containing all requirements for the downstream RNAseq analysis"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/qbicsoftware-r-analysis-templates/bin:$PATH
