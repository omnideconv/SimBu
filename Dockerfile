FROM rocker/verse:4.1.0
RUN apt-get update && apt-get install -y \
    libv8-dev \
    libudunits2-dev \
    liblzma-dev \
    libbz2-dev \
    git-all \
    python3-dev \
    python3-pip \
    libhdf5-dev


RUN R -e "install.packages(c('BiocManager'))"
RUN R -e "BiocManager::install(c('GEOquery', 'quantiseqr'))"
#RUN R -e "install.packages(c('Seurat'))"
RUN R -e "BiocManager::install(c('monocle'))"
RUN R -e "install.packages(c('hdf5r','ggplot2','data.table','Matrix','ggpubr','Hmisc','R.utils'))"
RUN R -e "install.packages(c('Seurat','Rcpp','remotes'))"
RUN R -e "remotes::install_github('mojaveazure/seurat-disk')"
RUN R -e "install.packages('anndata')"
RUN R -e "reticulate::install_miniconda()"
RUN R -e "anndata::install_anndata()"
RUN R -e "install.packages(c('tidyr'))"
RUN R -e "BiocManager::install('Biobase')"

RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir sfaira

CMD ["/bin/bash"]

#docker run --rm -v /home/alex/Documents/Studium_Bioinformatik/Master/Masterarbeit/Data:/home/Data -v /home/alex/Documents/Studium_Bioinformatik/Master/Masterarbeit/simulator:/home/simulator --name r_simulator r_simulator
