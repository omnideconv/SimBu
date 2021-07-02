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
RUN R -e "install.packages(c('ggplot2','data.table','Seurat','Matrix','ggpubr'))"
RUN R -e "BiocManager::install(c('monocle'))"
RUN R -e "install.packages(c('hdf5r'))"



#docker run --rm -d -p 127.0.0.1:9999:8787 -v /home/alex/Documents/Studium_Bioinformatik/Master/Masterarbeit/Data:/home/Data -v /home/alex/Documents/Studium_Bioinformatik/Master/Masterarbeit/simulator:/home/simulator -e ROOT=TRUE -e PASSWORD=123456 --name rstudio9999 r_simulator
