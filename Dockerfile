FROM davidliwei/seurat-mageck:0.5.9.2

RUN apt-get update  


COPY . /app

WORKDIR /app

RUN R --vanilla -e 'library(reticulate);use_python("/usr/local/bin/python");py_config()'
RUN R CMD build bioconductor & R CMD INSTALL scMAGeCK*.gz

ENTRYPOINT echo "Welcome to scMAGeCK Docker" &  /bin/bash 

