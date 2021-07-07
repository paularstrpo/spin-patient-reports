FROM rocker/verse:4.0.3

# libraries & dependencies
RUN apt-get update
RUN Rscript -e "install.packages(c('BiocManager','optparse', 'knitr', 'rmarkdown', 'knitr', 'igraph', 'tidyverse', 'RColorBrewer', 'jsonlite', 'svglite', 'RJSONIO', 'grid', 'circlize', 'formattable', 'kableExtra','survival', 'survminer', 'ggraph', 'tidygraph'), dependencies=TRUE, repos = 'http://cran.us.r-project.org')"
RUN Rscript -e "BiocManager::install(c('ComplexHeatmap'))"

ADD assets/ ./assets/
ADD refData/ ./refData/
ADD functions.r .
ADD spin_report.r .
ADD report_template.rmd .

CMD ["/bin/bash"]