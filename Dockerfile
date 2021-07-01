FROM rocker/verse:4.0.3

# libraries & dependencies
RUN apt-get update
RUN Rscript -e "install.packages(c('BiocManager','optparse', 'knitr', 'rmarkdown', 'igraph', 'tidyverse', 'RColorBrewer', 'jsonlite', 'svglite', 'RJSONIO', 'grid', 'circlize', 'formattable', 'kableExtra','survival', 'survminer', 'ggraph', 'tidygraph'), dependencies=TRUE, repos = 'http://cran.us.r-project.org')"
RUN Rscript -e "BiocManager::install(c('GenVisR', 'ComplexHeatmap'))"

# copy rendering script
COPY spinReportFromRmd.R /bin/
RUN chmod +x /bin/spinReportFromRmd.R

# copy markdown source code in (required for rendering!!)
COPY spinPatientReport.Rmd /bin/
# logos
COPY MS_RGB_Hrztl.png /bin/
COPY MS_KO_Hrztl.png /bin/
CMD ["/bin/bash"]