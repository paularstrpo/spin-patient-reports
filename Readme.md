# Report Visualization

## Required Dependencies

```
Rscript -e "install.packages(c('BiocManager','optparse', 'knitr', 'rmarkdown', 'knitr', 'igraph', 'tidyverse', 'RColorBrewer', 'jsonlite', 'svglite', 'RJSONIO', 'grid', 'circlize', 'formattable', 'kableExtra','survival', 'survminer', 'ggraph', 'tidygraph'), dependencies=TRUE, repos = 'http://cran.us.r-project.org')"
Rscript -e Rscript -e "BiocManager::install(c('ComplexHeatmap'))"
```

## General Usage:

### Example Command

```
#!/bin/bash

docker build -t daphni_report . # run with docker
docker run -v "/data1/users/restrp01/example_data:/exampleInputs" -v "/data2/mm/data/MSSM/sema4_clinical/21064632/processed:/processed" daphni_report \
Rscript spin_report.r \
 --tiers "1A,1B" \
 --sampleName "21064632" \
 --sexChrs "XX" \
 --dateSampleCollected "05/06/2020" \
 --dateSampleSequenced "04/03/2021" \
 --patientName "JANE DOE" \
 --MRN "XXXXXX" \
 --DOB "01/01/1960" \
 --IgClass "N/A" \
 --LightChainType "N/A" \
 --RevisedISS "I" \
 --DiseaseStatus "Unknown/TBD" \
 --drugDetailsRanked "/processed/prediction_engine/all.drug_details.ranked.prediction_engine.results.tsv" \
 --rxBucketSummary "/processed/prediction_engine/tier1.rxbucket_summary.ranked.prediction_engine.results.tsv" \
 --variantSupportSummary "/processed/prediction_engine/variant_summary.ranked.prediction_engine.results.tsv" \
 --prognosticMarkerSummary "/processed/prediction_engine/prognostic_marker.ranked.prediction_engine.results.tsv" \
 --mmPSNFile "/processed/secondary/Predicted_class.csv" \
 --tumorPurityPloidy "/processed/primary/tumor.facets_output.txt" \
 --facetsCNCFResultsFile "/processed/primary/tumor.facets_cncf.txt" \
 --treeFile "/processed/secondary/sample.summ.json.gz" \
 --zScoreTable "/exampleInputs/zscores_latest.csv" \
 --gep70File "/exampleInputs/gep70scores_latest.tsv" \
 --seliFile '/processed/secondary/selinescores_latest.tsv' \
 --scarFile "/exampleInputs/scarScore_latest.tsv" \
 --tmbFile "/exampleInputs/tumorMutationBurden_latest.tsv" \
 --geneCoordinates "refData/hg38_geneCoordinates_TxDB.tsv" \
 --psnReferenceTable "refData/mmrf_mmpsn_subgroup_survival.tsv"  \
 --cytobandCoordinates 'refData/hg38_cytoband_coordinates.tsv'
```

### Options

```
Usage: spin_report.r [options]
----------------------
[DAPHNI v2.0] Precision Medicine Report Generation
----------------------

Options:
        --tiers=TIERS
                comma-separated list, enclosed in quotes, of main tiers to include in the summary section of the report [REQUIRED]
        --patientName=PATIENTNAME
                patient name, in the format of first last, must be enclosed in quotes
        --sexChrs=SEXCHRS
                Gender of patient, specified as either XX or XY [REQUIRED]
        --MRN=MRN
                MRN of patient, must be enclosed in quotes
        --DOB=DOB
                patient date of birth, in format of mm/dd/yyyy, must be enclosed in quotes
        --IgClass=IGCLASS
        --LightChainType=LIGHTCHAINTYPE
        --RevisedISS=REVISEDISS
        --DiseaseStatus=DISEASESTATUS
        --sampleName=SAMPLENAME
                sample ID name for patient to run prediction engine on. [REQUIRED]
        --dateSampleCollected=DATESAMPLECOLLECTED
                date of sample collection, in format of mm/dd/yyyy, must be enclosed in quotes
        --dateSampleSequenced=DATESAMPLESEQUENCED
                date that sample sequencing was performed, in format of mm/dd/yyyy, must be enclosed in quotes
        --drugDetailsRanked=DRUGDETAILSRANKED
                prediction engine full drug details table (tab-delimited; e.g. all.drug_details.ranked.prediction_engine.results.tsv)
        --rxBucketSummary=RXBUCKETSUMMARY
                prediction engine tier 1 rx bucket summary table (tab-delimited; e.g. tier1.rxbucket_summary.ranked.prediction_engine.results.tsv)
        --variantSupportSummary=VARIANTSUPPORTSUMMARY
                prediction engine variant support summary table (tab-delimited; e.g. variant_summary.ranked.prediction_engine.results.tsv)
        --prognosticMarkerSummary=PROGNOSTICMARKERSUMMARY
                prediction engine prognostic marker results table (tab-delimited; e.g. prognostic_markers.prediction_engine.results.tsv)
        --zScoreTable=ZSCORETABLE
                Table with expression zscores for all patients, including the patient of interest (tab-delimited; e.g. zscores_latest.tsv)
        --gep70File=GEP70FILE
                Table containing gep70 scores at cohort level, including the patient of interest (tab-delimited; e.g. gep70scores_latest.tsv)
        --seliFile=SELIFILE
                Table containing selinexor signature data at cohort level, including the patient of interest (tab-delimited; e.g. selinescores_latest.tsv)
        --msiFile=MSIFILE
                Table containing MSI scores at cohort level, including the patient of interest (tab-delimited; e.g. msiScore_latest.tsv)
        --scarFile=SCARFILE
                Table containing scar scores at cohort level, including the patient of interest (tab-delimited; e.g. scarScore_latest.tsv)
        --tmbFile=TMBFILE
                Table containing tumor mutation burdens at cohort level, including the patient of interest (tab-delimited; e.g. tumorMutationBurden_latest.tsv)
        --somFile=SOMFILE
                prediction engine results for somatic mutations (tab-delimited; e.g. somatic_mutation.prediction_engine.results.tsv)
        --cnaFile=CNAFILE
                prediction engine results for CNVs (tab-delimited; e.g. cna.prediction_engine.results.tsv)
        --exprFile=EXPRFILE
                prediction engine results for expression (tab-delimited; e.g. expression.prediction_engine.results.tsv)
        --mmPSNFile=MMPSNFILE
                file containing predicted mm-psn subgroup information for the sample (comma-delimited; e.g. Predicted_class.csv
        --treeFile=TREEFILE
                gzip-compressed json file with clonal tree structure (e.g. sample.sum.json.gz)
        --facetsCNCFResultsFile=FACETSCNCFRESULTSFILE
                facets cncf file (tab-delimited; e.g. .tumor.facets_cncf.txt)
        --tumorPurityPloidy=TUMORPURITYPLOIDY
                facets output file containing purity and ploidy (tab-delimited; e.g. .tumor.facets_output.txt)
        --cytobandCoordinates=CYTOBANDCOORDINATES
                Table with full cytoband coordinates (tab-delimited; e.g. hg38_cytoband_coordinates.tsv) [REQUIRED]
        --psnReferenceTable=PSNREFERENCETABLE
                Table with MMRF PSN survival per subgroup information (tab-delimited; e.g. mmrf_mmpsn_subgroup_survival.tsv) [REQUIRED]
        --geneCoordinates=GENECOORDINATES
                Table with genomic start and end coordinates for each gene (tab-delimited; e.g hg38_geneCoordinates_TxDB.tsv) [REQUIRED]
        -h, --help
                Show this help message and exit
```

### Outputs:

- `*.main.daphni_report.html` the main report output
- `*.kable.daphni_report` formatted html tables that are included in the report
- `*.plot.dahni_report.svg` plots that are included in the report
