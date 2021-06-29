# Report Visualization

## Get Started

Install required dependencies
```
libs <- c('tidyverse', 'RColorBrewer', 'ComplexHeatmap', 
          'circlize', 'kableExtra', 'formattable', 'fontawesome',
          'grid', 'gridExtra', 'gridBase', 'gridGraphics',
          'svglite', 'knitr', 'rmarkdown')

install.packages(libs)
```

General Usage:

```
docker run -v "$PWD:/data" test_report \
Rscript spinReportFromRmd.R \
 --sampleName "MM-5130-DNA-T-03" \
 --sexChrs "XX" \
 --dateSampleCollected "05/06/2020" \
 --dateSampleSequenced "04/03/2021" \
 --patientName "JANE DOE" \
 --MRN "XXXXXXXXX" \
 --DOB "04/20/1969" \
 --IgClass "N/A" \
 --LightChainType "N/A" \
 --RevisedISS "I" \
 --DiseaseStatus "Unknown/TBD" \
 --drugDetailsRanked "exampleInputs/MM_5130_T_03/drug_details.ranked.prediction_engine.results.tsv" \
 --rxBucketSummary "exampleInputs/MM_5130_T_03/rxbucket_summary.ranked.prediction_engine.results.tsv" \
 --variantSupportSummary "exampleInputs/MM_5130_T_03/variant_summary.ranked.prediction_engine.results.tsv" \
 --treeFile "exampleInputs/MM_5130_T_03/sample.summ.json.gz" \
 --zScoreTable "exampleInputs/zscores_mssm.tsv" \
 --cnaFile "exampleInputs/MM_5130_T_03/cna.prediction_engine.results.tsv" \
 --somFile "exampleInputs/MM_5130_T_03/somatic_mutation.prediction_engine.results.tsv" \
 --snpFile "exampleInputs/MM_5130_T_03/germline_snp.prediction_engine.results.tsv" \
 --exprFile "exampleInputs/MM_5130_T_03/expression.prediction_engine.results.tsv" \
 --facetsCNCFResultsFile "exampleInputs/MM_5130_T_03/tumor.facets_cncf.txt" \
 --mmPSNPlotListData "refData/mmPSN_survivalCurvePlotList.rds"  \
 --scarFile "exampleInputs/scarScore_latest.tsv" \
 --tmbFile "exampleInputs/tumorMutationBurden_latest.tsv" \
 --gep70File "exampleInputs/gep70scores_latest.tsv" \
 --tumorPurityPloidy "exampleInputs/MM_5130_T_03/tumor.facets_output.txt" \
 --geneCoordinateBED "refData/gene_conversion_table_Grch38_with_bands.txt"
```