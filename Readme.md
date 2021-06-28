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

General Usage

```
# Patient 001
Rscript spinReportFromRmd.R \
 --reportTemplateFile "refData/spinPatientReport.Rmd" \
 --sampleName "ISMMS01" --sexChrs "XX" --dateSampleCollected "01-01-2020" \
 --dateSampleSequenced "02-02-2020" --patientName="Jane Doe" --MRN '1234567' \
 --DOB "04-20-1969" --IgClass "IgG" --LightChainType "Kappa" --RevisedISS "VI" \
 --DiseaseStatus "Relapsed" --mainCloneResults "exampleInputs/ISMMS01/mainCloneResults.tsv" \
 --treeFile "exampleInputs/ISMMS01/sample.summ.json.gz" --zScoreTable "exampleInputs/zscores_latest.csv" \
 --sampleNameList "exampleInputs/sampleNameList.csv" --msiFileList "exampleInputs/msiFileList.csv" \
 --tmbFileList "exampleInputs/tmbFileList.csv" --mmPSNFile "exampleInputs/ISMMS01/Predicted_class.csv" \
 --snvFile "exampleInputs/ISMMS01/snvMainCloneDetailsTable.tsv" --exprFile "exampleInputs/ISMMS01/expDetailsTable.tsv" \
 --immunoFile "exampleInputs/ISMMS01/immunoResults.tsv" --facetsCNCFResultsFile "exampleInputs/ISMMS01/tumor.facets_cncf.txt" \
 --mmPSNPlotListData "refData/mmPSN_survivalCurvePlotList.rds"  --gep70File "exampleInputs/gep70scores_latest.csv"

# Patient 003
 Rscript spinReportFromRmd.R \
 --reportTemplateFile "refData/spinPatientReport.Rmd" \
 --sampleName "ISMMS03" --sexChrs "XX" --dateSampleCollected "01-01-2020" \
 --dateSampleSequenced "02-02-2020" --patientName="Jane Doe" --MRN '1234567' \
 --DOB "04-20-1969" --IgClass "IgG" --LightChainType "Kappa" --RevisedISS "VI" \
 --DiseaseStatus "Relapsed" --mainCloneResults "exampleInputs/ISMMS03/mainCloneResults.tsv" \
 --treeFile "exampleInputs/ISMMS03/sample.summ.json.gz" --zScoreTable "exampleInputs/zscores_latest.csv" \
 --sampleNameList "exampleInputs/sampleNameList.csv" --msiFileList "exampleInputs/msiFileList.csv" \
 --tmbFileList "exampleInputs/tmbFileList.csv" --mmPSNFile "exampleInputs/ISMMS03/Predicted_class.csv" \
 --cnaFile "exampleInputs/ISMMS03/cnaSubCloneDetailsTable.tsv" \
 --facetsCNCFResultsFile "exampleInputs/ISMMS03/tumor.facets_cncf.txt" \
 --mmPSNPlotListData "refData/mmPSN_survivalCurvePlotList.rds"  --gep70File "exampleInputs/gep70scores_latest.csv" \
 --tumorPurityPloidy "exampleInputs/ISMMS03/tumor.facets_output.txt" --subCloneResults "exampleInputs/ISMMS03/subCloneResults.tsv" \
 --snvSubCloneFile "exampleInputs/ISMMS03/snvSubCloneDetailsTable.tsv" --cnaSubCloneFile "exampleInputs/ISMMS03/cnaSubCloneDetailsTable.tsv" 
```