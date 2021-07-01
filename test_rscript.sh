#!/bin/bash

docker build -t test_report .
docker run -v "$PWD:/data" test_report \
Rscript spinReportFromRmd.R \
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
 --drugDetailsRanked "/data2/mm/data/MSSM/sema4_clinical/21064632/processed/prediction_engine/all.drug_details.ranked.prediction_engine.results.tsv" \
 --rxBucketSummary "/data2/mm/data/MSSM/sema4_clinical/21064632/processed/prediction_engine/tier1.rxbucket_summary.ranked.prediction_engine.results.tsv" \
 --variantSupportSummary "/data2/mm/data/MSSM/sema4_clinical/21064632/processed/prediction_engine/variant_summary.ranked.prediction_engine.results.tsv" \
 --treeFile "/data2/mm/data/MSSM/sema4_clinical/21064632/processed/secondary/sample.summ.json.gz" \
 --zScoreTable "/data1/daphni/MSSM/generalRNA/zscores_latest.csv" \
 --cnaFile "/data2/mm/data/MSSM/sema4_clinical/21067200/processed/prediction_engine/cna.prediction_engine.results.tsv" \
 --exprFile "/data2/mm/data/MSSM/sema4_clinical/21064632/processed/prediction_engine/expression.prediction_engine.results.tsv" \
 --facetsCNCFResultsFile "/data2/mm/data/MSSM/sema4_clinical/21064632/processed/primary/tumor.facets_cncf.txt" \
 --mmPSNPlotListData "refData/mmPSN_survivalCurvePlotList.rds"  \
 --scarFile "exampleInputs/scarScore_latest.tsv" \
 --tmbFile "exampleInputs/tumorMutationBurden_latest.tsv" \
 --gep70File "exampleInputs/gep70scores_latest.tsv" \
 --tumorPurityPloidy "/data2/mm/data/MSSM/sema4_clinical/21064632/processed/primary/tumor.facets_output.txt" \
 --geneCoordinateBED "refData/hg38_geneCoordinates_TxDB.tsv"