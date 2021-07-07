#!/bin/bash

docker build -t spinreport .
docker run -v "$PWD:/data" -v "/data1/users/restrp01/example_data:/exampleInputs" -v "/data2/mm/data/MSSM/sema4_clinical/21064632/processed:/processed" spinreport \
Rscript spin_report.r \
 --outdir '/data' \
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
 --exprFile "/processed/expression.prediction_engine.results.tsv" \
 --geneCoordinates "refData/hg38_geneCoordinates_TxDB.tsv" \
 --psnReferenceTable "refData/mmrf_mmpsn_subgroup_survival.tsv"  \
 --cytobandCoordinates 'refData/hg38_cytoband_coordinates.tsv'