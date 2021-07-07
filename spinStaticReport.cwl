#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

requirements:
    DockerRequirement:
        dockerPull: "sinaiiidgst/spinreport:latest"
    InlineJavascriptRequirement: {}
  

inputs:

    # required arguments
    sampleName:
        type: string
        inputBinding:
            prefix: "--sampleName" 

    geneCoordinateBED:
        type: File
        inputBinding:
            prefix: "--geneCoordinateBED" 

    sexChrs:
        type: string
        inputBinding:
            prefix: "--sexChrs"

    drugDetailsRanked:
        type: File
        inputBinding:
            prefix: "--drugDetailsRanked"
            
    rxBucketSummary:
        type: File
        inputBinding:
            prefix: "--rxBucketSummary"
            
    variantSupportSummary:
        type: File
        inputBinding:
            prefix: "--variantSupportSummary"    
    
    # These parameters should be included if data is available
    dateSampleCollected:
        type: string?
        inputBinding:
            prefix: "--dateSampleCollected"
    
    dateSampleSequenced:
        type: string?
        inputBinding:
            prefix: "--dateSampleSequenced" 

    patientName:
        type: string?
        inputBinding:
            prefix: "--patientName" 

    DOB:
        type: string?
        inputBinding:
            prefix: "--DOB" 

    MRN:
        type: string?
        inputBinding:
            prefix: "--MRN" 

    IgClass:
        type: string?
        inputBinding:
            prefix: "--IgClass" 

    LightChainType:
        type: string?
        inputBinding:
            prefix: "--LightChainType" 

    RevisedISS:
        type: string?
        inputBinding:
            prefix: "--RevisedISS"

    DiseaseStatus:
        type: string?
        inputBinding:
            prefix: "--DiseaseStatus" 
    
    translocFile:
        type: File?
        inputBinding:
            prefix: "--translocFile"
            
    treeFile:
        type: File?
        inputBinding:
            prefix: "--treeFile"
    
    zScoreTable:
        type: File?
        inputBinding:
            prefix: "--zScoreTable" 

    mmPSNFile:
        type: File?
        inputBinding:
            prefix: "--mmPSNFile"

    somFile:
        type: File?
        inputBinding:
            prefix: "--somFile"

    snvFile:
        type: File?
        inputBinding:
            prefix: "--snvFile"

    cnaFile:
        type: File?
        inputBinding:
            prefix: "--cnaFile" 

    exprFile:
        type: File?
        inputBinding:
            prefix: "--exprFile"
    
    pathwayFile:
        type: File?
        inputBinding:
            prefix: "--pathwayFile" 

    drugRepuposingFile:
        type: File?
        inputBinding:
            prefix: "--drugRepurposingFile" 

    msiFile:
        type: File?
        inputBinding:
            prefix: "--msiFile"

    tmbFile:
        type: File?
        inputBinding:
            prefix: "--tmbFile"
    
    scarFile:
        type: File?
        inputBinding:
            prefix: "--scarFile"
            
    gep70File:
        type: File?
        inputBinding:
            prefix: "--gep70File" 


    facetsCNCFResultsFile:
        type: File?
        inputBinding:
            prefix: "--facetsCNCFResultsFile" 

    mmPSNPlotListData:
        type: File?
        inputBinding:
            prefix: "--mmPSNPlotListData" 

    tumorPurityPloidy:
        type: File?
        inputBinding:
            prefix: "--tumorPurityPloidy" 

baseCommand: [Rscript, /bin/spinReportFromRmd.R]

outputs:
    report:
        type: File
        outputBinding:
            glob: $(inputs.sampleID + "_report.html")
    figures:
        type: File[]
        outputBinding:
            glob: $(inputs.sampleID + "_report_files/figure-html/" + "*.svg")
     kables:
        type: File[]
        outputBinding:
            glob: "*formatted_kable.spinreport.html"
