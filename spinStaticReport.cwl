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

    geneCoordinates:
        type: File
        inputBinding:
            prefix: "--geneCoordinates"

    cytobandCoordinates:
        type: File
        inputBinding:
            prefix: "--cytobandCoordinates"
    
    psnReferenceTable:
        type: File?
        inputBinding:
            prefix: "--psnReferenceTable" 

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

    prognosticMarkerSummary:
        type: File
        inputBinding:
            prefix: "--prognosticMarkerSummary"    
    
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

    
    zScoreTable:
        type: File?
        inputBinding:
            prefix: "--zScoreTable"
    
    seliFile:
        type: File?
        inputBinding:
            prefix: "--seliFile"
            
    gep70File:
        type: File?
        inputBinding:
            prefix: "--gep70File" 

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

    facetsCNCFResultsFile:
        type: File?
        inputBinding:
            prefix: "--facetsCNCFResultsFile" 

    tumorPurityPloidy:
        type: File?
        inputBinding:
            prefix: "--tumorPurityPloidy"
            
    treeFile:
        type: File?
        inputBinding:
            prefix: "--treeFile"
    mmPSNFile:
        type: File?
        inputBinding:
            prefix: "--mmPSNFile"

baseCommand: [Rscript, spin_report.R]

outputs:
    report:
        type: File
        outputBinding:
            glob: "*.main.daphni_report.html"
    figures:
        type: File[]
        outputBinding:
            glob: "*.svg"
     kables:
        type: File[]
        outputBinding:
            glob: "*kable.daphni_report.html"
