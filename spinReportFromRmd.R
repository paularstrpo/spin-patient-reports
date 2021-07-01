options(stringsAsFactors=FALSE)
library(optparse)

# OPTION PARSER (INPUT FILE SPECS)
option_list <- list(
    make_option("--sampleName", type='character', help='[REQUIRED]', default=NA),
    make_option("--sexChrs", type='character', help='[REQUIRED]', default=NA),
    
    make_option("--drugDetailsRanked", type='character',help='[REQUIRED]', default=NA),
    make_option("--rxBucketSummary", type='character',help='[REQUIRED]', default=NA),
    make_option("--variantSupportSummary", type='character',help='[REQUIRED]', default=NA),
    
    make_option("--dateSampleCollected", type='character', default=NA),
    make_option("--dateSampleSequenced", type='character', default=NA),
    
    make_option("--patientName", type='character', default=NA),
    make_option("--MRN", type='character', default=NA),
    make_option("--DOB", type='character', default=NA),
    
    make_option("--IgClass", type='character', default=NA),
    make_option("--LightChainType", type='character', default=NA),
    make_option("--RevisedISS", type='character', default=NA),
    make_option("--DiseaseStatus", type='character', default=NA),
    
    make_option("--treeFile", type='character', default=NA),
    make_option("--zScoreTable", type='character', default=NA),
    
    make_option("--msiFile", type='character', default=NA),
    make_option("--scarFile", type='character', default=NA),
    make_option("--tmbFile", type='character', default=NA),
    make_option("--gep70File", type='character', default=NA),
    
    make_option("--somFile", type='character', default=NA),
    make_option("--snpFile", type='character', default=NA),
    make_option("--cnaFile", type='character', default=NA),
    make_option("--exprFile", type='character', default=NA),
    make_option("--pathwayFile", type='character', default=NA),
    make_option("--drugRepurposingFile", type='character', default=NA),
    
    make_option("--facetsCNCFResultsFile", type='character', default=NA),
    make_option("--tumorPurityPloidy", type='character', default=NA),
    
    make_option("--mmPSNFile", type='character', default=NA),
    make_option("--translocFile", type='character', default=NA),
    
    make_option("--mmPSNPlotListData", type='character', default=NA),
    make_option("--geneCoordinateBED", type='character', default=NA)
)

# Get command line options, if help option encountered print help and exit.
des <- "----------------------\n[DAPHNI v2.0] Precision Medicine Drug Report Generation\n----------------------"
opt <- parse_args(OptionParser(option_list=option_list, description = des))

# If required arguments are not present, exit with an error.
requiredArgs <- c('sampleName','sexChrs')
for (arg in requiredArgs)if(is.na(opt[arg])) stop(paste(arg, 'parameter must be provided. See script usage (--help)'))

reportPrefix <- paste0(opt$sampleName, '_report')
reportOutputFilename <- paste0(reportPrefix, ".html")

paramList <- list(
    sampleName = opt$sampleName,
    
    dateSampleCollected = opt$dateSampleCollected, 
    dateSampleSequenced = opt$dateSampleSequenced,
    dateReportGenerated = Sys.Date(),
    
    patientName = opt$patientName,
    sexChrs = opt$sexChrs,
    MRN = opt$MRN,
    DOB = opt$DOB,
    IgClass = opt$IgClass,
    LightChainType = opt$LightChainType,
    RevisedISS = opt$RevisedISS,
    DiseaseStatus = opt$DiseaseStatus,
    
    drugDetailsRanked = opt$drugDetailsRanked,
    rxBucketSummary = opt$rxBucketSummary,
    variantSupportSummary = opt$variantSupportSummary,
    
    
    treeFile = opt$treeFile,
    zScoreTable = opt$zScoreTable,
    
    msiFile = opt$msiFile,
    tmbFile = opt$tmbFile,
    scarFile = opt$scarFile,
    gep70File = opt$gep70File,
    
    mmPSNFile = opt$mmPSNFile, 
    translocFile = opt$translocFile,
    
    somFile = opt$somFile,
    snpFile = opt$snpFile,
    cnaFile = opt$cnaFile, 
    exprFile = opt$exprFile,
    pathwayFile = opt$pathwayFile, 
    drugRepurposingFile = opt$drugRepurposingFile, 
    
    tumorPurityPloidy = opt$tumorPurityPloidy, 
    facetsCNCFResultsFile = opt$facetsCNCFResultsFile, 
    
    mmPSNPlotListData = opt$mmPSNPlotListData,
    geneCoordinateBED = opt$geneCoordinateBED
    )

# grab the current working directory to force all the 
# working / temp directories in knitr to cooperate with docker
# and knit the darn document
wd = getwd()
rmarkdown::render('spinPatientReport.Rmd',
    params =  paramList,
    output_file = paste0(wd, '/', reportOutputFilename),
    knit_root_dir = wd, 
    intermediates_dir = wd)