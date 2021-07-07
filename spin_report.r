options(stringsAsFactors=FALSE)
library(optparse)

# OPTION PARSER (INPUT FILE SPECS)
option_list <- list(
    make_option("--patientName", help='patient name, in the format of first last, must be enclosed in quotes', type='character', default=NA),
    make_option("--sexChrs", help='Gender of patient, specified as either XX or XY [REQUIRED]', type='character', default=NA),
    make_option("--MRN", help = 'MRN of patient, must be enclosed in quotes', type='character', default=NA),
    make_option("--DOB", help = 'patient date of birth, in format of mm/dd/yyyy, must be enclosed in quotes', type='character', default=NA),
    make_option("--IgClass", type='character', default=NA),
    make_option("--LightChainType", type='character', default=NA),
    make_option("--RevisedISS", type='character', default=NA),
    make_option("--DiseaseStatus", type='character', default=NA),

    make_option("--sampleName", help='sample ID name for patient to run prediction engine on. [REQUIRED]', type='character', default=NA),
    make_option("--dateSampleCollected", help = 'date of sample collection, in format of mm/dd/yyyy, must be enclosed in quotes', type='character', default=NA),
    make_option("--dateSampleSequenced", help = 'date that sample sequencing was performed, in format of mm/dd/yyyy, must be enclosed in quotes', type='character', default=NA),

    make_option("--drugDetailsRanked", help = 'prediction engine full drug details table (tab-delimited; e.g. all.drug_details.ranked.prediction_engine.results.tsv)', type='character', default=NA),
    make_option("--rxBucketSummary",help = 'prediction engine tier 1 rx bucket summary table (tab-delimited; e.g. tier1.rxbucket_summary.ranked.prediction_engine.results.tsv)', type='character', default=NA),
    make_option("--variantSupportSummary", help = 'prediction engine variant support summary table (tab-delimited; e.g. variant_summary.ranked.prediction_engine.results.tsv)', type='character', default=NA),
    make_option('--prognosticMarkerSummary', help = 'prediction engine prognostic marker results table (tab-delimited; e.g. prognostic_markers.prediction_engine.results.tsv)', type='character', default=NA),
    
    make_option("--zScoreTable", help='Table with expression zscores for all patients, including the patient of interest (tab-delimited; e.g. zscores_latest.tsv)', type='character', default=NA),
    make_option("--gep70File", help='Table containing gep70 scores at cohort level, including the patient of interest (tab-delimited; e.g. gep70scores_latest.tsv)',type='character',  default=NA),
    make_option("--seliFile", help='Table containing selinexor signature data at cohort level, including the patient of interest (tab-delimited; e.g. selinescores_latest.tsv)', type='character', default=NA),

    make_option("--msiFile",  help='Table containing MSI scores at cohort level, including the patient of interest (tab-delimited; e.g. msiScore_latest.tsv)', type='character', default=NA),
    make_option("--scarFile", help='Table containing scar scores at cohort level, including the patient of interest (tab-delimited; e.g. scarScore_latest.tsv)', type='character', default=NA),
    make_option("--tmbFile", help='Table containing tumor mutation burdens at cohort level, including the patient of interest (tab-delimited; e.g. tumorMutationBurden_latest.tsv)', type='character', default=NA),

    make_option("--somFile", help = 'prediction engine results for somatic mutations (tab-delimited; e.g. somatic_mutation.prediction_engine.results.tsv)', type='character', default=NA),
    make_option("--cnaFile", help = 'prediction engine results for CNVs (tab-delimited; e.g. cna.prediction_engine.results.tsv)', type='character', default=NA),
    make_option("--exprFile",help = 'prediction engine results for expression (tab-delimited; e.g. expression.prediction_engine.results.tsv)', type='character', default=NA),

    make_option("--snpFile", type='character', default=NA),
    make_option("--pathwayFile", type='character', default=NA),
    make_option("--drugRepurposingFile", type='character', default=NA),

    make_option("--mmPSNFile", help='file containing predicted mm-psn subgroup information for the sample (comma-delimited; e.g. Predicted_class.csv' , type='character', default=NA),
    make_option("--treeFile", help='gzip-compressed json file with clonal tree structure (e.g. sample.sum.json.gz)', type='character', default=NA),
    make_option("--facetsCNCFResultsFile", help='facets cncf file (tab-delimited; e.g. .tumor.facets_cncf.txt)', type='character', default=NA),
    make_option("--tumorPurityPloidy", help = 'facets output file containing purity and ploidy (tab-delimited; e.g. .tumor.facets_output.txt)', type='character', default=NA),
    
    make_option('--cytobandCoordinates', help='Table with full cytoband coordinates (tab-delimited; e.g. hg38_cytoband_coordinates.tsv) [REQUIRED]', type='character', default=NA),
    make_option("--psnReferenceTable", help='Table with MMRF PSN survival per subgroup information (tab-delimited; e.g. mmrf_mmpsn_subgroup_survival.tsv) [REQUIRED]', type='character', default=NA),
    make_option("--geneCoordinates", help='Table with genomic start and end coordinates for each gene (tab-delimited; e.g hg38_geneCoordinates_TxDB.tsv) [REQUIRED]', type='character', default=NA)
)

# Get command line options, if help option encountered print help and exit.
des <- "----------------------\n[DAPHNI v2.0] Precision Medicine Report Generation\n----------------------"
opt <- parse_args(OptionParser(option_list=option_list, description = des))

# If required arguments are not present, exit with an error.
requiredArgs <- c('sampleName','sexChrs', 'geneCoordinates', 'cytobandCoordinates', 'psnReferenceTable')
for (arg in requiredArgs)if(is.na(opt[arg])) stop(paste(arg, 'parameter must be provided. See script usage (--help)'))

# libraries
libs <- c(
    'tidyverse', 'RColorBrewer', 
    'rmarkdown', 'knitr',
    'RJSONIO', 'jsonlite', 'svglite',  
    'igraph', 'ggraph', 'tidygraph',
    'survival', 'survminer',
    'ComplexHeatmap', 'circlize', 'grid',
    'formattable', 'kableExtra'
)

lapply(libs, library, character.only=TRUE)
rm(libs)
# ------------------- #

# load in functions
source('functions.r')

# ---- LOAD DATA ---- #
# Reference data
chromosomes <- paste0('chr', c(1:22, uniqchars(toupper(opt$sexChrs))))
cyto <- read_tsv(opt$cytobandCoordinates)
cyto <- cyto[cyto$chrom %in% chromosomes, ]

geneCoords <- read_tsv(opt$geneCoordinates)
surv_psn <- read.delim(opt$psnReferenceTable, sep='\t')

# Patient Clincal Info
patientName <- opt$patientName
sexChrs <- opt$sexChrs
MRN <- opt$MRN
DOB <- opt$DOB
IgClass <- opt$IgClass
LightChainType <- opt$LightChainType
RevisedISS <- opt$RevisedISS
DiseaseStatus <- opt$DiseaseStatus

# Sample Name & Info
sampleName <- opt$sampleName
rnaSampleID <- gsub('-', '_', gsub('DNA', 'RNA', sampleName), fixed=TRUE)
dnaSampleID <- gsub('RNA', 'DNA', sampleName)

dateSampleCollected = opt$dateSampleCollected 
dateSampleSequenced = opt$dateSampleSequenced
dateReportGenerated = as.character(Sys.Date())

# per-alteration results (from prediction engine)
exprRes <-loadIfExists(opt$exprFile)
cnaRes <- loadIfExists(opt$cnaFile)
somRes <- loadIfExists(opt$somFile)

snpRes <- loadIfExists(opt$snpFile)
pathwayRes <- loadIfExists(opt$pathwayFile)
drugRepurposingRes <- loadIfExists(opt$drugRepurposingFile)

# mm-psn
mmPSNTable <- loadIfExists(opt$mmPSNFile, sep=',')

# Rx Tables
drugDetailsRankedTable <- loadIfExists(opt$drugDetailsRanked)
rxBucketSummaryTable <- loadIfExists(opt$rxBucketSummary)
variantSupportSummaryTable <- loadIfExists(opt$variantSupportSummary)

# cohort-level files:
# - RNA
zscoredf <- loadIfExists(opt$zScoreTable)
gep70df <- loadIfExists(opt$gep70File)
seliMat <- loadIfExists(opt$seliFile)
# - DNA
tmbMat <- loadIfExists(opt$tmbFile)
msiMat <- loadIfExists(opt$msiFile)
scarMat <- loadIfExists(opt$scarFile)


# tumor purity and ploidy
tumorPurityPloidy <- loadIfExists(opt$tumorPurityPloidy)

if(!is.na(tumorPurityPloidy)){
  purity <- tumorPurityPloidy$purity[1]
  ploidy <- tumorPurityPloidy$ploidy[1]
  
}else{
  purity <- NA
  ploidy <- NA
}

# facets CNV
facetsRes <- loadIfExists(opt$facetsCNCFResultsFile)
# ------------------- #

# -- MISC OPTIONS / SETUP -- #

# report output options
outputPrefix <- paste(sampleName, dateReportGenerated, sep='.')

# set up color code for plots etc.
colors <- c(`Somatic Mutation`="black",
            `Germline/Non-Somatic SNP`='gray70',
            `CNA Deletion`="#de64b1",
            `CNA Amplification`="#d808BC",
            `Gene under-expression`="#15D6D6", 
            `Gene over-expression`="#00AEEF", 
            `Downregulated Immunotherapy Target`="#b2e80e",
            `Upregulated Immunotherapy Target`="#07ad60"
            )

colordf <- data.frame(color=colors, type=names(colors))

# which gene expression markers to look at?
expMarkerGenes <- c('BCL2', 'TNFRSF17', 'GPRC5D', 'MCL1')

# ggplot theme
theme <- theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
               axis.line = element_blank(), text = element_text(size=9, color='black'),
               plot.title = element_text(hjust = 0.5, face='bold'))

actBED <- list() # list object for circos plot
# ------------------- #

# -- PROCESS  DATA -- #

# PREPARE RX RECC TABLES & VARIANT STATEMENTS
bucketSummaryTable <- drugDetailsRankedTable %>% 
  left_join(rxBucketSummaryTable) %>%
  left_join(variantSupportSummaryTable) %>%
  filter(bucket.summary_score > 0 & rx.summary_score > 0) %>%
  arrange(-bucket.summary_score) %>% 
  dplyr::select(Rx.Bucket, Rx.Reccomendation, variant_name, variant_statement, Tier) %>%
  unique() %>% 
  mutate(variant_statement = paste0('**', variant_name,':** ', variant_statement ),
  drug_reccomendations = paste0('&#9733; ', Rx.Reccomendation)
  ) %>% group_by(Rx.Bucket, Tier) %>%
  summarise(variant_statement = paste(unique(variant_statement), collapse='<br> '),
            drug_reccomendations = paste(unique(drug_reccomendations), collapse='<br> ')) %>%
  ungroup() %>% unique()


# somatic mutations
if(!is.na(somRes)){
  somDetailTable <- somRes %>%
    dplyr::select(Hugo_Symbol, Chromosome, Start_Position, End_Position, protein_change, cDNA_Change, tumor_f, t_alt_count, Variant_Classification) %>%
    mutate(ID=paste(Hugo_Symbol, protein_change)) %>%
    unique()
  
  somBED <- data.frame(chr=somDetailTable$Chromosome,
                      start=somDetailTable$Start_Position,
                      end=somDetailTable$Start_Position,
                      gene=somDetailTable$ID,
                      type='Somatic Mutation')
  
  actBED[['SOM']] <- unique(somBED)
  
  # Also: make the pretty somatic mutation chart
  if(nrow(somDetailTable) > 0) {
    donutPlotList <- sapply(somDetailTable$tumor_f, vaf.donut)
    vaftxt <-  paste0(round(somDetailTable$tumor_f*100, 2), "%")
    somDetailTable$VAF.format <- vaftxt
    somDetailTable$donutPlot <- donutPlotList
    
    somResKable <- tibble(`Gene` = somDetailTable$Hugo_Symbol, 
          `Variant Classification` = somDetailTable$Variant_Classification,
          `cDNA Change` = somDetailTable$cDNA_Change,
          `Protein Change` =  somDetailTable$protein_change, 
          `VAF` = donutPlotList,
          `Tumor Read Depth` = formattable::color_bar(color="#00AEEF")(somDetailTable$t_alt_count)) %>% 
      kableExtra::kable(escape=FALSE) %>%
      kableExtra::kable_styling(full_width = TRUE, html_font = "Arial", font_size=14,
                                bootstrap_options = c("hover", "condensed", "responsive")) %>%
      kableExtra::row_spec(1:nrow(somDetailTable), extra_css = "border-top: 1px solid lightgray;")
    
    # save kable
    save_kable(somResKable, paste0(outputPrefix, '.somatic_mutations.kable.spinreport.html'))
  }
}

# non-somatic snvs
if(!is.na(snpRes)){
  snpDetailTable <- snpRes %>%
    dplyr::select(Hugo_Symbol, Chromosome, Start_Position, End_Position, protein_change, cDNA_Change, tumor_f, t_alt_count, Variant_Classification, ID) %>%
    mutate(ID=paste(Hugo_Symbol, protein_change, ifelse(startsWith(ID, 'T'), '(tumor)', '(normal)'))) %>%
    unique()
  
  snpBED <- data.frame(chr=snpDetailTable$Chromosome,
                       start=snpDetailTable$Start_Position,
                       end=snpDetailTable$Start_Position,
                       gene=snpDetailTable$ID,
                       type='Germline/Non-Somatic SNP')
  
  actBED[['SNP']] <- unique(snpBED)
  
  # Also: make the pretty somatic mutation chart
  if(nrow(snpDetailTable) > 0) {
    donutPlotList <- sapply(snpDetailTable$tumor_f, vaf.donut)
    vaftxt <-  paste0(round(snpDetailTable$tumor_f*100, 2), "%")
    snpDetailTable$VAF.format <- vaftxt
    snpDetailTable$donutPlot <- donutPlotList
    
    snpResKable <- tibble(`Gene` = snpDetailTable$Hugo_Symbol,
                          `Source Tissue` = ifelse(grepl(pattern='tumor', x=snpDetailTable$ID), 'Tumor', 'Normal'),
                          `Variant Classification` = snpDetailTable$Variant_Classification,
                          `cDNA Change` = snpDetailTable$cDNA_Change,
                          `Protein Change` =  snpDetailTable$protein_change, 
                          `VAF` = donutPlotList,
                          `Read Depth` = formattable::color_bar(color="#00AEEF")(snpDetailTable$t_alt_count)) %>% 
      kableExtra::kable(escape=FALSE) %>%
      kableExtra::kable_styling(full_width = TRUE, html_font = "Arial", font_size=14,
                                bootstrap_options = c("hover", "condensed", "responsive")) %>%
      kableExtra::row_spec(1:nrow(snpDetailTable), extra_css = "border-top: 1px solid lightgray;")
    
    # save kable
    save_kable(snpResKable,  paste0(outputPrefix, '.germline_snp.kable.spinreport.html'))
  }
}

# copy number variants
if(!is.na(cnaRes)){
  cnaTable <- cnaRes %>%
    dplyr::select(Hugo_Symbol, Chromosome, Start_Position, End_Position, band, Variant_Classification) %>%
    mutate(type=ifelse(grepl('amp', Variant_Classification, ignore.case=TRUE),  "Amplification", "Deletion"),
           ID=paste(Hugo_Symbol, type),
           typeFormatted=sapply(Variant_Classification, function(x){paste(ifelse(grepl('amp', x, ignore.case=TRUE),  "&#8593;", "&#8595;"), x)})) %>%
    unique()
  
  cnvBED <- data.frame(chr=cnaTable$Chromosome,
                       start=cnaTable$Start_Position,
                       end=cnaTable$Start_Position,
                       gene=cnaTable$ID,
                       type=paste('CNA', cnaTable$type))

  actBED[['CNA']] <- unique(cnvBED)

  cnvDetailTable <- tibble(Locus = paste(cnaTable$Chromosome, as.character(cnaTable$band)), Gene = cnaTable$Hugo_Symbol, Type = cnaTable$typeFormatted) %>%
    unique() %>% group_by(Locus, Type) %>% summarise(Genes = paste(paste('&#9733;', unique(na.omit(Gene))), collapse="<br> ")) %>%
    ungroup() %>% unique()

  cnvResKable <-  cnvDetailTable %>% 
    kableExtra::kable(escape = FALSE) %>%
    kableExtra::kable_styling(full_width = TRUE, html_font = "Arial", font_size=14,
                            bootstrap_options = c("hover", "condensed", "responsive")) %>%
    kableExtra::row_spec(1:nrow(cnvDetailTable), extra_css = "border-top: 1px solid lightgray;")
  
  # save kable
  save_kable(cnvResKable, paste0(outputPrefix, '.cnv.kable.spinreport.html'))
}


# expression
if(!is.na(exprRes) && length(exprRes) > 0){
  exprsTable <- data.frame(unique(exprRes[, c('Hugo_Symbol', 'variant_name', 'Zscore', 'Variant_Classification', 'variant_statement')]))
  exprsTable <- exprsTable %>% left_join(geneCoords, by=c('Hugo_Symbol'='Gene'))
  
  exprBED <- data.frame(chr=exprsTable$seqnames,
                        start=exprsTable$start,
                        end=exprsTable$start,
                        gene=exprsTable$variant_name,
                        type=paste('Gene', gsub('ed', 'ion', exprsTable$Variant_Classification))
                        )
  
  actBED[['Expression']] <- exprBED
  
  exprPlotList <- sapply(exprsTable$Hugo_Symbol, expr.dist)
  exprResKable <- tibble(`Gene` = exprsTable$Hugo_Symbol,
                         `Variant` = exprsTable$Variant_Classification,
                         `Expression in MSSM Cohort` = exprPlotList) %>%
    kableExtra::kable(escape = FALSE) %>% # note to self: use column_spec to limit width of plots explicitly
    kableExtra::kable_styling(full_width = FALSE, html_font = "Arial", font_size=14,
                                bootstrap_options = c("hover", "condensed", "responsive")) %>%
    kableExtra::row_spec(1:nrow(exprsTable), extra_css = "border-top: 1px solid lightgray;")
  
  # save kable
  save_kable(exprResKable, paste0(outputPrefix, '.expression.kable.spinreport.html'))
}

# ------------------- #

# MAKE  COHORT  PLOTS #

# tumor mutation burden
if(!is.na(tmbMat)){
  tmbMat$isPt <- ifelse(grepl(dnaSampleID, tmbMat$sampleID), 'yes', 'no')
  mutBurden <- round(tmbMat[tmbMat$isPt=='yes','nonsilentMutationalBurden'], 2)
  
  tmbPlot <- ggplot(tmbMat) + 
    aes(x=reorder(sampleID, nonsilentMutationalBurden), y=nonsilentMutationalBurden, fill=isPt, alpha=isPt, width=ifelse(isPt=='yes', 5, 1)) + 
    geom_bar(stat='identity',position='identity', color='white', size=0, show.legend=FALSE) + 
    scale_fill_manual(values=c(yes='#00AEFF', no='gray80')) +
    scale_alpha_manual(values=c(yes=1, no=0.5)) + 
    theme_classic() +  theme +
    labs(x='Patient Samples', y = 'Tumor Mutation Burden\n(Nonsilent Muts/MB)')

  svglite(paste0(outputPrefix, '.tmb.plot.daphni_report.svg'), height=2.5, width=3)
  tmbPlot
  dev.off()
}

# scar score
if(!is.na(scarMat)){
  scarMat$isPt <-  ifelse(grepl(dnaSampleID, scarMat$sampleID), 'yes', 'no')
  scarScore <- round(scarMat[scarMat$isPt=='yes','HRD.sum'], 2)
  
  scarPlot <- ggplot(scarMat) + aes(x=reorder(sampleID, HRD.sum), y = HRD.sum, fill=isPt, alpha=isPt, width=ifelse(isPt=='yes', 5, 1)) + 
    geom_bar(stat='identity',position='identity', color='white', size=0, show.legend=FALSE) + 
    scale_fill_manual(values=c(yes='#00AEFF', no='gray80')) + 
    scale_alpha_manual(values=c(yes=1, no=0.5)) + 
    theme_classic() + theme +
    labs(x='Patient Samples', y = 'scarHRD Score')

  svglite(paste0(outputPrefix, '.scar.plot.daphni_report.svg'), height=2.5, width=3)
  scarPlot
  dev.off()
}

# msi
if(!is.na(msiMat)){
  msiMat$isPt <- ifelse(grepl(dnaSampleID, msiMat$sampleID), 'yes', 'no')
  msiMat$msi_category <- factor(ifelse(msiMat$msi < 3, 'MSS', ifelse(msiMat$msi < 10, 'MSI-Low', 'MSI-High')),
                                levels=c('MSS', 'MSI-Low', 'MSI-High'), ordered=TRUE)
  
  msiVal <- round(as.numeric(msiMat[msiMat$isPt=='yes','msi']), 2)
  msiCat <- as.character(msiMat[msiMat$isPt=='yes','msi_category'])
  
  msiPlot <- ggplot(msiMat) + aes(x=reorder(sampleName, msi), y=msi, fill=isPt, alpha=isPt, width=ifelse(isPt=='yes', 5, 1)) +
    geom_hline(yintercept = 10, linetype='dashed', color='black', size=0.25) +
    geom_hline(yintercept = 3, linetype='dashed', color='black', size=0.25) +
    geom_bar(stat='identity',position='identity', color='white', size=0, show.legend=FALSE) +
    geom_rug(aes(color=msi_category, width=NULL, alpha=NULL), sides='b', show.legend=FALSE) + 
    scale_fill_manual(values=c(yes='#00AEFF', no='gray80')) + 
    scale_alpha_manual(values=c(yes=1, no=0.5)) + 
    theme_classic() + theme +
    labs(x='Patient Samples', y = 'MSI')

  svg(paste0(outputPrefix, '.msi.plot.daphni_report.svg'), height=2.5, width=3)
  msiPlot
  dev.off()
}

# gep 70
if(!is.na(gep70df)){
  gep70df$sampleName <- rownames(gep70df)
  gep70df$score <- as.numeric(gep70df$x)
  gep70df$isPt <- ifelse(grepl(rnaSampleID, gep70df$sampleName), 'yes', 'no')
  gep70df$risk_category <- factor(ifelse(gep70df$score < 40, 'low', 
                                    ifelse(gep70df$score < 45, 'intermediate low', 
                                           ifelse(gep70df$score < 50, 'intermediate high', 'high'))),
                    levels=c('low', 'intermediate low', 'intermediate high', 'high'), ordered=TRUE)
   
  riskScore <- round(gep70df$score[gep70df$isPt == 'yes'], 2)
  riskCat <- as.character(gep70df$risk_category[gep70df$isPt == 'yes'])
  
  gepPlot <- ggplot(gep70df) + 
    aes(x=reorder(sampleName, score), y=score, fill=isPt, alpha=isPt, width=ifelse(isPt=='yes', 5, 1)) +
    geom_bar(stat='identity',position='identity', color='white', size=0, show.legend=FALSE) + 
    geom_rug(aes(color=risk_category, width=NULL, alpha=NULL), sides='b', show.legend=FALSE) + 
    scale_fill_manual(values=c(yes='#00AEFF', no='gray80')) + 
    scale_alpha_manual(values=c(yes=1, no=0.5)) + 
    theme_classic() + theme +
    labs(x='Patient Samples', y = 'GEP70 Score')

  svglite(paste0(outputPrefix, '.gep70.plot.daphni_report.svg'), height=2.5, width=3)
  gepPlot
  dev.off()

}

# selinexor signature
if(!is.na(seliMat)){
  seliMat$sampleID <- rownames(seliMat)
  seliMat$isPt <-  ifelse(grepl(rnaSampleID, seliMat$sampleID), 'yes', 'no')
  seliScore <- round(seliMat[seliMat$isPt=='yes','score'], 2)
  
  seliPlot <- ggplot(seliMat) + aes(x=reorder(sampleID, score), y = score, fill=isPt, alpha=isPt, width=ifelse(isPt=='yes', 5, 1)) + 
    geom_bar(stat='identity',position='identity', color='white', size=0, show.legend=FALSE) + 
    scale_fill_manual(values=c(yes='#00AEFF', no='gray80')) + 
    scale_alpha_manual(values=c(yes=1, no=0.5)) + 
    theme_classic() + theme +
    labs(x='Patient Samples', y = 'SelineScore Gene Signature Expression')


  svglite(paste0(outputPrefix, '.selinescore.plot.daphni_report.svg'), height=2.5, width=3)
  seliPlot
  dev.off()
}

# expression plots for genes of interest
expressionPlotList <- list()
if(!is.na(zscoredf)){
  zscoreMelt <- zscoredf %>%
    dplyr::select(-X) %>%
    filter(Gene %in% expMarkerGenes) %>%
    pivot_longer(cols = colnames(zscoredf)[3:ncol(zscoredf)], 
                 names_to='sampleName', values_to='zscore') %>%
    pivot_wider(names_from='Gene', values_from='zscore') %>%
    data.frame()

  for(gene in expMarkerGenes){
    df <- data.frame(sampleName=zscoreMelt$sampleName,
                     zscore=zscoreMelt[, gene],
                     geneName=gene,
                     isPt = ifelse(grepl(rnaSampleID, zscoreMelt$sampleName), 'yes', 'no'))
    
    
    plot <- ggplot(df) +  aes(x=reorder(sampleName, zscore), y=zscore, fill=isPt, alpha=isPt, width=ifelse(isPt=='yes', 5, 1)) +
      geom_bar(stat='identity', position='identity', color='white', size=0, show.legend=FALSE) + 
      scale_fill_manual(values=c(yes='#00AEFF', no='gray80')) + theme_classic() + 
      scale_alpha_manual(values=c(yes=1, no=0.5)) + theme +
      labs(x='Patient Samples', y = paste(gene, 'Expression'), 
           title=paste0(gene, ': z-score = ', round(df$zscore[df$isPt=='yes'], 2)))
    

    svglite(paste0(outputPrefix,'.', gene, '_expression.plot.daphni_report.svg'), height=2.5, width=3)
    plot
    dev.off()
    
    expressionPlotList[[gene]] <- plot
    rm(df, plot)
  }
}

# ------------------- #


# -- SAMPLE  PLOTS -- #

# mm-psn
if(!is.na(mmPSNTable)){
    group <- mmPSNTable$Subgroup[1]
    mmPSNPredictedClass <- group
    subgroups <- unique(surv_psn$SubGroup)

    # move factor level for group of interest to very end so it's on top
    surv_psn$SubGroup <- factor(surv_psn$SubGroup, c(subgroups[subgroups != group], group))
    fit <- survfit(Surv(ttcos, censos) ~ SubGroup, data = surv_psn)

    psn_plt <- ggsurvplot(fit, data = surv_psn, risk.table = FALSE, pval = FALSE, censor=FALSE,
                    palette = c(rep("gray80",length(subgroups)-1), "#00AEEF"),
                    ylab="Overall Survival", legend.title="SubGroup", pval.size = 18, legend="none",
                    ggtheme = theme_classic())
    # divide by average number of days in a month (approx 365/12 =~30.417) and round to two decimals for presentation purposes
    psnSubgroupMedianOS <- round(median(as.numeric(surv_psn[surv_psn$SubGroup == group,'ttcos'])) / (365/12), 2)
}

# clonal tree
if(!is.na(opt$treeFile)) {
  # get top clonal tree from structure file
  # based on the llh value closest to zero
  structFile <- RJSONIO::fromJSON(opt$treeFile)
  llh <- sapply(structFile$trees, function(x){abs(x$llh)})
  cloneTree <- structFile$trees[[which(llh == min(llh))]]
  treeMatrix <- as.matrix(tree_mat(cloneTree$structure))
  
  rm(structFile, llh) # clean-up
  
  # From the tree structure matrix, annotate a clone table
  # with tree information, and clone specific ccf, cnv burden, snv burden,
  # and each clone's downstream descendants
  cloneData <- data.frame(
    from=as.character(treeMatrix[,1]),
    to=as.character(treeMatrix[,2])
    ) %>%
    mutate(clone=from) %>% 
    group_by(clone) %>%
    summarise(descendant_clones = paste(unique(to), collapse='|')) %>%
    right_join(get_clones(cloneTree$populations)) %>%
    mutate(
      tree.linearity_index = cloneTree$linearity_index,
      tree.llh = cloneTree$llh,
      tree.clustering_index = cloneTree$clustering_index,
      tree.branching_index = cloneTree$branching_index
    )
  
  clones <- cloneData$clone
  
  tree <- igraph::graph_from_edgelist(treeMatrix+1, directed = FALSE) %>% 
    as_tbl_graph() %>% create_layout(layout = 'tree', root=1, flip.y=FALSE) %>% 
    mutate(cloneName=.ggraph.index-1)
  
  clonalTreePlot <- ggraph(tree) +
    geom_edge_link(size=8) + geom_node_point(size = 12) +
    geom_node_text(aes(label=cloneName), color='white', size=9) +
    coord_flip() + theme_void() + theme(legend.position = 'none') +
    guides(color = FALSE, size = FALSE)

  svglite(paste0(outputPrefix, '.clonal_tree.plot.daphni_report.svg'), height=6, width=7)
  clonalTreePlot
  dev.off()
}

# global copy number profile
if(!is.na(facetsRes)){
  # Obtain cytogenetic information for the genome of interest
  facetsResLong <- facetsRes %>% pivot_longer(cols=c("tcn.em", "lcn.em"),names_to='allele', values_to='copynumber') %>%
    mutate(chrom=factor(paste0('chr', ifelse(chrom==23, 'X', chrom)), levels=paste0('chr', c(1:22, 'X', 'Y', 'M')))) %>%
    filter(chrom %in% chromosomes) %>%
    filter(!is.na(copynumber))
  
  cnProfilePlot <- ggplot(data = facetsResLong) + 
    geom_hline(yintercept=2, linetype='dashed', color='gray80', size=0.5) +
    geom_segment(size=1.5, aes(x = as.numeric(start), xend = as.numeric(end), y=as.numeric(copynumber), yend=as.numeric(copynumber), color=allele, fill=NULL) ) +
    geom_rect(data=cyto, color='black',size=0.075, show.legend=FALSE,  aes(xmin=as.numeric(start), xmax=as.numeric(end), 
                  ymin=-0.25, ymax = -0.5, color=NULL, fill=gieStain)) +
    scale_color_manual(values=c(tcn.em='#00AEFF', lcn.em='gray40'), labels=c(tcn.em='Total Copy Number', lcn.em='Minor Allele Copy Number'), name='') +
    scale_fill_manual(values=c(acen='firebrick', gneg='white', gpos='black', 
                               gpos25='gray70', gpos50='gray50', gpos75='gray30', gvar='gray20', stalk='tomato')) +
    facet_wrap(~chrom, scales='free_x', nrow=1, strip.position = 'bottom') + labs(x='Genomic Coordinate', y='Copy Number') + theme_minimal() + 
    theme(title=element_text(color='black', size=22), text=element_text(color='black', size=22),
          axis.text.x=element_blank(), panel.grid=element_blank(),
          strip.background=element_rect(fill='white', color='white'),  strip.placement='outside',
          legend.position='top', legend.justification='center', 
          axis.line.y.left=element_line(color='black'), 
          axis.line.x.bottom=element_line(color='black'))
  
  svglite(paste0(outputPrefix, '.global_cna_profile.plot.daphni_report.svg'), height=5, width=22)
  cnProfilePlot
  dev.off()
}

if(length(actBED) > 0){
  actBED <- do.call(rbind, actBED) %>% inner_join(colordf)
  actBED <- actBED[!is.na(actBED$start)&!is.na(actBED$end)&!is.na(actBED$chr),]
  
  print_circos <- function(){

    par(mar=c(0,0,0,0), oma=c(0,0,0,0))
    circos.par(track.height = 0.05, track.margin=c(0,0), cell.padding=c(0,0,0,0), canvas.xlim = c(-1.25, 1.25), canvas.ylim = c(-1.25, 1.25))
    circos.initializeWithIdeogram(species="hg38", chromosome.index=chromosomes, plotType = NULL)
    circos.genomicLabels(actBED, labels.column = 4, side = "outside", cex=1, niceFacing=TRUE, col = actBED[[6]], line_col = actBED[[6]])

    circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
        chr = gsub('chr', '', CELL_META$sector.index)
        xlim = CELL_META$xlim
        ylim = CELL_META$ylim
        circos.rect(xlim[1], 0, xlim[2], 1, col = 'black')
        circos.text(mean(xlim), mean(ylim), chr, cex = 0.9, col = "white", facing = "inside", font=2, niceFacing = TRUE)
    }, track.height = 0.12, bg.border = NA)
    circos.clear()

  }

  svglite(paste0(outputPrefix, '.circos_landscape.plot.daphni_report.svg'), height=8, width=8)
  print_circos()
  legend("bottomleft", names(colors), pch=19, col=colors, cex=0.9, bty = "n")
  dev.off()
}
# -------------------- #

# - RENDER  MARKDOWN - #

# grab the current working directory to force all the 
# working / temp directories in knitr to cooperate with docker
# and knit the darn document
wd = getwd()

rmarkdown::render('report_template.rmd',
    output_file = paste0(outputPrefix, '.main.daphni_report.html'),
    knit_root_dir = wd, 
    intermediates_dir = wd)