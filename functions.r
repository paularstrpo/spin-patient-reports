# misc. functions
df.empty <- function(df){is.data.frame(df) && nrow(df)==0}
unlist.cell <- function(x){unique(unlist(strsplit(x, '|', fixed=TRUE)))}
uniqchars <- function(x) unique(strsplit(x, "")[[1]]) 

# colorize text in-line
colorize <- function(x, color) {
  if(knitr::is_latex_output()) {
    sprintf("\\textcolor{%s}{%s}", color, x)
  } else if(knitr::is_html_output()) {
    sprintf("<span style='color: %s;'>%s</span>", color, 
      x)
  } else x
}

# check if a file is empty
file.empty <- function(file){
  con <- file(file, open="r")
  res <- ! grepl('[^[:space:]]', readLines(con, n=1))
  close.connection(con)
  return(res)
}

# Read table in only if the file exists and is not empty (otherwise return NA)
loadIfExists <- function(inputFile, sep='\t', ...){
  if(!is.na(inputFile)){
    if(file.exists(inputFile) && !file.empty(inputFile)){
      return(read.delim(inputFile, sep = sep, header = TRUE, ...))
    }else{
      return(NA)
    }
  } else{
    return(NA)
  }
}

# make vaf donut as svg text
vaf.donut <- function(x){ # x=vaf
  s <- svglite::svgstring(); s()
  dat <- data.frame(fraction=c(x, (1-x)), category=c('ALT', 'REF'))
  dat$ymax = cumsum(dat$fraction)
  dat$ymin = c(0, head(dat$ymax, n=-1))
  
  pie <- ggplot(dat) +
    aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category) +
    geom_rect(show.legend=FALSE) + scale_fill_manual(values=c(ALT='#D80B8C', REF='gray80')) +
    coord_polar(theta="y") + xlim(c(1, 4)) + theme_void() +
    labs(title = paste0(round(x*100, 2), '%')) +
    theme(panel.background = element_rect(fill = "transparent", colour = NA),
          plot.background = element_rect(fill = "transparent", colour = NA),
          plot.title = element_text(hjust = 0.5, size=120))
  
  pietxt <- {print(pie); s()}
  invisible(dev.off())
  return(pietxt)
}

# make expression distribution plot as svg text
expr.dist <- function(gene){ # x = gene symbol
  s <- svglite::svgstring(); s()
  
  dat <- zscoredf %>%
    dplyr::select(-X) %>%
    filter(Gene == gene) %>%
    pivot_longer(cols = colnames(zscoredf)[3:ncol(zscoredf)], 
                 names_to='sampleName', values_to='zscore') %>%
    pivot_wider(names_from='Gene', values_from='zscore') %>%
    data.frame()
  colnames(dat) <- make.names(colnames(dat))

  df <- data.frame(
    sampleName=dat$sampleName, 
    zscore=dat[, make.names(gene)],
    geneName=gene,
    isPt=ifelse(grepl(rnaSampleID, dat$sampleName), 'yes', 'no')
    )
  
  plt <- ggplot(df) +
    aes(x=reorder(sampleName, zscore), y=zscore, fill=isPt, alpha=isPt, width=ifelse(isPt=='yes', 5, 1)) +
    geom_bar(stat='identity',position='identity', color='white', size=0, show.legend=FALSE) + 
    scale_fill_manual(values = c(yes='#00AEFF', no='gray80')) + 
    scale_alpha_manual(values = c(yes=1, no=0.5)) +
    labs(x = NULL, y = NULL, title = paste0('z-score = ', round(df$zscore[df$isPt=='yes'], 2)) ) +
    theme_void() +
    theme(panel.background = element_rect(fill = "transparent", colour = NA),
          plot.background = element_rect(fill = "transparent", colour = NA),
          plot.title = element_text(hjust = 0.5, size=46))

  plttxt <- {print(plt); s()}
  invisible(dev.off())
  return(plttxt)
}

# convert clonal tree json to matrix format
tree_mat <- function(z){
  b <- unlist(z)
  nm <- names(z)
  nx <- sapply(z, length)
  names(b) <- rep(nm, nx)
  c <- cbind(row.names(as.matrix(b)), as.matrix(b))
  row.names(c) <- NULL
  c <- apply(c, 2, as.integer)
  return(c)
}

get_clones <- function(populations){
  clone_df <- data.frame()
  for(clone in names(populations)){
    pop <- populations[[clone]]
    row <- data.frame(clone=clone, cellular_prevalence=pop$cellular_prevalence, num_cnvs=pop$num_cnvs, num_ssms=pop$num_ssms)
    clone_df <- rbind(clone_df, row)
    rm(pop, row)
  }
  return(clone_df)
}

list_variants <- function(alt, df){
    cols <- colnames(df)[grepl(paste0(alt, '.variant_names'), colnames(df), fixed=TRUE)]
    mat <- as.matrix(df[, c(cols)])
    
    apply(mat, 1, function(x){paste(as.character(na.omit(unique(as.character(unlist(strsplit(as.character(x), split=';')))))), collapse="<br> ")})
    }