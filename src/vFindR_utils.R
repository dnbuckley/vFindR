require(Rsamtools)
require(rtracklayer)
require(parallel)
require(dplyr)
require(openxlsx)


# A function for collapsing the list of lists into a single list
# as per the Rsamtools vignette
.unlist <- function (x){
  x1 <- x[[1L]]
  if (is.factor(x1)){
    structure(unlist(x), class = "factor", levels = levels(x1))
  } else {
    do.call(c, x)
  }
}
.readBAM <- function(bamFile){
  bam <- scanBam(bamFile)
  bam_field <- names(bam[[1]])
  list <- lapply(bam_field, function(y) .unlist(lapply(bam, "[[", y)))
  bam_df <- do.call("DataFrame", list)
  names(bam_df) <- bam_field
  #return a list that can be called as a data frame
  return(as.data.frame(bam_df))
}
.bamDF2GR <- function(df) {
  GRanges(seqnames = df$rname, ranges = IRanges(start = df$pos,
                                                end = df$pos))
}
.vChr2Species <- function(vchrs, ref = NULL) {
  if (is.null(ref)) {
    ref <- read.xlsx("~/Desktop/salhia_lab/vFindR/viral_seqnames.xlsx")
  }
  x <- data.frame(seqnames = vchrs)
  x <- x %>% 
    dplyr::left_join(ref, by = "seqnames")
  return(x$species)
} 
.species2vChr <- function(species, ref = NULL) {
  if (is.null(ref)) {
    ref <- read.xlsx("~/Desktop/salhia_lab/vFindR/viral_seqnames.xlsx")
  }
  x <- data.frame(species = species)
  x <- x %>% 
    dplyr::left_join(ref, by = "species")
  return(x$seqnames)
}












