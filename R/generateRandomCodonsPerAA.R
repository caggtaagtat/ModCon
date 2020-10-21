## Generate random sequence version
generateRandomCodonsPerAA <- function(aaVector) {

  ## Check if codonsInsert table is present
  if (!any(ls() == "codonsInsert"))
    codonsInsert <- Codons

  ## Paste together the codon sequences to form a new CDS
  paste(as.character(lapply(aaVector, function(x) {
    a <- c(codonsInsert$seq[codonsInsert$AA == x])
    a[sample(seq_len(length(a)), 1)]
  })), collapse = "")

}
