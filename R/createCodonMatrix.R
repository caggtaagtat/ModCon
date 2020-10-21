createCodonMatrix <- function(cds){

  ## Error messages for wrong entries cds
  if ((!all(strsplit(cds, "")[[1]] %in% c("a", "c",
        "g", "t", "G", "C", "T", "A")))|(nchar(cds)%%3!=0))
    stop(paste0("ERROR during setting of variable 'cds'.",
         "The entered sequence must be a character string of A C G and T",
         " with a length of a multiple by 3."))

  ## Create Codon Matrix by splitting up the sequence by 3nt
  sst <- strsplit(cds, "")[[1]]
  sst <- paste0(sst[c(TRUE, FALSE, FALSE)], sst[c(FALSE, TRUE, FALSE)],
                sst[c(FALSE, FALSE, TRUE)])

  ## Initiaing the Codons matrix plus corresponding amino acids
  fan <- matrix(nrow=2, ncol=nchar(cds)/3)
  fan[1, ] <- sst
  fan[2, ] <- sst

  ## Return the matrix
  return(fan)

}
