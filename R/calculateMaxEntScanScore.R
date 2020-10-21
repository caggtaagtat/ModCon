## Function to calculate MaxEntScan score
calculateMaxEntScanScore <- function(seqVector, ssType) {

  #Convert seqVector to character in case its a factor
  seqVector <- as.character(seqVector)

  ## Error messages for wrong entries seq
  if(!ssType %in% c(3,5))
    stop(paste0("ERROR during setting of variable 'ssType'.",
    "The entered ss type must be numeric value 3 or 5, for SA or SD sites."))

  if(((ssType == 3) & (unique(nchar(seqVector)) != 23)))
    stop(paste0("ERROR during setting of variable 'seqVector'.",
    "3'ss (SA) sequences must be of length 23 for MaxEntScan score calculation."))

  if((ssType == 5) & (unique(nchar(seqVector)) != 9))
    stop(paste0("ERROR during setting of variable 'seqVector'.",
    "3'ss (SA) sequences must be of length 9 for MaxEntScan score calculation."))

  if (!all(strsplit(paste(seqVector, collapse = ""), "")[[1]] %in% c("a",
                                               "c", "g", "t", "G", "C", "T", "A")))
    warning(paste("WARNING during setting of variable 'seqVector'.",
    "One or more sequences showed characters apart from A C G and T.",
    " MaxEntScan score calculation not possible for the affected sequences."))


  ## Save current working directory
  curwd <- getwd()

  ## Go to ModCon package directory to access the required data
  fpath <- system.file("extdata", "me2x5", package = "ModCon")
  setwd(dirname(fpath))

  ## Delete inputSequences file in case it already exists
  if (file.exists("inputSequences"))
    file.remove("inputSequences")

  ## Create new text file with the sequences saved
  cat(seqVector, file = "inputSequences", sep = "\n", append = TRUE)

  ## Execute the respective Perl script with the respective Sequence file
  if (ssType == 3)
    cmd <- paste("score3.pl", "inputSequences")
  if (ssType == 5)
    cmd <- paste("score5.pl", "inputSequences")

  ## Save the calculated Maxent score file
  x <- system2("perl", args=cmd, stdout = TRUE)

  if (length(x) != 0) {

    ## Substracting the Scores from the Maxent Score Table
    x <- substr(x, (regexpr("\t", x)[[1]] + 1), nchar(x))

  }

  ## Go to back to original working directory
  setwd(curwd)

  ## Returning the maxent table
  return(x)

}


# use_data(Codons, hbg, hex, hex2, cds, overwrite = TRUE)


