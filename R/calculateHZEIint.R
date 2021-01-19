## HZEI integral calculation
calculateHZEIint <- function(ntSequence) {

  ## Error messages for wrong entries ntSequence
  if ((!all(strsplit(ntSequence, "")[[1]] %in% c("a", "c",
          "g", "t", "G", "C", "T", "A")))|(nchar(ntSequence)<11))
    stop("ERROR during setting of variable 'ntSequence'.",
    "The entered sequence must be a character string of A C G and T.")

  out <- getOverlappingVectorsFromVector(strsplit(toupper(ntSequence),
                                                  "")[[1]], 6, 5)
  out <- lapply(out, function(x) paste(x, collapse = ""))
  out <- out[nchar(unlist(out)) == 6]
  out <- as.character(hex$value[match(out, hex$seq)])
  out9 <- getOverlappingVectorsFromVector(out, 6, 5)
  out9 <- out9[lapply(out9, length) == 6]
  out9 <- lapply(out9, function(x) {
    mean(as.numeric(x))
  })
  out9 <- sum(unlist(out9))

  return(out9)
}

