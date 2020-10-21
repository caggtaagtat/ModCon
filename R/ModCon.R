## ModCon function
ModCon <- function(cds, sdSeqStartPosition, upChangeCodonsIn=16,
                   downChangeCodonsIn=16, modconMode="optimalContext",
                   sdMaximalHBS=10, acMaximalMaxent=4, optiRate=100,
                   nGenerations=30, parentSize=300, startParentSize=1000,
                   bestRate=40, semiLuckyRate=20, luckyRate=5,
                   mutationRate=1e-04, nCores=-1) {

  ## Error messages for wrong entries cds
  if (!all(strsplit(cds, "")[[1]] %in% c("a", "c", "g", "t", "G", "C", "T", "A")))
    stop(paste0("ERROR during setting of variable 'cds'. The entered sequence",
                " must be a character string of A C G and T."))

  ## sdSeqStartPosition
  if (sdSeqStartPosition < 1 | sdSeqStartPosition > (nchar(cds) - 10)) {
    stop(paste0("ERROR during setting of variable 'sdSeqStartPosition'. SD",
                " sequence not within cds."))
  }

  ## upChangeCodonsIn
  if (!is.numeric(upChangeCodonsIn) | (upChangeCodonsIn < 4)) {
    stop(paste0("ERROR during setting of variable 'upChangeCodonsIn'. Must",
                " be numeric and greater 4."))
  }

  ## downChangeCodonsIn
  if (!is.numeric(downChangeCodonsIn) | (downChangeCodonsIn < 4)) {
    stop(paste0("ERROR during setting of variable 'downChangeCodonsIn'. Must",
                " be an numeric and greater 4."))
  }

  ## modconMode
  if ((!modconMode %in% c("optimalContext", "suboptimalContext"))) {
    stop(paste0("Setting of the variable 'modconMode' not correct. Please",
                " either choose the value 'optimalContext' or 'suboptimalContext'."))
  }

  ## sdMaximalHBS
  if (sdMaximalHBS < 1.8 | sdMaximalHBS > 23.8) {
    stop(paste0("ERROR during setting of variable 'sdMaximalHBS'. HBS limit",
                " must be within 1.8 to 23.8."))
  }

  ## sdMaximalHBS
  if (!is.numeric(sdMaximalHBS)) {
    stop(paste0("ERROR during setting of variable 'acMaximalMaxent'.",
                " MaxEntScan score limit must be numeric."))
  }

  ## optiRate
  if (!is.numeric(optiRate) | optiRate < 0) {
    stop(paste0("ERROR during setting of variable 'optiRate'. optiRate",
                " must be numeric and not lower than 0 ."))
  }

  ## nGenerations
  if (!is.numeric(nGenerations) | nGenerations < 2) {
    stop(paste0("ERROR during setting of variable 'nGenerations'. nGenerations",
                " must be numeric and greater 1."))
  }

  ## parentSize
  if (!is.numeric(parentSize) | parentSize < 10) {
    stop(paste0("ERROR during setting of variable 'parentSize'. parentSize",
                " must be numeric and greater 9."))
  }

  ## startParentSize
  if (!is.numeric(startParentSize) | startParentSize < 10) {
    stop(paste0("ERROR during setting of variable 'startParentSize'. startParentSize",
                " must be numeric and greater 9."))
  }

  ## bestRate
  if (!is.numeric(bestRate) | bestRate <= 0) {
    stop(paste0("ERROR during setting of variable 'bestRate'. bestRate must",
                " be numeric and greater 0."))
  }

  ## semiLuckyRate
  if (!is.numeric(semiLuckyRate) | semiLuckyRate <= 0) {
    stop(paste0("ERROR during setting of variable 'semiLuckyRate'. semiLuckyRate",
                " must be numeric and greater 0."))
  }

  ## luckyRate
  if (!is.numeric(luckyRate) | luckyRate <= 0) {
    stop(paste0("ERROR during setting of variable 'luckyRate'. luckyRate",
                " must be numeric and greater 0."))
  }

  ## mutationRate
  if (!is.numeric(mutationRate) | mutationRate < 0 | mutationRate > 1) {
    stop(paste0("ERROR during setting of variable 'mutationRate'. mutationRate",
                " must be numeric and range from 0 to 1."))
  }

  ## Define region of interest
  upChangeCodons <- (upChangeCodonsIn * 3) - 1
  downChangeCodons <- (downChangeCodonsIn * 3)

  ## Define codonsInsert object, to allow differences between lookup table Codons
  ## and insertion table codonsInsert
  codonsInsert <- Codons

  ## If SD context will be optimized to enhance usage
  if (modconMode == "optimalContext") {

    ## Find position of first nucleotid in first codon upstream of SD
    endPos <- sdSeqStartPosition
    while (endPos%%3 != 0) endPos <- endPos - 1

    ## Select subsequence to optimize upstream of the SD (within ORF)
    inputSeq <- substr(cds, endPos - upChangeCodons, endPos)
    insertStartPos <- endPos - upChangeCodons

    ## Report next step
    message("Maximizing HZEI integral of upstream sequence...")

    ## Maximize HZEI of upstream sequence
    res <- changeSequenceHZEI(inputSeq, increaseHZEI=TRUE, nGenerations=nGenerations,
                              parentSize=parentSize, startParentSize=startParentSize,
                              bestRate=bestRate, semiLuckyRate=semiLuckyRate,
                              luckyRate=luckyRate, mutationRate=mutationRate,
                              optiRate=optiRate, sdMaximalHBS=sdMaximalHBS,
                              acMaximalMaxent=acMaximalMaxent, nCores=nCores)

    ## Saving wildtype and optimized sequence
    before <- inputSeq
    after <- res[[3]]
    upstreamPositiveSeq <- after

    ## Get downstream sequence for adjustment
    endPos <- sdSeqStartPosition + 10
    while (endPos%%3 != 0) endPos <- endPos + 1
    inputSeq <- substr(cds, endPos + 1, endPos + downChangeCodons)
    insertEndPos <- endPos + downChangeCodons

    ## Report next step
    message("Minimizing HZEI integral of downstream sequence ...")

    ## Minimize HZEI of downstream sequence
    res <- changeSequenceHZEI(inputSeq, increaseHZEI=FALSE, nGenerations=nGenerations,
                              parentSize=parentSize, startParentSize=startParentSize,
                              bestRate=bestRate, semiLuckyRate=semiLuckyRate,
                              luckyRate=luckyRate, mutationRate=mutationRate,
                              optiRate=optiRate, sdMaximalHBS=sdMaximalHBS,
                              acMaximalMaxent=acMaximalMaxent, nCores=nCores)

    ## Saving wildtype and optimized sequence
    before <- inputSeq
    after <- res[[3]]
    downstreamPositiveSeq <- after

    ## Get SD sequence
    sdSeq <- substr(cds, (insertStartPos + upChangeCodons + 1), (insertEndPos -
                                                                   (downChangeCodons)))

    # Fuse all sequences to final sequence
    completeSequence <- paste0(upstreamPositiveSeq, sdSeq, downstreamPositiveSeq)
    completeSequence <- strsplit(completeSequence, "")[[1]]
    cdsInput <- cds
    cds <- strsplit(cds, "")[[1]]
    cds[insertStartPos:insertEndPos] <- completeSequence
    cds <- paste(cds, collapse = "")

    ## Report important information about adjustment
    message(paste("Position of first nt of SD seq:", sdSeqStartPosition))
    message(paste("Length of substituted sequences upstream and dowstream:",
                  (upChangeCodons + 1), "and", downChangeCodons, "nt"))

    ## If SD context will be optimized to silence usage
  } else {

    ## Find position of first nucleotid in first codon upstream of SD
    endPos <- sdSeqStartPosition + 10
    while (endPos%%3 != 0) endPos <- endPos + 1

    ## Select subsequence to optimize upstream of the SD (within ORF)
    inputSeq <- substr(cds, endPos + 1, endPos + downChangeCodons)
    insertEndPos <- endPos + downChangeCodons

    ## Report next step
    message("Maximize HZEI integral of downstream sequence...")

    ## Maximize HZEI of downstream sequence
    res <- changeSequenceHZEI(inputSeq, increaseHZEI=TRUE, nGenerations=nGenerations,
                              parentSize=parentSize, startParentSize=startParentSize,
                              bestRate=bestRate, semiLuckyRate=semiLuckyRate,
                              luckyRate=luckyRate, mutationRate=mutationRate,
                              optiRate=optiRate, sdMaximalHBS=sdMaximalHBS,
                              acMaximalMaxent=acMaximalMaxent, nCores=nCores)

    ## Saving wildtype and optimized sequence
    before <- inputSeq
    after <- res[[3]]
    downstreamPositiveSeq <- after

    ## Get downstream sequence for adjustment
    endPos <- sdSeqStartPosition
    while (endPos%%3 != 0) endPos <- endPos - 1
    inputSeq <- substr(cds, endPos - upChangeCodons, endPos)
    insertStartPos <- endPos - upChangeCodons

    ## Report next step
    message("Minimizing HZEI integral of upstream sequence...")

    ## Minimize HZEI of downstream sequence
    res <- changeSequenceHZEI(inputSeq, increaseHZEI=FALSE, nGenerations=nGenerations,
                              parentSize=parentSize, startParentSize=startParentSize,
                              bestRate=bestRate, semiLuckyRate=semiLuckyRate,
                              luckyRate=luckyRate, mutationRate=mutationRate,
                              optiRate=optiRate, sdMaximalHBS=sdMaximalHBS,
                              acMaximalMaxent=acMaximalMaxent, nCores=nCores)

    ## Saving wildtype and optimized sequence
    before <- inputSeq
    after <- res[[3]]
    upstreamPositiveSeq <- after

    ## Get SD sequence
    sdSeq <- substr(cds, (insertStartPos + upChangeCodons + 1), (insertEndPos -
                                                                   (downChangeCodons)))

    # Fuse all sequences to final sequence
    completeSequence <- paste0(upstreamPositiveSeq, sdSeq, downstreamPositiveSeq)
    completeSequence <- strsplit(completeSequence, "")[[1]]
    cdsInput <- cds
    cds <- strsplit(cds, "")[[1]]
    cds[insertStartPos:insertEndPos] <- completeSequence
    cds <- paste(cds, collapse = "")

    ## Report important information about adjustment
    message(paste("Position of first nt of SD seq:", sdSeqStartPosition))
    message(paste("Length of substituted sequences upstream and dowstream:",
                  (upChangeCodons + 1), "and", downChangeCodons, "nt"))

  }

  ## Return full cds with the SD and its adjusted sequence context
  return(cds)
}

