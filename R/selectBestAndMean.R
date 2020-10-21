## Select maximal and mean HZEI integral wthin a sequenceVector
selectBestAndMean <- function(sequenceVector,
                              clusterName,
                              increaseHZEI=TRUE) {

  sequenceVectorDf <- data.frame(seq=as.character(sequenceVector),
                                 int=0)
  sequenceVectorDf$int <- as.numeric(parLapply(clusterName,
                                               as.character(sequenceVectorDf$seq),
                                               calculateHZEIint))

  if (increaseHZEI)
    return(c(max(sequenceVectorDf$int), mean(sequenceVectorDf$int)))
  if (!increaseHZEI)
    return(c(min(sequenceVectorDf$int), mean(sequenceVectorDf$int)))
}
