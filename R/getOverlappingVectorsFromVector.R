## Overlap Vector Function
getOverlappingVectorsFromVector <- function(largeVector,
                                            subvectorLength,
                                            subvectorOverlap) {
  startPositions = seq(1, length(largeVector),
                       by = subvectorLength - subvectorOverlap)
  end_positions = startPositions + subvectorLength - 1
  end_positions[end_positions > length(largeVector)] = length(largeVector)
  lapply(seq_len(length(startPositions)), function(x) largeVector[startPositions[x]:
                                                                    end_positions[x]])
}
