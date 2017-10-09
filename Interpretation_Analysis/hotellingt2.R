calchot <- function(pcaobj,samplenum) {
  return(sum(pcaobj$x[samplenum,]/pcaobj$sdev)^2)
}