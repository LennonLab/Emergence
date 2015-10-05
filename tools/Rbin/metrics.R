#### A function to generate observed richness
S.obs <- function(x = ""){ rowSums(x > 0) * 1}


sp.turnover <- function(site1, site2){
  
  site1 <- site1[site1 != '']
  site1 <- site1[!is.na(site1)]
  site2 <- site2[site2 != '']
  site2 <- site2[!is.na(site2)]
  #if(length(site1) | length(site2) == 0){
  #  return -1
  #  }
  gamma = union(site1, site2)         # Gamma species pool
  s     = length(gamma)                                   # Gamma richness
  a.bar = mean(c(length(site1), length(site2)))   # Mean sample richness
  b.w   = round(s/a.bar - 1, 4)
  return(b.w)
  }



get.vectors <- function(inputFile){
  vectors <- c()
  i <- 1
  con  <- file(inputFile, open = "r")
  while (length(oneLine <- readLines(con, n = 1)) > 0) {
    myLine <- unlist((strsplit(oneLine, ",")))
    myLine[myLine == '[]'] <- -1
    myLine[is.na(myLine)] <- -1
    myLine <- myLine[which(myLine!=-1)]
    myLine <- myLine[which(myLine!="")]
  
    if (length(myLine > 0)){
      myLine <- as.vector(myLine)
      vectors[i] <- myLine
      i <- i + 1
    }
  }
  close(con)
  return(vectors)
}


sem <- function(x) return(sqrt(var(x)/length(x)))
