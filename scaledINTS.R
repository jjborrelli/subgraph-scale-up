load("C:/Users/jjborrelli/Desktop/motQSS.RData")
library(igraph)

conversion <- function(tm){
  for(i in 1:nrow(tm)){
    for(j in 1:ncol(tm)){
      if(tm[i,j] == 1 & tm[j,i] == 0){tm[j,i] <- -1}
    }
  }
  return(tm)
}

qssSCALED <- function(motMATS){
  mots <- length(motMATS)
  
  motifS <- motif_counter(lapply(motMATS, graph.adjacency))
  
  rem.l <- list()
  mchange <- list()
  for(i in 1:mots){
    mat <- motMATS[[i]]
    pairs1 <- unlist(apply(combn(1:nrow(mat), 2), 2, list), recursive = F)
    
    lmat <- lapply(pairs1,function(x){
      mat2 <- mat
      if(mat2[x[1],x[2]] == 1){mat2[x[1],x[2]] <- 0}
      if(mat2[x[2],x[1]] == 1){mat2[x[2],x[1]] <- 0};return(mat2)})
  
    removed <- sapply(lmat,sum) != sum(mat)
    lmat <- lmat[removed]
    
    rem.l[[i]] <- which(removed)
    mchange[[i]] <- abs(apply(do.call(rbind, apply(motif_counter(lapply(lmat, graph.adjacency)), 1, function(x){x - motifS[i,]})), c(1,2), 
                          function(x){if(x > 0){x <- 0}else{x}}))
  }
  
  
  class1 <- lapply(mchange, function(x){rowSums(x[,c(1,4,5)]) > 0 & rowSums(x[,-c(1,4,5)] == 0)})
  class2 <- lapply(mchange, function(x){x[,2] > 0 & rowSums(x[,-c(1,2,3,5)]) == 0})
  class3 <- lapply(mchange, function(x){rowSums(x[,-c(1,4,5)]) > 0})
  
  lam.max <- matrix(nrow = 1000, ncol = mots)
  for(n in 1:1000){
    matCon <- lapply(motMATS[1:mots], conversion)
    
    for(i in 1:mots){
      for(j in 1:length(rem.l[[i]])){
        row1 <- pairs1[[rem.l[[i]][j]]][1]
        col1 <- pairs1[[rem.l[[i]][j]]][2]
        
        if(class1[[i]][j]){if(matCon[[i]][row1, col1]==1){matCon[[i]][row1, col1] <- runif(1, 0, 10)}
                           if(matCon[[i]][col1, row1]==-1){matCon[[i]][col1, row1] <- runif(1, -1, 0)}}
        
        if(class2[[i]][j]){if(matCon[[i]][row1, col1]==1){matCon[[i]][row1, col1] <- runif(1, 0, 1)}
                           if(matCon[[i]][col1, row1]==-1){matCon[[i]][col1, row1] <- runif(1, -.1, 0)}}
        
        if(class3[[i]][j]){if(matCon[[i]][row1, col1]==1){matCon[[i]][row1, col1] <- runif(1, 0, .1)}
                           if(matCon[[i]][col1, row1]==-1){matCon[[i]][col1, row1] <- runif(1, -.01, 0)}else if(matCon[[i]][col1, row1]==1){
                             matCon[[i]][col1, row1] <- runif(1, 0, .1)}}
      }
    }
    
    matCon <- lapply(matCon, function(x){diag(x) <- runif(nrow(x), -1, 0); return(x)})
    
    lam.max[n,] <- sapply(matCon, function(x){max(Re(eigen(x)$values))})
  }
  
  qss.2 <- apply(lam.max, 2, function(x){sum(x < 0)/length(x)})
  
  return(qss.2)
}

system.time(
qss4.2 <- qssSCALED(fournode.am)
)
#112 sec

plot(qss4.2~qss4)
abline(a = 0, b = 1, xpd = FALSE)

system.time(
  qss5.2 <- qssSCALED(fivenode.am)
)
#625 sec

plot(qss5.2~qss5)
abline(a = 0, b = 1, xpd = FALSE)

system.time(
  qss10.2 <- qssSCALED(lapply(tens, get.adjacency, sparse = F))
)
#1228 sec

plot(qss10.2~qss10)
abline(a = 0, b = 1, xpd = FALSE)

NULL_qssSCALED <- function(motMATS){
  mots <- length(motMATS)
  
  motifS <- motif_counter(lapply(motMATS, graph.adjacency))
  
  rem.l <- list()
  mchange <- list()
  for(i in 1:mots){
    mat <- motMATS[[i]]
    pairs1 <- unlist(apply(combn(1:nrow(mat), 2), 2, list), recursive = F)
    
    lmat <- lapply(pairs1,function(x){
      mat2 <- mat

      if(mat2[x[1],x[2]] == 1){mat2[x[1],x[2]] <- 0}
      if(mat2[x[2],x[1]] == 1){mat2[x[2],x[1]] <- 0};return(mat2)})
    
    removed <- sapply(lmat,sum) != sum(mat)
    lmat <- lmat[removed]
    
    rem.l[[i]] <- which(removed)
    mchange[[i]] <- abs(apply(do.call(rbind, apply(motif_counter(lapply(lmat, graph.adjacency)), 1, function(x){x - motifS[i,]})), c(1,2), 
                              function(x){if(x > 0){x <- 0}else{x}}))
  }
  
  
  class1 <- lapply(mchange, function(x){rowSums(x[,c(1,4,5)]) > 0 & rowSums(x[,-c(1,4,5)]) == 0})
  class2 <- lapply(mchange, function(x){x[,2] > 0 & rowSums(x[,-c(1,2,3,5)]) == 0})
  class3 <- lapply(mchange, function(x){rowSums(x[,-c(1,4,5)]) > 0})
  
  matCon.N <- list()
  lam.max <- matrix(0, nrow = 1000, ncol = mots)
  lam.max.N <- matrix(0, nrow = 1000, ncol = mots)
  for(n in 1:1000){
    matCon <- lapply(motMATS[1:mots], conversion)
    
    for(i in 1:mots){
      for(j in 1:length(rem.l[[i]])){
        row1 <- pairs1[[rem.l[[i]][j]]][1]
        col1 <- pairs1[[rem.l[[i]][j]]][2]
        
        if(class1[[i]][j]){if(matCon[[i]][row1, col1]==1){matCon[[i]][row1, col1] <- runif(1, 0, 10)}
                           if(matCon[[i]][col1, row1]==-1){matCon[[i]][col1, row1] <- runif(1, -1, 0)}}
        
        if(class2[[i]][j]){if(matCon[[i]][row1, col1]==1){matCon[[i]][row1, col1] <- runif(1, 0, 1)}
                           if(matCon[[i]][col1, row1]==-1){matCon[[i]][col1, row1] <- runif(1, -.1, 0)}}
        
        if(class3[[i]][j]){if(matCon[[i]][row1, col1]==1){matCon[[i]][row1, col1] <- runif(1, 0, .1)}
                           if(matCon[[i]][col1, row1]==-1){matCon[[i]][col1, row1] <- runif(1, -.01, 0)}else if(matCon[[i]][col1, row1]==1){
                             matCon[[i]][col1, row1] <- runif(1, 0, .1)}}
      }
      links <- do.call(rbind, pairs1[rem.l[[i]]])
      same = TRUE
      while(same){
        s <- sample(1:nrow(links))
        same <- sum(s == 1:nrow(links)) == nrow(links)
      }
      
      matCon.N[[i]] <- matCon[[i]]
      matCon.N[[i]][links] <- matCon[[i]][links[s,]]
      matCon.N[[i]][links[,c(2,1)]] <- matCon[[i]][links[s,c(2,1)]]
    }
    
    matCon <- lapply(matCon, function(x){diag(x) <- runif(nrow(x), -1, 0); return(x)})
    matCon.N <- lapply(matCon.N, function(x){diag(x) <- runif(nrow(x), -1, 0); return(x)})
    
    lam.max[n,] <- sapply(matCon, function(x){max(Re(eigen(x)$values))})
    lam.max.N[n,] <- sapply(matCon.N, function(x){max(Re(eigen(x)$values))})
  }
  
  qss.2 <- apply(lam.max, 2, function(x){sum(x < 0)/length(x)})
  qss.2.N <- apply(lam.max.N, 2, function(x){sum(x < 0)/length(x)})
  
  return(list(qss.2, qss.2.N))
}

system.time(
  qss4.2.test <- NULL_qssSCALED(fournode.am)
)
#112 sec

plot(qss4.2~qss4)
abline(a = 0, b = 1, xpd = FALSE)

system.time(
  qss5.2 <- qssSCALED(fivenode.am)
)
#625 sec

plot(qss5.2~qss5)
abline(a = 0, b = 1, xpd = FALSE)

system.time(
  qss10.2 <- qssSCALED(lapply(tens, get.adjacency, sparse = F))
)
#1228 sec

plot(qss10.2~qss10)
abline(a = 0, b = 1, xpd = FALSE)



links <- do.call(rbind, pairs1[rem.l])
same = TRUE
while(same){
  s <- sample(1:nrow(links))
  same <- sum(s == 1:nrow(links)) == nrow(links)
}


matCon <- conversion(fivenode.am[[1]])
matCon <- apply(matCon, c(1,2), function(x){if(x == 1){abs(rnorm(1))}else if(x == -1){-abs(rnorm(1))}else{0}})

matCon.N <- matCon
matCon.N[links] <- matCon[links[s,]]
matCon.N[links[,c(2,1)]] <- matCon[links[s,c(2,1)]]
matCon.N == matCon

