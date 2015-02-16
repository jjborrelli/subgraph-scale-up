library(igraph)

# motif_counter function from "web_functions.R"
motif_counter <- function(graph.lists){
  require(igraph)
  
  if(!is.list(graph.lists)){
    stop("The input should be a list of graph objects")
  }
  
  triad.count <- lapply(graph.lists, triad.census)
  triad.matrix <- matrix(unlist(triad.count), nrow = length(graph.lists), ncol = 16, byrow = T)
  colnames(triad.matrix) <- c("empty", "single", "mutual", "s5", "s4", "s1", "d4",
                              "d3", "s2", "s3","d8", "d2", "d1", "d5", "d7", "d6")
  
  triad.df <- as.data.frame(triad.matrix)
  
  motif.data.frame <- data.frame(s1 = triad.df$s1, s2 = triad.df$s2, s3 = triad.df$s3, s4 = triad.df$s4, 
                                 s5 = triad.df$s5, d1 = triad.df$d1, d2 = triad.df$d2, d3 = triad.df$d3, d4 = triad.df$d4,
                                 d5 = triad.df$d5, d6 = triad.df$d6, d7 = triad.df$d7, d8 = triad.df$d8)
  
  return(motif.data.frame)
}

# Stability analysis functions from "motifAnalysis.R"
conversion <- function(tm){
  for(i in 1:nrow(tm)){
    for(j in 1:ncol(tm)){
      if(tm[i,j] == 1 & tm[j,i] == 0){tm[j,i] <- -1}
    }
  }
  return(tm)
}

ran.unif <- function(motmat){
  newmat <- apply(motmat, c(1,2), function(x){
    if(x==1){runif(1, 0, 10)}else if(x==-1){runif(1, -1, 0)} else{0}
  })
  diag(newmat) <- runif(nrow(newmat), -1, 0)
  return(newmat)
}

maxRE <- function(rmat){
  lam.max <- max(Re(eigen(rmat)$values))
  return(lam.max)
}

eig.analysis <- function(n, matrices){
  cols <- length(matrices)
  rows <- n
  eigenMATRIX <- matrix(0, nrow = rows, ncol = cols)
  for(i in 1:n){
    ranmat <- lapply(matrices, ran.unif)
    
    eigs <- sapply(ranmat, maxRE)
    eigenMATRIX[i,] <- eigs
  }
  return(eigenMATRIX)
}



# Motif Dictionary
fournode <- list(
  id14 = matrix(c(2,1,3,1,4,1), ncol = 2, byrow = T),
  id28 = matrix(c(3,1,4,1,1,2), ncol = 2, byrow = T),
  id30 = matrix(c(3,1,4,1,1,2,2,1), ncol = 2, byrow = T),
  id74 = matrix(c(2,1,3,2,4,1), ncol = 2, byrow = T),
  id76 = matrix(c(3,1,3,2,4,1), ncol = 2, byrow = T),
  id78 = matrix(c(2,1,3,1,4,1,3,2), ncol = 2, byrow = T),
  id90 = matrix(c(2,1,1,2,3,2,4,1), ncol = 2, byrow = T),
  id92 = matrix(c(1,2,3,1,3,2,4,1), ncol = 2, byrow = T), 
  id94 = matrix(c(1,2,2,1,4,1,3,1,3,2), ncol = 2, byrow = T),
  id204 = matrix(c(3,1,3,2,4,1,4,2), ncol = 2, byrow = T),
  id206 = matrix(c(2,1,3,1,4,1,3,2,4,2), ncol = 2, byrow = T),
  id222 = matrix(c(2,1,3,1,4,1,3,2,4,2,1,2), ncol = 2, byrow = T),
  id280 = matrix(c(1,2,1,3,4,1), ncol = 2, byrow = T),
  id282 = matrix(c(1,2,1,3,4,1,2,1), ncol = 2, byrow = T),
  id286 = matrix(c(1,2,1,3,4,1,2,1,3,1), ncol = 2, byrow = T),
  id328 = matrix(c(1,3,4,1,3,2), ncol = 2, byrow = T),
  id330 = matrix(c(2,1,1,3,4,1,3,2), ncol = 2, byrow = T),
  id332 = matrix(c(4,1,3,2,3,1,1,3), ncol = 2, byrow = T),
  id334 = matrix(c(4,1,3,2,3,1,1,3,2,1), ncol = 2, byrow = T),
  id344 = matrix(c(4,1,3,2,3,1,1,3,1,2), ncol = 2, byrow = T),
  id346 = matrix(c(2,1,1,3,4,1,3,2,1,2), ncol = 2, byrow = T),
  id348 = matrix(c(4,1,3,2,3,1,1,3,1,2), ncol = 2, byrow = T),
  id350 = matrix(c(4,1,3,2,3,1,1,3,2,1,1,2), ncol = 2, byrow = T),
  id390 = matrix(c(2,1,1,3,3,1,4,2), ncol = 2, byrow = T),
  id392 = matrix(c(1,3,4,1,4,2), ncol = 2, byrow = T),
  id394 = matrix(c(2,1,1,3,4,1,4,2), ncol = 2, byrow = T),
  id396 = matrix(c(3,1,1,3,4,1,4,2), ncol = 2, byrow = T),
  id398 = matrix(c(2,1,1,3,3,1,4,2,4,1), ncol = 2, byrow = T),
  id404 = matrix(c(1,3,3,1,1,2,4,2), ncol = 2, byrow = T),
  id406 = matrix(c(1,3,3,1,1,2,2,1,4,2), ncol = 2, byrow = T),
  id408 = matrix(c(1,3,2,1,4,1,4,2), ncol = 2, byrow = T),
  id410 = matrix(c(2,1,1,3,4,2,4,1,1,2), ncol = 2, byrow = T),
  id412 = matrix(c(1,2,3,1,1,3,4,1,4,2), ncol = 2, byrow = T),
  id414 = matrix(c(1,2,3,1,1,3,4,1,4,2,2,1), ncol = 2, byrow = T),
  id454 = matrix(c(2,1,1,3,3,1,3,2,4,2), ncol = 2, byrow = T),
  id456 = matrix(c(1,3,4,1,2,3,4,2), ncol = 2, byrow = T),
  id458 = matrix(c(1,3,2,1,4,1,4,2,3,2), ncol = 2, byrow = T),
  id460 = matrix(c(1,3,3,1,4,1,4,2,3,2), ncol = 2, byrow = T),
  id462 = matrix(c(2,1,1,3,3,1,3,2,4,1,4,2), ncol = 2, byrow = T),
  id468 = matrix(c(1,2,1,3,3,1,3,2,4,2), ncol = 2, byrow = T),
  id470 = matrix(c(1,3,3,1,1,2,2,1,3,2,4,1,4,2), ncol = 2, byrow = T),
  id472 = matrix(c(1,3,1,2,4,1,4,2,3,2), ncol = 2, byrow = T),
  id474 = matrix(c(1,3,1,2,2,1,3,2,4,1,4,2), ncol = 2, byrow = T),
  id476 = matrix(c(1,3,3,1,1,2,3,2,4,1,4,2), ncol = 2, byrow = T),
  id478 = matrix(c(1,3,3,1,1,2,2,1,3,2,4,1,4,2), ncol = 2, byrow = T),
  id856 = matrix(c(1,3,1,2,2,3,3,2,4,1), ncol = 2, byrow = T),
  id858 = matrix(c(1,3,1,2,2,1,2,3,3,2,4,1), ncol = 2, byrow = T),
  id862 = matrix(c(1,3,3,1,1,2,2,1,2,3,3,2,4,1), ncol = 2, byrow = T),
  id904 = matrix(c(1,3,2,3,4,1,4,2), ncol = 2, byrow = T),
  id906 = matrix(c(1,3,2,1,2,3,4,1,4,2), ncol = 2, byrow = T),
  id908 = matrix(c(1,3,3,1,4,1,4,2,2,3), ncol = 2, byrow = T),
  id910 = matrix(c(1,3,3,1,2,1,2,3,4,1,4,2), ncol = 2, byrow = T),
  id922 = matrix(c(1,3,1,2,2,1,2,3,4,1,4,2), ncol = 2, byrow = T),
  id924 = matrix(c(1,3,3,1,1,2,2,3,4,1,4,2), ncol = 2, byrow = T),
  id926 = matrix(c(1,3,3,1,1,2,2,1,4,1,4,2,2,3), ncol = 2, byrow = T),
  id972 = matrix(c(1,3,3,1,2,3,3,2,4,1,4,2), ncol = 2, byrow = T),
  id974 = matrix(c(1,3,3,1,2,1,2,3,3,2,4,1,4,2), ncol = 2, byrow = T),
  id990 = matrix(c(1,3,3,1,1,2,2,1,2,3,3,2,4,1,4,2), ncol = 2, byrow = T),
  id2184 = matrix(c(4,1,4,2,4,3), ncol = 2, byrow = T),
  id2186 = matrix(c(4,1,4,2,4,3,2,1), ncol = 2, byrow = T),
  id2190 = matrix(c(4,1,4,2,4,3,2,1,3,1), ncol = 2, byrow = T),
  id2202 = matrix(c(1,2,2,1,4,1,4,2,4,3), ncol = 2, byrow = T),
  id2204 = matrix(c(3,1,1,2,4,1,4,2,4,3), ncol = 2, byrow = T),
  id2206 = matrix(c(1,2,2,1,3,1,4,1,4,2,4,3), ncol = 2, byrow = T),
  id2252 = matrix(c(3,1,3,2,4,1,4,2,4,3), ncol = 2, byrow = T),
  id2254 = matrix(c(3,1,2,1,3,2,4,1,4,2,4,3), ncol = 2, byrow = T),
  id2270 = matrix(c(1,2,2,1,3,1,3,2,4,1,4,2,4,3), ncol = 2, byrow = T),
  id2458 = matrix(c(1,2,2,1,1,3,4,1,4,2,4,3), ncol = 2, byrow = T),
  id2462 = matrix(c(1,3,3,1,1,2,2,1,4,1,4,2,4,3), ncol = 2, byrow = T),
  id2506 = matrix(c(1,3,2,1,3,2,4,1,4,2,4,3), ncol = 2, byrow = T),
  id2510 = matrix(c(1,3,3,1,2,1,3,2,4,1,4,2,4,3), ncol = 2, byrow = T),
  id2524 = matrix(c(1,2,1,3,3,1,3,2,4,1,4,2,4,3), ncol = 2, byrow = T),
  id2526 = matrix(c(1,2,2,1,1,3,3,1,3,2,4,1,4,2,4,3), ncol = 2, byrow = T),
  id3038 = matrix(c(1,2,2,1,1,3,3,1,3,2,2,3,4,1,4,2,4,3), ncol = 2, byrow = T),
  id4370 = matrix(c(1,3,1,2,2,1,4,1), ncol = 2, byrow = T),
  id4374 = matrix(c(1,3,3,1,1,2,2,1,4,1), ncol = 2, byrow = T),
  id4382 = matrix(c(1,2,2,1,1,3,3,1,1,4,4,1), ncol = 2, byrow = T),
  id4418 = matrix(c(1,3,2,1,1,4,3,2), ncol = 2, byrow = T),
  id4420 = matrix(c(1,3,3,1,3,2,1,4), ncol = 2, byrow = T),
  id4422 = matrix(c(1,3,3,1,3,2,2,1,1,4), ncol = 2, byrow = T),
  id4424 = matrix(c(1,3,3,2,4,1,1,4), ncol = 2, byrow = T),
  id4426 = matrix(c(1,3,2,1,3,2,1,4,4,1), ncol = 2, byrow = T),
  id4428 = matrix(c(1,3,3,1,1,4,4,1,3,2), ncol = 2, byrow = T),
  id4430 = matrix(c(1,3,3,1,2,1,1,4,4,1,3,2), ncol = 2, byrow = T),
  id4434 = matrix(c(3,1,1,2,2,1,1,4,3,2), ncol = 2, byrow = T),
  id4436 = matrix(c(1,3,3,1,1,2,1,4,3,2), ncol = 2, byrow = T),
  id4438 = matrix(c(1,2,2,1,1,3,3,1,3,2,1,4), ncol = 2, byrow = T),
  id4440 = matrix(c(1,3,1,2,1,4,4,1,3,2), ncol = 2, byrow = T),
  id4442 = matrix(c(1,3,1,2,2,1,1,4,4,1,3,2), ncol = 2, byrow = T),
  id4444 = matrix(c(1,3,3,1,1,4,4,1,1,2,3,2), ncol = 2, byrow = T),
  id4446 = matrix(c(1,3,3,1,1,4,4,1,1,2,2,1,3,2), ncol = 2, byrow = T),
  id4546 = matrix(c(1,3,2,1,1,4,4,2,3,2), ncol = 2, byrow = T),
  id4548 = matrix(c(1,3,3,1,1,4,3,2,4,2), ncol = 2, byrow = T),
  id4550 = matrix(c(1,3,3,1,2,1,1,4,3,2,4,2), ncol = 2, byrow = T),
  id4556 = matrix(c(1,3,3,1,1,4,3,2,4,2,4,1), ncol = 2, byrow = T),
  id4558 = matrix(c(1,3,3,1,1,4,4,1,4,2,2,1,3,2), ncol = 2, byrow = T),
  id4562 = matrix(c(1,3,1,2,2,1,3,2,4,2,1,4), ncol = 2, byrow = T),
  id4564 = matrix(c(1,3,3,1,1,2,1,4,3,2,4,2), ncol = 2, byrow = T),
  id4566 = matrix(c(1,3,3,1,1,2,2,1,3,2,4,2,1,4), ncol = 2, byrow = T),
  id4572 = matrix(c(1,3,3,1,1,2,3,2,4,2,1,4,4,1), ncol = 2, byrow = T),
  id4574 = matrix(c(1,2,2,1,1,3,3,1,1,4,4,1,3,2,4,2), ncol = 2, byrow = T),
  id4678 = matrix(c(3,1,2,1,1,4,3,2,2,3), ncol = 2, byrow = T),
  id4682 = matrix(c(2,1,1,4,4,1,2,3,3,2), ncol = 2, byrow = T),
  id4686 = matrix(c(3,1,2,1,1,4,4,1,2,3,3,2), ncol = 2, byrow = T),
  id4692 = matrix(c(3,1,1,2,1,4,2,3,3,2), ncol = 2, byrow = T),
  id4694 = matrix(c(1,2,2,1,1,4,2,3,3,2,3,1), ncol = 2, byrow = T),
  id4698 = matrix(c(1,2,2,1,2,3,3,2,1,4,4,1), ncol = 2, byrow = T),
  id4700 = matrix(c(3,1,1,2,1,4,4,1,2,3,3,2), ncol = 2, byrow = T),
  id4702 = matrix(c(3,1,1,2,2,1,1,4,4,1,2,3,3,2), ncol = 2, byrow = T),
  id4740 = matrix(c(3,1,4,2,1,4,2,3), ncol = 2, byrow = T),
  id4742 = matrix(c(3,1,4,2,1,4,2,3,2,1), ncol = 2, byrow = T),
  id4748 = matrix(c(3,1,4,2,1,4,2,3,4,1), ncol = 2, byrow = T),
  id4750 = matrix(c(3,1,2,1,1,4,4,1,2,3,4,2), ncol = 2, byrow = T),
  id4758 = matrix(c(1,2,2,1,3,1,4,2,2,3,1,4), ncol = 2, byrow = T),
  id4764 = matrix(c(1,2,3,1,2,3,1,4,4,1,4,2), ncol = 2, byrow = T),
  id4766 = matrix(c(1,2,2,1,1,4,4,1,3,1,2,3,4,2), ncol = 2, byrow = T),
  id4812 = matrix(c(3,1,4,2,1,4,4,1,2,3,3,2), ncol = 2, byrow = T),
  id4814 = matrix(c(3,1,4,2,1,4,4,1,2,3,3,2,2,1), ncol = 2, byrow = T),
  id4830 = matrix(c(3,1,4,2,1,4,4,1,2,3,3,2,2,1,1,2), ncol = 2, byrow = T),
  id4946 = matrix(c(1,2,2,1,1,3,1,4,2,3,3,2), ncol = 2, byrow = T),
  id4950 = matrix(c(1,2,2,1,1,3,3,1,2,3,3,2,1,4), ncol = 2, byrow = T),
  id4952 = matrix(c(1,2,1,3,1,4,4,1,2,3,3,2), ncol = 2, byrow = T),
  id4954 = matrix(c(1,2,1,3,1,4,4,1,2,3,3,2,2,1), ncol = 2, byrow = T),
  id4958 = matrix(c(1,2,1,3,1,4,4,1,2,3,3,2,2,1,3,1), ncol = 2, byrow = T),
  id4994 = matrix(c(2,1,1,3,1,4,2,3,4,2), ncol = 2, byrow = T),
  id4998 = matrix(c(2,1,1,3,1,4,2,3,4,2,3,1), ncol = 2, byrow = T),
  id5002 = matrix(c(2,1,1,3,1,4,2,3,4,2,4,1), ncol = 2, byrow = T),
  id5004 = matrix(c(1,3,3,1,1,4,4,1,2,3,4,2), ncol = 2, byrow = T),
  id5006 = matrix(c(1,3,3,1,2,1,1,4,4,1,2,3,4,2), ncol = 2, byrow = T),
  id5010 = matrix(c(1,3,1,2,2,1,1,4,2,3,4,2), ncol = 2, byrow = T),
  id5012 = matrix(c(1,3,3,1,1,2,1,4,2,3,4,2), ncol = 2, byrow = T),
  id5014 = matrix(c(1,3,3,1,1,2,2,1,1,4,2,3,4,2), ncol = 2, byrow = T),
  id5016 = matrix(c(1,3,1,2,1,4,4,1,2,3,4,2), ncol = 2, byrow = T),
  id5018 = matrix(c(1,2,2,1,1,3,1,4,4,1,2,3,4,2), ncol = 2, byrow = T),
  id5020 = matrix(c(1,3,3,1,1,4,4,1,1,2,2,3,4,2), ncol = 2, byrow = T),
  id5022 = matrix(c(1,2,2,1,1,3,3,1,1,4,4,1,2,3,4,2), ncol = 2, byrow = T),
  id5058 = matrix(c(1,3,2,1,1,4,2,3,3,2,4,2), ncol = 2, byrow = T),
  id5062 = matrix(c(1,3,3,1,2,1,1,4,2,3,3,2,4,2), ncol = 2, byrow = T),
  id5064 = matrix(c(1,3,1,4,2,3,3,2,4,1,4,2), ncol = 2, byrow = T),
  id5066 = matrix(c(1,3,1,4,4,1,2,3,3,2,2,1,4,2), ncol = 2, byrow = T),
  id5068 = matrix(c(1,3,3,1,1,4,4,1,2,3,3,2,4,2), ncol = 2, byrow = T),
  id5070 = matrix(c(1,3,3,1,1,4,4,1,2,3,3,2,4,2,2,1), ncol = 2, byrow = T),
  id5074 = matrix(c(1,3,1,4,4,1,2,3,3,2,4,2,2,1,1,2), ncol = 2, byrow = T),
  id5076 = matrix(c(1,3,3,1,2,3,3,2,1,4,1,2,4,2), ncol = 2, byrow = T),
  id5078 = matrix(c(1,3,3,1,1,2,2,1,2,3,3,2,4,2,1,4), ncol = 2, byrow = T),
  id5080 = matrix(c(1,3,1,2,2,3,3,2,4,2,1,4,4,1), ncol = 2, byrow = T),
  id5082 = matrix(c(1,3,1,2,2,3,3,2,4,2,1,4,4,1,2,1), ncol = 2, byrow = T),
  id5084 = matrix(c(1,3,1,2,2,3,3,2,4,2,1,4,4,1,3,1), ncol = 2, byrow = T),
  is5086 = matrix(c(1,3,3,1,1,2,2,1,2,3,3,2,4,2,1,4,4,1), ncol = 2, byrow = T),
  id6342 = matrix(c(2,1,3,1,1,4,3,2,4,3,4,2), ncol = 2, byrow = T),
  id6348 = matrix(c(3,1,3,2,1,4,4,1,4,3,4,2), ncol = 2, byrow = T),
  id6350 = matrix(c(2,1,1,4,4,1,4,2,4,3,3,2,3,1), ncol = 2, byrow = T),
  id6356 = matrix(c(1,2,1,4,4,2,4,3,3,2,3,1), ncol = 2, byrow = T),
  id6358 = matrix(c(1,2,2,1,1,4,4,2,4,3,3,2,3,1), ncol = 2, byrow = T),
  id6364 = matrix(c(1,2,1,4,4,2,4,3,3,2,3,1,4,1), ncol = 2, byrow = T),
  id6366 = matrix(c(1,2,1,4,4,2,4,3,3,2,3,1,4,1,2,1), ncol = 2, byrow = T),
  id6550 = matrix(c(1,3,3,1,1,2,2,1,1,4,4,3,4,2), ncol = 2, byrow = T),
  id6552 = matrix(c(1,3,1,2,1,4,4,1,4,2,4,3), ncol = 2, byrow = T),
  id6554 = matrix(c(1,3,1,2,1,4,4,1,4,2,4,3,2,1), ncol = 2, byrow = T),
  id6558 = matrix(c(1,3,1,2,1,4,4,1,4,2,4,3,2,1,3,1), ncol = 2, byrow = T),
  id6598 = matrix(c(1,3,3,1,2,1,1,4,3,2,4,3,4,2), ncol = 2, byrow = T),
  id6602 = matrix(c(1,3,2,1,1,4,4,1,3,2,4,3,4,2), ncol = 2, byrow = T),
  id6604 = matrix(c(1,3,3,1,1,4,4,1,3,2,4,3,4,2), ncol = 2, byrow = T),
  id6606 = matrix(c(1,3,3,1,1,4,4,1,2,1,3,2,4,3,4,2), ncol = 2, byrow = T),
  id6614 = matrix(c(1,3,3,1,1,2,2,1,3,2,1,4,4,2,4,3), ncol = 2, byrow = T),
  id6616 = matrix(c(1,3,1,2,1,4,4,1,3,2,4,2,4,3), ncol = 2, byrow = T),
  id6618 = matrix(c(1,2,2,1,1,3,1,4,4,1,3,2,4,2,4,3), ncol = 2, byrow = T),
  id6620 = matrix(c(1,3,3,1,1,2,1,4,4,1,3,2,4,2,4,3), ncol = 2, byrow = T),
  id6622 = matrix(c(1,2,2,1,1,3,3,1,1,4,4,1,4,2,4,3), ncol = 2, byrow = T),
  id6854 = matrix(c(3,1,2,1,1,4,2,3,3,2,4,2,4,3), ncol = 2, byrow = T),
  id6858 = matrix(c(2,1,2,3,3,2,1,4,4,1,4,2,4,3), ncol = 2, byrow = T),
  id6862 = matrix(c(3,1,2,1,1,4,4,1,2,3,3,2,4,2,4,3), ncol = 2, byrow = T),
  id6870 = matrix(c(1,2,2,1,3,1,1,4,2,3,3,2,4,2,4,3), ncol = 2, byrow = T),
  id6874 = matrix(c(1,2,2,1,1,4,4,1,2,3,3,2,4,2,4,3), ncol = 2, byrow = T),
  id6876 = matrix(c(1,2,3,1,1,4,4,1,2,3,3,2,4,2,4,3), ncol = 2, byrow = T),
  id6878 = matrix(c(3,1,1,2,2,1,1,4,4,1,2,3,3,2,4,2,4,3), ncol = 2, byrow = T),
  id7126 = matrix(c(1,3,3,1,1,2,2,1,2,3,3,2,1,4,4,2,4,3), ncol = 2, byrow = T),
  id7128 = matrix(c(1,3,1,2,1,4,4,1,2,3,3,2,4,2,4,3), ncol = 2, byrow = T),
  id7130 = matrix(c(1,2,2,1,1,3,1,4,4,1,2,3,3,2,4,2,4,3), ncol = 2, byrow = T),
  id7134 = matrix(c(1,2,2,1,1,3,3,1,1,4,4,1,2,3,3,2,4,2,4,3), ncol = 2, byrow = T),
  id13142 = matrix(c(1,2,2,1,1,3,3,1,2,3,3,2,1,4,2,4), ncol = 2, byrow = T),
  id13146 = matrix(c(1,2,2,1,1,3,1,4,4,1,2,3,3,2,2,4), ncol = 2, byrow = T),
  id13148 = matrix(c(1,3,3,1,1,2,2,4,1,4,4,1,2,3,3,2), ncol = 2, byrow = T),
  id13150 = matrix(c(1,3,3,1,1,2,2,1,1,4,4,1,2,3,3,2,2,4), ncol = 2, byrow = T),
  id13260 = matrix(c(1,3,3,1,2,4,4,2,1,4,4,1,2,3,3,2), ncol = 2, byrow = T),
  id13262 = matrix(c(1,3,3,1,2,4,4,2,1,4,4,1,2,3,3,2,2,1), ncol = 2, byrow = T),
  id13278 = matrix(c(1,3,3,1,2,4,4,2,1,4,4,1,2,3,3,2,1,2,2,1), ncol = 2, byrow = T),
  id14678 = matrix(c(1,3,3,1,2,4,1,4,3,2,1,2,2,1,4,3), ncol = 2, byrow = T),
  id14686 = matrix(c(1,2,2,1,1,3,3,1,1,4,4,1,3,2,2,4,4,3), ncol = 2, byrow = T),
  id14790 = matrix(c(1,3,3,1,2,1,1,4,3,2,2,4,4,2,4,3), ncol = 2, byrow = T),
  id14798 = matrix(c(1,3,3,1,1,4,4,1,3,2,2,1,2,4,4,2,4,3), ncol = 2, byrow = T),
  id14810 = matrix(c(1,2,2,1,1,3,1,4,4,1,3,2,2,4,4,2,4,3), ncol = 2, byrow = T),
  id14812 = matrix(c(1,3,3,1,1,2,1,4,4,1,3,2,4,3,4,2,2,4), ncol = 2, byrow = T),
  id14814 = matrix(c(1,3,3,1,1,2,1,4,4,1,3,2,4,3,4,2,2,4,2,1), ncol = 2, byrow = T),
  id15258 = matrix(c(1,3,1,4,4,1,1,2,2,1,2,3,2,4,4,2,4,3), ncol = 2, byrow = T),
  id15262 = matrix(c(1,3,3,1,1,4,4,1,1,2,2,1,2,3,2,4,4,2,4,3), ncol = 2, byrow = T),
  id15310 = matrix(c(1,3,3,1,2,1,1,4,4,1,2,3,3,2,4,3), ncol = 2, byrow = T),
  id15326 = matrix(c(1,3,3,1,1,2,2,1,1,4,4,1,2,3,3,2,2,4,4,2,4,3), ncol = 2, byrow = T),
  id31710 = matrix(c(1,3,3,1,1,2,2,1,1,4,4,1,2,3,3,2,2,4,4,2,4,3,3,4), ncol = 2, byrow = T)
)


# Correct an error in manual entering of the graphs, switching the columns
fournode2 <- list()
for(i in 1:length(fournode)){
  m <- matrix(ncol = 2, nrow = nrow(fournode[[i]]))
  m[,1] <- fournode[[i]][,2]
  m[,2] <- fournode[[i]][,1]
  fournode2[[i]] <- m
}
rm(fournode)

# Convert adjacency lists into graph objects

fournode.gr <- lapply(fournode2, graph.edgelist)

# Create adjacency matrices
fournode.am <- lapply(fournode.gr, get.adjacency, sparse = F)

# count motifs
mot4 <- motif_counter(fournode.gr)
rownames(mot4) <- names(fournode2)


# run the stability analysis on four node subgraphs
fourN.co <- lapply(fournode.am, conversion)
system.time(
  emat <- eig.analysis(1000, fourN.co)
)
#2.33 min

nedges <- sapply(fournode.am, sum)
qss4 <- apply(emat, 2, function(x){sum(x < 0)/length(x)})
names(qss4) <- names(fournode.am)
plot(qss4[order(nedges)])


m4 <- matrix(nrow = nrow(mot4), ncol = 13)
for(i in 1:nrow(mot4)){
  m4[i,] <- unlist((mot4[i,]-mean(unlist(mot4[i,])))/sd(mot4[i,])) 
}
colnames(m4) <- colnames(mot4)
m4 <- as.data.frame(m4)
fit4.full <- glm(cbind(1000*qss4[nedges == 7], 1000-(1000*qss4[nedges == 7]))~m4[nedges == 7,], family = "binomial")
summary(fit4.full)
fitr <- glm(cbind(1000*qss4, 1000-(1000*qss4))~as.matrix(mot4), family = "binomial")
summary(fitr)


boxplot(qss4~nedges)

spmot4 <- split(mot4, nedges)
spqss4 <- split(qss4, nedges)

sapply(spmot4, colSums)

summary(lm(spqss4[["7"]]~as.matrix(spmot4[["7"]])[,1:5]))
# those with decent sample size
sub <- which(nedges >=4 & nedges <=9)

## Five Nodes

fivenode <- list()
system.time(
for(i in 1:1000){
  connect = FALSE
  while(!connect){
    e <- sample(4:20, 1, prob = c(rep(.1, 7), rep(.03, 10)))
    fivenode[[i]] <- erdos.renyi.game(5, p.or.m = e, type = "gnm", directed = T)  
    
    comp <- c()
    for(j in 1:length(fivenode)){
      if(length(fivenode) == 1){break}
      if(j > 1){j <- j-1}
      test <- sum(get.adjacency(fivenode[[i]]) == get.adjacency(fivenode[[j]])) == 25
      comp <- c(comp, test)
    }
    if(sum(comp) >= 1){connect = FALSE}else{connect <- is.connected(fivenode[[i]])}
  }
  if(i%%100 == 0) print(i)
}
)
#32 min

# Create adjacency matrices
fivenode.am <- lapply(fivenode, get.adjacency, sparse = F)


# count motifs
mot5 <- motif_counter(fivenode)


# five node stability
fiveN.co <- lapply(fivenode.am, conversion)
system.time(
  emat5 <- eig.analysis(1000, fiveN.co)
)
#9.1 min for 750 matrices

qss5 <- apply(emat5, 2, function(x){sum(x < 0)/length(x)})
hist(qss5)

apply(mot5[qss5 >= .8,], 2, mean)
apply(mot5[qss5 < .8 & qss5 >= .3,], 2, mean)
apply(mot5[qss5 < .3,], 2, mean)


apply(mot4[qss4 >= .8,], 2, mean)
apply(mot4[qss4 < .8 & qss4 >= .3,], 2, mean)
apply(mot4[qss4 < .3,], 2, mean)

# larger webs

tens <- list()
for(i in 1:1000){
  connect = FALSE
  while(!connect){
    con <- rbeta(1, 1, 3)
    tens[[i]] <- erdos.renyi.game(10, p.or.m = con, type = "gnp", directed = T)  
    
    comp <- c()
    for(j in 1:length(tens)){
      if(length(tens) == 1){break}
      if(j > 1){j <- j-1}
      test <- sum(get.adjacency(tens[[i]]) == get.adjacency(tens[[j]])) == 100
      comp <- c(comp, test)
    }
    if(sum(comp) >= 1){connect = FALSE}else{connect <- is.connected(tens[[i]])}
  }
  #print(i)
  if(i%%100 == 0) print(i)
}


mat10 <- lapply(tens, get.adjacency, sparse = F)
ned <- sapply(mat10, sum)
cmat10 <- lapply(mat10, conversion)


system.time(
  eig10 <- eig.analysis(1000, cmat10)
)
# 4.47 min for 200 matrices

qss10 <- apply(eig10, 2, function(x){sum(x < 0)/length(x)})
hist(qss10)

mot10 <- motif_counter(tens)

apply(mot10[qss10 >= .8,], 2, mean)
apply(mot10[qss10 < .8 & qss10 >= .3,], 2, mean)
apply(mot10[qss10 < .3,], 2, mean)




qs <- c(qss4, qss5, qss10)
rs <- c(r, r5, r10)

plot(qs~rs)
fits <- glm(cbind(1000*qs, 1000-(1000*qs))~rs, family = "binomial")
points(sort(fits$fitted.values, decreasing = F)~sort(rs, decreasing = F), col = "blue", typ = "o")



niche.model<-function(S,C){
  require(igraph)
  connected = FALSE
  while(!connected){  
    new.mat<-matrix(0,nrow=S,ncol=S)
    ci<-vector()
    niche<-runif(S,0,1)
    r<-rbeta(S,1,((1/(2*C))-1))*niche
    
    for(i in 1:S){
      ci[i]<-runif(1,r[i]/2,niche[i])
    }
    
    r[which(niche==min(niche))]<-.00000001
    
    for(i in 1:S){
      
      for(j in 1:S){
        if(niche[j]>(ci[i]-(.5*r[i])) && niche[j]<(ci[i]+.5*r[i])){
          new.mat[j,i]<-1
        }
      }
    }
    
    new.mat<-new.mat[,order(apply(new.mat,2,sum))]
    
    connected <- is.connected(graph.adjacency(new.mat))
  }
  return(new.mat)
}

niche_maker <- function(n, S, C){
  niche.list <- list()
  for (i in 1:n){
    niche.list[[i]]<- niche.model(S, C)
  }
  return(niche.list)
}


n1 <- niche_maker(1000, 15, .1)
N.co <- lapply(n1, conversion)
system.time(
  ematN <- eig.analysis(1000, N.co)
)


qssN <- apply(ematN, 2, function(x){sum(x < 0)/length(x)})
motN <- motif_counter(lapply(n1, graph.adjacency))


apply(motN[qssN >= .8,], 2, mean)
apply(motN[qssN < .8 & qss10 >= .3,], 2, mean)
apply(motN[qssN < .3,], 2, mean)



#############################################
fw.mot <- read.csv("https://raw.githubusercontent.com/jjborrelli/Ecological-Networks/master/FoodWebs/Tables/motifCOUNTS.csv", row.names = 2)[,-1]

fw.r <- rowSums(fw.mot[,c(1,4,5)]+1)/rowSums(fw.mot[,-c(1,4,5)]+1)

fw.pred <- inv.logit(fitN$coefficients[1] + fitN$coefficients[2]*fw.r )
points(fw.pred~fw.r, col = "green", pch = 17)

load("C:/Users/jjborrelli/Desktop/GitHub/Subgraph-Stability/webGRAPHS.Rdata")
web.adj <- lapply(web.graphs, get.adjacency, sparse = F)
fw.eig <- eig.analysis(10000, web.adj)
fw.qss <- apply(fw.eig, 2, function(x){sum(x < 0)})
points(fw.qss~fw.r, col = "darkgreen", pch = 17)

cor.test(qssN, rN)
nc <- c()
for(i in 1:1000){
  nqss <- sample(qssN)
  nc[i] <- cor.test(nqss, rN)$estimate
}
hist(nc)

Ndf <- data.frame(s = 10000*qssN, f = 10000-(10000*qssN), rat = rN)

each.rN<- rep(rN, each = 10000)
bin = c()
for(i in 1:nrow(Ndf)){
  bin <- c(bin, rep(1, Ndf[i,1]), rep(0, Ndf[i,2]))
}

ndf <- data.frame(bin, each.rN)
ggplot(ndf, aes(x = each.rN, y = bin)) + geom_point(data = Ndf, aes(x = rN, y = s/10000), alpha = .5) + stat_smooth(method = "glm", family = "binomial", size = 1.5) + theme_bw() + xlab("ratio") + ylab("Quasi sign-stability")

#qssN[which(qssN == 0)] <- qssN[which(qssN == 0)] + 10^-7
#qssN[which(qssN == 1)] <- qssN[which(qssN == 1)] - 10^-7

breg <- betareg(qssN~rN)
summary(breg)
summary(fitN)

pcM <- princomp(motN)
pcM$loadings
biplot(pcM)

fitN <- glm(cbind(10000*qssN, 10000-(10000*qssN))~con+as.matrix(motN), family = "binomial")
summary(fitN)

m2 <- matrix(nrow = nrow(motN), ncol = 13)
for(i in 1:nrow(motN)){
  m2[i,] <- unlist((motN[i,] - mean(unlist(motN[i,])))/sd(motN[i,]))
  #m2[i,] <- unlist(motN[i,]/nume[i])
  #m2[i,] <- unlist(motN[i,]/sum(sqrt(motN[i,])))
}
colnames(m2) <- colnames(motN)
fitN <- glm(cbind(10000*qssN, 10000-(10000*qssN))~m2, family = "binomial")
summary(fitN)

par(mfrow = c(3,5))
for(i in 1:13){
  plot(qssN~m2[,i], xlab = colnames(motN)[i])
}

fitN <- glm(cbind(10000*qssN, 10000-(10000*qssN))~motN$s1+motN$s2+motN$s3+motN$s4+motN$s5+motN$d1+motN$d2+motN$d3+motN$d4+motN$d5+motN$d6+motN$d7+motN$d8, family = "binomial")
m2 <- as.data.frame(m2)
fitN <- glm(cbind(10000*qssN, 10000-(10000*qssN))~m2$s1+m2$s2+m2$s3+m2$s4+m2$s5+m2$d1+m2$d2+m2$d3+m2$d4+m2$d5+m2$d6+m2$d7+m2$d8, family = "binomial")
summary(fitN)
hist(m2$d8)


testN <- t(apply(as.matrix(motN), 1, function(x){x/sum(x)})) 
fit.test <- glm(cbind(10000*qssN, 10000-(10000*qssN))~testN, family = "binomial")



##############################################################################
curve_ball<-function(m){
  RC=dim(m)
  R=RC[1]
  C=RC[2]
  hp=list()
  for (row in 1:dim(m)[1]) {hp[[row]]=(which(m[row,]==1))}
  l_hp=length(hp)
  for (rep in 1:5*l_hp){
    AB=sample(1:l_hp,2)
    a=hp[[AB[1]]]
    b=hp[[AB[2]]]
    ab=intersect(a,b)
    l_ab=length(ab)
    l_a=length(a)
    l_b=length(b)
    if ((l_ab %in% c(l_a,l_b))==F){
      tot=setdiff(c(a,b),ab)
      l_tot=length(tot)
      tot=sample(tot, l_tot, replace = FALSE, prob = NULL)
      L=l_a-l_ab
      hp[[AB[1]]] = c(ab,tot[1:L])
      hp[[AB[2]]] = c(ab,tot[(L+1):l_tot])}
    
  }
  rm=matrix(0,R,C)
  for (row in 1:R){rm[row,hp[[row]]]=1}
  rm
}


the.nm <- niche.model(10, .2)
mynm <- the.nm
mymc <- motif_counter(graph.lists = list(graph.adjacency(mynm)))
myqss <- sum(eig.analysis(1000, list(conversion(mynm))) < 0)/1000

iter = 1000
for(j in 1:iter){
  mynm <- curve_ball(mynm)
  new1 <- conversion(mynm)
  mymc <- rbind(mymc, motif_counter(list(graph.adjacency(mynm))))
  meig <- c()
  rho <- c()
  for(i in 1:1000){
    ru <- ran.unif(new1)
    meig[i] <- maxRE(ru)
    rho[i] <-cor.test() 
  }
  myqss <- c(myqss, sum(meig < 0)/1000)
  print(j)
}

fit1 <- glm(cbind(1000*myqss, 1000-(1000*myqss))~as.matrix(mymc), family = "binomial")
summary(fit1)

fit2 <- glm(cbind(1000*myqss, 1000-(1000*myqss))~mymc$s4+mymc$d4+mymc$s1+mymc$s5, family = "binomial")
summary(fit2)

fit.all <- glm(cbind(1000*myqss, 1000-(1000*myqss))~mymc$s1+mymc$s2+mymc$s3+mymc$s4+mymc$s5+mymc$d1+mymc$d2+mymc$d3+mymc$d4+mymc$d5+mymc$d6+mymc$d7+mymc$d8, family = "binomial")
summary(fit.all)
stepAIC(fit.all, direction = "both")


testN <- apply(motN, c(1,2), function(x){if(!x == 0){log(x)}else{0}}) 
fit.test <- glm(cbind(1000*qssN, 1000-(1000*qssN))~as.matrix(testN), family = "binomial")
