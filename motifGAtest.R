library(igraph)
library(magrittr)
library(pheatmap)
library(compiler)
library(ggplot2)
library(reshape2)

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


library(compiler)
ran.unif.cmp <- cmpfun(ran.unif)

maxRE <- function(rmat){
  lam.max <- max(Re(eigen(rmat)$values))
  return(lam.max)
}

eig.analysis <- function(n, matrices){
  cols <- length(matrices)
  rows <- n
  eigenMATRIX <- matrix(0, nrow = rows, ncol = cols)
  for(i in 1:n){
    ranmat <- lapply(matrices, ran.unif.cmp)
    
    eigs <- sapply(ranmat, maxRE)
    eigenMATRIX[i,] <- eigs
  }
  return(eigenMATRIX)
}

qss <- function(x){sum(x < 0)/ length(x)}

spec <- function(erg.fill, n){
  cond <- TRUE
  while(cond){
    sp <- sample(1:10, n)
    c1 <- erg.fill[sp[1], sp[2]] != 0 & erg.fill[sp[2], sp[1]] != 0
    c2 <- erg.fill[sp[3], sp[4]] == 0 & erg.fill[sp[4], sp[3]] == 0
    cond <- !c1 | !c2
  }
  return(sp)
}

swap_links <- function(mat, s){
  mat[s[3], s[4]] <- mat[s[1], s[2]] 
  mat[s[4], s[3]] <- mat[s[2], s[1]]
  
  mat[s[1], s[2]] <- 0
  mat[s[2], s[1]] <- 0
  
  return(mat)
}

spec.cmp <- cmpfun(spec)
swap_links.cmp <- cmpfun(swap_links)

rev.conv <- function(erg.fill){
  return(apply(erg.fill, c(1,2), function(x){if(x > 0){x <- 1}else if(x < 0){x <- 0}else{x <- 0}}))
}

eigenTRACE <- function(x, erg.fill){
  time <- 0
  mo <- motif_counter(list(graph.adjacency(rev.conv(erg.fill))))
  links <- sum(rev.conv(erg.fill))
  erg.ei <- maxRE(erg.fill)
  for(i in 1:x){
    # Choose links to swap
    sp <- spec.cmp(erg.fill, 4)
    # Swap links
    erg.perm <- swap_links.cmp(erg.fill, sp)
    # Get eigenvalue
    ei <- maxRE(erg.perm)
    
    # Check condition, and keep matrix that is better
    if(ei < tail(erg.ei, 1)){
      erg.fill <- erg.perm
      # keep eigenvalue
      erg.ei <- c(erg.ei, ei)
      # get binary adj matrix
      erg.m <- rev.conv(erg.fill)
      # get motif structure of new matrix
      mo <- rbind(mo, motif_counter(list(graph.adjacency(erg.m))))
      # what iteration are we at
      time <- c(time, i) 
      # plot eigenvalue change through time
      #plot(time, erg.ei)
      links <- c(links, sum(erg.m)) 
      if(sum(erg.m) < links[1]){return(list(erg.m, erg.fill));stop("you lost a link")}
    }
  }
  return(list(time = time, eigen = erg.ei, motif = as.matrix(mo), links = links))
}

init_cond <- function(x){
  require(igraph)
  connect <- FALSE
  while(!connect){
    erg.g <- erdos.renyi.game(10, .2, "gnp", directed = TRUE)
    connect <- is.connected(erg.g)
  }
  
  ## Get adjacency matrix of random network
  erg.m <- get.adjacency(erg.g, sparse = F)
  ## Get sign matrix
  erg.c <- conversion(erg.m)
  
  # Fill sign matrix x times to make a list
  filled <- list()
  for(i in 1:x){
    filled[[i]] <- ran.unif.cmp(erg.c)
  }
  return(filled)
}

####
Get the lowest leading eigenvalue
####
erg <- init_cond(1000)

system.time(
sim.results <- lapply(erg, eigenTRACE, x = 1000)
)
# 458 for 1000, 1000, 10sp .2C

times <- sapply(sim.results, "[[", 1) 
eigens <- sapply(sim.results, "[[", 2)
motifs <- sapply(sim.results, "[[",3)
links <- sapply(sim.results, "[[", 4)

plot(times[[1]], eigens[[1]], xlim = c(0, 1000), ylim = c(0, 15))
for(i in 2:1000){
  points(times[[i]], eigens[[i]])
}

motifs
boxplot(do.call(rbind, lapply(motifs, tail, 1)))

ggplot(melt(do.call(rbind, lapply(motifs, tail, 1))), aes(x = Var2, y = value)) + geom_point(position ="jitter", alpha = .15) + stat_summary(fun.y = "mean", fun.ymin = function(x){quantile(x, prob = .025)}, fun.ymax = function(x){quantile(x, prob = .975)}, geom = "errorbar", col = "darkgreen", lwd = 1.5) + stat_summary(fun.y = "mean", geom = "point", col = "darkgreen", size = 5) + theme_bw()

#meanMOT <- matrix(nrow = max(sapply(times, length)), ncol = 13)
#meanEIG <- c()
#for(i in 1:max(sapply(times, length))){
#  timestep.mot <- do.call(rbind, lapply(motifs, function(x){if(i <= nrow(x)){x[i,]}else{NA}}))
#  timestep.eig <- do.call(rbind, lapply(eigens, function(x){if(i <= length(x)){x[i]}else{NA}}))
#  meanMOT[i,] <- colMeans(timestep.mot, na.rm = T)
#  meanEIG[i] <- mean(timestep.eig, na.rm = T)
#}

#pheatmap(meanMOT, cluster_rows = F, cluster_cols = F)

last.struct <- do.call(rbind, lapply(motifs, tail, 1))
lower <- apply(last.struct, 2, quantile, prob = .025)
lower2 <- apply(do.call(rbind, lapply(motifs, tail, 1)), 2, quantile, prob = .025)
upper <- apply(last.struct, 2, quantile, prob = .975)
upper2 <- apply(do.call(rbind, lapply(motifs, tail, 1)), 2, quantile, prob = .975)

final.struct <- melt(do.call(rbind, lapply(motifs, tail, 1)))

sub.names <- c("s1", "s2", "s3", "s4", "s5", "d1", "d2", "d3", "d4", "d5", "d6", "d7", "d8")
df <- data.frame(Subgraph = factor(names(lower), levels = sub.names), 
                 means = colMeans(last.struct), 
                 means2 = colMeans(do.call(rbind, lapply(motifs, tail, 1))), 
                 upper, upper2, 
                 lower, lower2)

df2 <- data.frame(trial = rep(1:2, each = 13), Subgraph = rep(df$Subgraph, 2) ,
                  mean = c(df$means, df$means2), upper = c(df$upper, df$upper2), 
                  lower = c(df$lower, df$lower2))
ggplot(df2, aes(x = Subgraph, y = mean, fill = factor(trial))) + geom_point(stat = "identity", position = "dodge") + geom_errorbar(aes(ymax = upper, ymin = lower), position = "dodge") + theme_bw()

