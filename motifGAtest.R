library(igraph)
library(magrittr)
library(pheatmap)
library(compiler)
library(ggplot2)
library(reshape2)
library(data.table)

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
    if(x==1){rbeta(1, 2, 5)*10}else if(x==-1){runif(1, -1, 0)} else{0}
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

swap_values <- function(mat){
  new.mat <- mat
  d <- diag(mat)
  diag(new.mat) <- 0
  
  pos <- which(new.mat > 0)
  neg <- which(new.mat < 0)
  
  new.mat[pos] <- mat[sample(pos)]
  new.mat[neg] <- mat[sample(neg)]
  
  diag(new.mat) <- d
  
  return(new.mat)
}

spec.cmp <- cmpfun(spec)
swap_links.cmp <- cmpfun(swap_links)
swap_values.cmp <- cmpfun(swap_values)

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

eigenTRACE.v2 <- function(x, erg.fill){
  time <- 0
  erg.ei <- maxRE(erg.fill)
  erg.list <- list(erg.fill)
  for(i in 1:x){
    # Swap values
    erg.perm <- swap_values.cmp(erg.fill)
    # Get eigenvalue
    ei <- maxRE(erg.perm)
    
    # Check condition, and keep matrix that is better
    if(ei < tail(erg.ei, 1)){
      erg.fill <- erg.perm
      # keep eigenvalue
      erg.ei <- c(erg.ei, ei)
      
      # what iteration are we at
      time <- c(time, i) 
      
      # save the matrix
      erg.list[[length(erg.list)+1]] <- erg.fill
    }
  }
  return(list(ergs = erg.list, eigens = erg.ei, time = time))
}

init_cond <- function(x, S = 10, C = 0.2){
  require(igraph)
  connect <- FALSE
  while(!connect){
    erg.g <- erdos.renyi.game(S, C, "gnp", directed = TRUE)
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

link_structure <- function(motMATS){
  mots <- length(motMATS)
  # mots is how many networks there are
  
  class1 <- list()
  class2 <- list()
  class3 <- list()
  for(n in 1:mots){
    mat <- motMATS[[n]]
    diag(mat) <- 0
    mo1 <- unlist(motif_counter(list(graph.adjacency(mat))))
    t1 <- which(motMATS[[n]] > 0)
    tst.lst <- list()
    for(i in 1:length(t1)){
      t2 <- motMATS[[n]]
      t2[t1[i]] <- 0
      tst.lst[[i]] <- t2
    }
    gr <- lapply(tst.lst, graph.adjacency)
    mo <- as.matrix(motif_counter(gr))
    diff1 <- t(apply(mo, 1, function(x){x - mo1}))
    participation <- apply(diff1, c(1,2), function(x){if(x > 0){x <- 0}else if(x < 0){x <- x*-1}else{x <- 0}})
    
    with.unst <- apply(participation, 1, function(x){sum(x[3], x[6:13]) > 0})
    st <- apply(participation, 1, function(x){sum(x[1], x[4], x[5]) > 0})
    intr <- apply(participation, 1, function(x){x[2] > 0})
    
    class1[[n]] <- which(st & !with.unst & !intr)
    class2[[n]] <- which(st & !with.unst & intr)
    class3[[n]] <- which(with.unst)
  }
  return(list(very = class1, mod = class2, unst = class3))
}


classify <- function(fm, l1){
  m <- matrix(0, nrow(fm), ncol(fm))
  m[which(fm > 0)][l1$very] <- "VERY"
  m[which(fm > 0)][l1$mod] <- "MOD"
  m[which(fm > 0)][l1$unst] <- "UNST"
  
  df <- data.frame(a.ij = factor(m[upper.tri(m)], levels = c("0", "MOD", "UNST", "VERY")),
                   a.ji = factor(t(m)[upper.tri(m)], levels = c("0", "MOD", "UNST", "VERY")))
  
  df2 <- factor(levels = c("0", "MOD", "UNST", "VERY"))
  for(i in 1:nrow(df)){
    if(df$a.ij[i] == 0 & df$a.ji[i] != 0){
      df2[i] <- factor(df$a.ji[i], levels = c("0", "MOD", "UNST", "VERY"))
    }else if(df$a.ij[i] != 0 & df$a.ji[i] == 0){
      df2[i] <- factor(df$a.ij[i], levels = c("0", "MOD", "UNST", "VERY"))
    }else if(df$a.ij[i] != 0 & df$a.ji[i] != 0 & df$a.ij[i] == df$a.ji[i]){
      df2[i] <- factor(df$a.ij[i], levels = c("0", "MOD", "UNST", "VERY"))
    }else{df2[i] <- factor(0, levels = c("0", "MOD", "UNST", "VERY"))}
  }
  return(df2)
}

####
# Get the lowest leading eigenvalue
####
n.web = 10
mag.class <- list()
for(i in 1:10){  
  erg <- init_cond(n.web, S = 150, C = .1)
  
  #system.time(
  sim.results <- lapply(erg, eigenTRACE.v2, x = 2000)
  #)
  # 458 for 1000, 1000, 10sp .2C
  
  # Get all matrices
  allmats <- sapply(sim.results, "[[", 1) 
  # Get eigenvalues for each matrix
  eigens <- sapply(sim.results, "[[", 2)
  # Get time stamp
  times <- sapply(sim.results, "[[",3)
  
  # Get final matrix structure
  fin.mats <- unlist(lapply(allmats, tail, 1), recursive = F) 
  lstr <- link_structure(fin.mats)
  
  pairs <- lapply(allmats, function(x){lapply(x, function(z){data.frame(aij = z[upper.tri(z)], aji = t(z)[upper.tri(z)])})})
  allpair <- rbindlist(unlist(lapply(pairs, function(x){tail(x, 1)}), recursive = F))
  
  d.all <- list()
  for(x in 1:n.web){
    class.l <- list(very = lstr$very[[x]], mod = lstr$mod[[x]], unst = lstr$unst[[x]])
    d.all[[x]] <- classify(fin.mats[[x]], class.l)
    #print(x)
  }
  
  mag.class[[i]] <- cbind(allpair, class = factor(unlist(d.all), levels = c("0", "UNST", "MOD", "VERY")))
  
  print(i)
}


mag.class.r <- rbindlist(mag.class)

mmag <- melt(mag.class.r)
mmag$sign <- NA
mmag$sign[which(mmag$value < 0)] <- "neg"
mmag$sign[which(mmag$value > 0)] <- "pos"
mmag$sign[which(mmag$value == 0)] <- "zero"

mmag <- mmag[which(mmag$sign != "zero" & mmag$sign == "pos" & mmag$class != 0)]

ggplot(mmag, aes(x = value, y = ..density..)) + geom_density(alpha = .5, aes(fill = class)) + facet_grid(class~sign)
ggsave("C:/Users/jjborrelli/Dropbox/distros3.jpg", width = 5, height = 5)

aggregate(mmag$value, list(mmag$class, mmag$sign), function(x){c(mean(x), median(x), sd(x))})
