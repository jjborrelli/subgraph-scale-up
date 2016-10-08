library(igraph)
library(igraphdata)
library(deSolve)
library(rend)
library(NetIndices)
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



data("foodwebs")

mich <- get.adjacency(foodwebs$Michigan, sparse = F)
dim(mich)

mm <- combn(39, 3)
dim(mm)

timing <- c()
webs <- c(1, 4, 7, 9, 10, 11, 12, 13, 15, 17, 19)
for(i in 1:length(webs)){
  strt <- Sys.time()
  mm <- combn(length(V(foodwebs[[webs[i]]])), 3)
  mich <- get.adjacency(foodwebs[[webs[i]]], sparse = F)
  tbt <- lapply(1:ncol(mm), function(x){mich[mm[,x], mm[,x]]})
  conn <- sapply(tbt, function(x) is.connected(graph.adjacency(x)))
  mots <- t(sapply(tbt[conn], function(x) motifs(graph.adjacency(x))))
  end <- Sys.time()
  
  timing[i] <- end - strt
  
  print(i)
}
timing
plot(timing~sapply(foodwebs[webs], function(z) length(V(z))))
foodwebs[webs]


getmot <- function(web){
  com <- combn(length(V(web)), 3)
  adj <- get.adjacency(web, sparse = F)
  
  tbt <- lapply(1:ncol(com), function(x){adj[com[,x], com[,x]]})
  conn <- sapply(tbt, function(x) is.connected(graph.adjacency(x)))
  mots <- motif_counter(lapply(tbt[conn], graph.adjacency))
  
  m1 <- melt(lapply(1:sum(conn), function(x) com[,x]))
  m2 <- mots[m1$L1,]
  m3 <- aggregate(m2, list(m1$value), sum)
  
  
  return(m3)
}

ext1 <- function(times, states, parms){
  with(as.list(states),{
    states[states < (10^(-10))] <- 0
    return(c(states))
  })
}

n.vals <- function(S,C){
  niche<-runif(S,0,1)
  r<-rbeta(S,1,((1/(2*C))-1))*niche
  
  ci <- runif(S, r/2, niche)
  
  return(cbind(niche = niche, ci = ci, r = r))
}

n.mat <- function(nv){
  S <- nrow(nv)
  new.mat<-matrix(0,nrow=S,ncol=S)
  
  for(i in 1:S){
    
    for(j in 1:S){
      if(nv[j,1]>(nv[i,2]-(.5*nv[i,3])) && nv[j,1]<(nv[i,2]+.5*nv[i,3])){
        new.mat[j,i]<-1
      }
    }
  }
  
  #new.mat<-new.mat[order(apply(new.mat,2,sum)),order(apply(new.mat,2,sum))]
  return(new.mat)
}


system.time(tmot <- getmot(foodwebs$ChesLower))
test <- rowSums(tmot)
sum(test == rowSums(tmot, na.rm = T))

system.time(tmot <- getmot(graph.adjacency(nm1)))


niches <- list()
adjs <- list()
webs <- list()
for(i in 1:250){
  cond <- FALSE
  while(!cond){
    niches[[i]] <- n.vals(60, .1)
    adjs[[i]] <- n.mat(niches[[i]])
    webs[[i]] <- graph.adjacency(adjs[[i]])
    
    cond <- is.connected(webs[[i]])
  }
}


#niches <- lapply(1:250, function(x) n.vals(50, .15))
#adjs <- lapply(niches, n.mat)
#webs <- lapply(adjs, graph.adjacency)

system.time(
roles <- lapply(webs, getmot)
)

tind <- lapply(adjs, TrophInd)
states1 <- t(sapply(1:250, function(x) runif(60, .5, 1)))

### xpar = 0.2
x2 <- Sys.time()
allpar <- lapply(1:250, function(x){
  xi <- (((10^2)^tind[[x]]$TL)/100)^-0.25*0.314
  ri <- as.numeric(colSums(adjs[[x]]) == 0)
  
  par <- list(K = 1, x.i = xi, yij = 8, eij = 0.85, xpar = 0.2, 
              B.o = 0.5, r.i = ri, A = adjs[[x]], G.i = Gi, FR = Fij)
})


dyna.2 <- list()
for(i in 1:250){
  dyna.2[[i]] <- ode(y = states1[i,], times = 1:2000, func = CRmod, parms = allpar[[i]],
              events = list(func = goExtinct, time = 1:2000))
  cat(i, "--- ")
}
x2end <- Sys.time()
x2end-x2

### xpar = 0

x0 <- Sys.time()

allpar0 <- lapply(1:250, function(x){
  xi <- (((10^2)^tind[[x]]$TL)/100)^-0.25*0.314
  ri <- colSums(adjs[[x]]) == 0
  
  par <- list(K = 1, x.i = xi, yij = 8, eij = 0.85, xpar = 0, 
              B.o = 0.5, r.i = ri, A = adjs[[x]], G.i = Gi, FR = Fij)
})

dyna0 <- list()
for(i in 1:250){
  dyna0[[i]] <- ode(y = states1[i,], times = 1:2000, func = CRmod, parms = allpar0[[i]],
                   events = list(func = ext1, time = 1:2000))
  cat(i, "--- ")
}
x0end <- Sys.time()
x0end-x0

### xpar = 1

x1 <- Sys.time()

allpar1 <- lapply(1:250, function(x){
  xi <- (((10^2)^tind[[x]]$TL)/100)^-0.25*0.314
  ri <- colSums(adjs[[x]]) == 0
  
  par <- list(K = 1, x.i = xi, yij = 8, eij = 0.85, xpar = 1, 
              B.o = 0.5, r.i = ri, A = adjs[[x]], G.i = Gi, FR = Fij)
})

dyna1 <- lapply(1:250, function(x) matrix(NA, nrow = 2000, ncol = 51))
for(i in 1:250){
  dyna1[[i]] <- ode(y = states1[i,], times = 1:2000, func = CRmod, parms = allpar1[[i]],
                   events = list(func = ext1, time = 1:2000), atol = 0)
  cat(i, "--- ")
}

x1end <- Sys.time()
x1end-x1


sapply(dyna, function(x) sum(is.nan(x)) > 1)
sapply(dyna0, function(x) nrow(x))

dyna.2[[1]]
