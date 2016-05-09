#################################
## Modified Fast Marching Code ##
#################################
### PARAMETERS:
## seeds    :: matrix of NA values everywhere, 
##             except at seed locations where the values should be equal to inception time t0
## n        :: baseline speed
## res      :: matrix resolution, in spatial units
## boosts   :: matrix of boost values, with 0's for unreacheable cells
#################################

ModifiedFastMarch <- function(seeds, n, res, boosts)
{
#### Parameter initialization
  time.step <- 1000;
  v <- time.step/res;
  speeds <- n*v;
  
  ## Create key matrices from boosts layer
  Map <- boosts;
  gridsize <- dim(Map)
  
  ## Set up seeds
  aux <- which(!is.na(seeds), arr.ind=T);
  diff.t <- array(0,c(1,length(aux[,1])));
  for(i in 1:length(aux[,1])) { diff.t[i] <- seeds[aux[i,1],aux[i,2]] }
  delta.t <- diff.t/time.step;
  coords <- rbind(t(aux[,1]),t(aux[,2]));
  seeds <- rbind(coords,delta.t,speeds);
  
  
#### Fast Marching initialization
  # sort seeds by inception time
  sources <- seeds[,order(seeds[3,])]; dim(sources) <- dim(seeds);
  V <- sources[4,];
  
  # retrieve seeds that birth at t = 0
  t0.sources <- which(sources[3,]==0);
  source.points <- sources[1:2,t0.sources];  dim(source.points) <- c(2,dim(seeds)[2]);
  
  # other initializations
  F <- array(1,c(gridsize[1],gridsize[2]));
  proc <- array(0,c(gridsize[1],gridsize[2]));  # process name/label
  T <- array(0,dim(F));
  Tm2 <- array(0,c(1,4));
  Frozen <- Map==0;
  
  # allocate memory to store neighbours of the (segmented) region
  neg_free <- 100000; neg_pos <- 0; neg_list <- array(0,c(4,neg_free)); 
  ne <- rbind(c(-1,0),c(1,0),c(0,-1),c(0,1));
  source.points <- floor(source.points);
  
  
#### Fast Marching Run
  # insert seeds
  x <- source.points[1,]; y <- source.points[2,];
  for (kk in 1:length(t0.sources)) {
    Frozen[x[kk],y[kk]] <- 1; T[x[kk],y[kk]] <- 0;
    proc[x[kk],y[kk]] <- length(t0.sources)-dim(source.points)[2]+kk;
    
    # add neighbours of seeds to narrow list
    for (k in 1:4) {
      i <- x[kk] + ne[k,1]; j <- y[kk] + ne[k,2];
      
      if ((i>0)&&(j>0)&&(i<=dim(F)[1])&&(j<=dim(F)[2])&&(Frozen[i,j]==0)) {
        proc[i,j] <- proc[x[kk],y[kk]];
  
        Tt <- (1/(V[proc[i,j]]+.Machine$double.eps));
        
        if (T[i,j]>0) {
          neg_list[1,T[i,j]] <- min(Re(Tt),neg_list[1,T[i,j]]);
          neg_list[4,T[i,j]] <- proc[x[kk],y[kk]];
        } else {
          neg_pos <- neg_pos+1;
          if (neg_pos>neg_free) {
            neg_free <- neg_free + 100000;
            neg_list[1,neg_free] <- 0;
          }
          neg_list[,neg_pos] <- rbind(Re(Tt),i,j,proc[x[kk],y[kk]]);
          T[i,j] <- neg_pos;
        }
      }
    }
    
  }
  
  # loop through all points
  for (itt in 1:length(F)) {
    if(neg_pos==0) { break }
    
    # get the closest pixel to wavefront and make it the current pixel
    t <- min(neg_list[1,1:neg_pos]); index <- which.min(neg_list[1,1:neg_pos]);
    x <- neg_list[2,index]; y <- neg_list[3,index];
    Frozen[x,y] <- 1; T[x,y] <- neg_list[1,index];
    proc[x,y] <- neg_list[4,index];
    
    if (length(index)<1) {browser()}
    
    # replace min value with the last value in the array
    if(index < neg_pos) {
      neg_list[,index] <- neg_list[,neg_pos];
      x2 <- neg_list[2,index]; y2 <-  neg_list[3,index];
      T[x2,y2] <- index;
    }
    neg_pos <- neg_pos-1;
    
    # loop through all neighbours of current pixel
    for (k in 1:4) {
      i <- x + ne[k,1]; j <- y + ne[k,2];
      
      # check if neighbour is not yet frozen
      if((i>0)&&(j>0)&&(i<=dim(F)[1])&&(j<=dim(F)[2])&&(Frozen[i,j]==0)&&(proc[i,j]==0)) {
        Tpatch <- matrix(Inf,5,5)
        for (nx in -2:2) {
          for (ny in -2:2) {
            i.n <- i + nx; j.n <- j + ny;
            if((i.n>0)&&(j.n>0)&&(i.n<=dim(F)[1])&&(j.n<=dim(F)[2])&&(Frozen[i.n,j.n]==1)&&(proc[i.n,j.n]==proc[x,y])) { Tpatch[nx+3,ny+3] <- T[i.n,j.n]}
          }
        }
        
        # this will store the order of derivative to use
        Order <- array(0,c(1,4));
        
        # 1st order derivatives in x-y and cross directions
        Tm <- array(0,c(1,4));
        Tm[1] <- min(c(Tpatch[2,3],Tpatch[4,3])); if(is.finite(Tm[1])) {Order[1] <- 1}
        Tm[2] <- min(c(Tpatch[3,2],Tpatch[3,4])); if(is.finite(Tm[2])) {Order[2] <- 1}
        Tm[3] <- min(c(Tpatch[2,2],Tpatch[4,4])); if(is.finite(Tm[3])) {Order[3] <- 1}
        Tm[4] <- min(c(Tpatch[2,4],Tpatch[4,2])); if(is.finite(Tm[4])) {Order[4] <- 1}
        
        # make 2nd derivatives
        # pixels with a pixel distance 2 from the center must be lower in value otherwise use other side or first order
        Tm2 <- array(0,c(1,4));
        ch1 <- (Tpatch[1,3]<Tpatch[2,3])&&is.finite(Tpatch[2,3])
        ch2 <- (Tpatch[5,3]<Tpatch[4,3])&&is.finite(Tpatch[4,3])
        if(ch1) { Tm2[1] <- (4*Tpatch[2,3]-Tpatch[1,3])/3; Order[1] <- 2 }
        if(ch2) { Tm2[1] <- (4*Tpatch[4,3]-Tpatch[5,3])/3; Order[1] <- 2 }
        if(ch1&&ch2) { Tm2[1] <- min(c((4*Tpatch[2,3]-Tpatch[1,3])/3,(4*Tpatch[4,3]-Tpatch[5,3])/3)); Order[1] <- 2;}
        
        ch1 <- (Tpatch[3,1]<Tpatch[3,2])&&is.finite(Tpatch[3,2])
        ch2 <- (Tpatch[3,5]<Tpatch[3,4])&&is.finite(Tpatch[3,4])
        if(ch1) { Tm2[2] <- (4*Tpatch[3,2]-Tpatch[3,1])/3; Order[2] <- 2 }
        if(ch2) { Tm2[2] <- (4*Tpatch[3,4]-Tpatch[3,5])/3; Order[2] <- 2 }
        if(ch1&&ch2) { Tm2[2] <- min(c((4*Tpatch[3,2]-Tpatch[3,1])/3,(4*Tpatch[3,4]-Tpatch[3,5])/3)); Order[2] <- 2;}
        
        ch1 <- (Tpatch[1,1]<Tpatch[2,2])&&is.finite(Tpatch[2,2])
        ch2 <- (Tpatch[5,5]<Tpatch[4,4])&&is.finite(Tpatch[4,4])
        if(ch1) { Tm2[3] <- (4*Tpatch[2,2]-Tpatch[1,1])/3; Order[3] <- 2 }
        if(ch2) { Tm2[3] <- (4*Tpatch[4,4]-Tpatch[5,5])/3; Order[3] <- 2 }
        if(ch1&&ch2) { Tm2[3] <- min(c((4*Tpatch[2,2]-Tpatch[1,1])/3,(4*Tpatch[4,4]-Tpatch[5,5])/3)); Order[3] <- 2;}
        
        ch1 <- (Tpatch[1,5]<Tpatch[2,4])&&is.finite(Tpatch[2,4])
        ch2 <- (Tpatch[5,1]<Tpatch[4,2])&&is.finite(Tpatch[4,2])
        if(ch1) { Tm2[4] <- (4*Tpatch[2,4]-Tpatch[1,5])/3; Order[4] <- 2 }
        if(ch2) { Tm2[4] <- (4*Tpatch[4,2]-Tpatch[5,1])/3; Order[4] <- 2 }
        if(ch1&&ch2) { Tm2[4] <- min(c((4*Tpatch[2,4]-Tpatch[1,5])/3,(4*Tpatch[4,2]-Tpatch[5,1])/3)); Order[4] <- 2;}
        
        # calculates the distance using x-y and cross directions only
        Coeff <- c(0,0,-1/((V[proc[x,y]]*Map[x,y])^2));
        for (t in 1:2) { if(Order[t]>0) { Coeff <- switch(Order[t],Coeff+c(1,-2*Tm[t],Tm[t]^2),Coeff+c(1,-2*Tm2[t],Tm2[t]^2)*9/4); }}
        Tt <- polyroot(rev(Coeff)); Tt <- max(Re(Tt));
        Coeff <- c(0,0,-1/((V[proc[x,y]]*Map[x,y])^2));
        for (t in 3:4) { if(Order[t]>0) { Coeff <- switch(Order[t],Coeff+0.5*c(1,-2*Tm[t],Tm[t]^2),Coeff+0.5*c(1,-2*Tm2[t],Tm2[t]^2)*9/4); }} 
        Tt2 <- polyroot(rev(Coeff));
        
        # selects minimum distance and check for upwind condition
        if(length(Tt2)>0) {Tt2 <- max(Re(Tt2)); Tt <- min(c(Tt,Tt2));}
        DirectNeigbInSol <- Tm[is.finite(Tm)];
        if(length(which(DirectNeigbInSol>=Tt))>0) { Tt <- min(DirectNeigbInSol) + (1/((V[proc[x,y]]*Map[x,y]))) }
        
        # updates neighbour list
        if(T[i,j]>0) {
          mT <- min(c(Tt,neg_list[1,T[i,j]])); ind <- which.min(c(Tt,neg_list[1,T[i,j]]));
          neg_list[1,T[i,j]] <- Re(mT);
          if (ind==1) {neg_list[4,T[i,j]] <- proc[x,y];}
        } else {
          neg_pos <- neg_pos+1;
          if(neg_pos>neg_free) {neg_free <- neg_free + 100000; neg_list[1,neg_free]<- 0}
          neg_list[,neg_pos] <- Re(c(Tt,i,j,proc[x,y]));
          T[i,j] <- neg_pos;
        }
      }
    }
  }
  
#### Output
  T <- T*time.step
  proc[proc==0] <- NA
  out <- list(seeds = seeds, boosts = boosts, res = res, speeds = n, T = T, proc = proc)
  class(out) <- 'FastMarch'
  
  return(out)
}



##############################
## Compile Code (for speed) ##
##############################
require(compiler)
enableJIT(3)
MFM <- cmpfun(ModifiedFastMarch)



##########################################
## Quick Plot for Fast Marching Results ##
#########################################
plot.FastMarch <- function(FM) {
  
  aux <- unique(c(FM$proc))
  aux <- aux[which(is.na(aux)==F)]
  aux <- sort(aux)
  lev <- diff(aux)/2+aux[1:NROW(aux)-1]
  
  image(FM$T, axes=F); box()
  contour(FM$proc, levels=lev, add=T, drawlabels = F)
  
}
