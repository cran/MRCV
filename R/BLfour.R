#####################################################################################
#                      Functions for Bilder & Loughin (2004)                        #
#####################################################################################

# item.response.table() 
# A function that converts a raw data set to an item response table 
# or data frame (create.dataframe = TRUE)

item.response.table<-function(data, I, J, K=NULL, create.dataframe = FALSE) {
  if (is.numeric(K)) {
    if (K==0) {
      K<-NULL
    }
  }
  nvars<-2+is.numeric(K)
  nrows<-(2^nvars)*I*max(1,J)*max(1,K)
  if (nvars==2) {
    wy.count<-as.data.frame(matrix(data = NA, nrow = nrows, ncol = (2*nvars+1)))
    counter<-0
    for (i in 1:I) {
      for (j in 1:J) {
        counter<-counter+2^nvars
        cell.count<-table(data[,i],data[,I+j])
        wy.count[(counter-3),]<-c(names(data)[i], names(data)[(I+j)], 0, 0, 
                                  cell.count[1,1])
        wy.count[(counter-2),]<-c(names(data)[i], names(data)[(I+j)], 0, 1, 
                                  cell.count[1,2])
        wy.count[(counter-1),]<-c(names(data)[i], names(data)[(I+j)], 1, 0, 
                                  cell.count[2,1])
        wy.count[(counter),]  <-c(names(data)[i], names(data)[(I+j)], 1, 1, 
                                  cell.count[2,2])
      }
    }
    wy.count[,1:nvars]<-lapply(wy.count[,1:nvars], factor)
    wy.count[,(nvars+1):(2*nvars+1)]<-lapply(wy.count[,(nvars+1):(2*nvars+1)], 
                                             as.numeric)
    colnames(wy.count) <- c("W", "Y", "wi", "yj", "count")
    ir.table<-tabular(count*(mean)*W*Heading()*Factor(wi, wi, c(0,1))
                      ~Heading()*Y*Heading()*Factor(yj, yj, c(0,1)), data=wy.count, 
                      suppressLabels=3)
    output<-ir.table
    #Create a summary data set
    if (create.dataframe) {
      output<-wy.count
    }
  }
  if (nvars==3) {
    wyz.count<-as.data.frame(matrix(data = NA, nrow = nrows, ncol = (2*nvars+1)))
    counter<-0
    for (i in 1:I) {
      for (j in 1:J) {
        for (k in 1:K) {
          counter<-counter+2^nvars
          cell.count<-table(data[,i],data[,(I+j)],data[,(I+J+k)])
          wyz.count[(counter-7),]<-c(names(data)[i],names(data)[(I+j)],
                                     names(data)[(I+J+k)],0,0,0,cell.count[1,1,1])
          wyz.count[(counter-6),]<-c(names(data)[i],names(data)[(I+j)],
                                     names(data)[(I+J+k)],0,0,1,cell.count[1,1,2])
          wyz.count[(counter-5),]<-c(names(data)[i],names(data)[(I+j)],
                                     names(data)[(I+J+k)],0,1,0,cell.count[1,2,1])
          wyz.count[(counter-4),]<-c(names(data)[i],names(data)[(I+j)],
                                     names(data)[(I+J+k)],0,1,1,cell.count[1,2,2])
          wyz.count[(counter-3),]<-c(names(data)[i],names(data)[(I+j)],
                                     names(data)[(I+J+k)],1,0,0,cell.count[2,1,1])
          wyz.count[(counter-2),]<-c(names(data)[i],names(data)[(I+j)],
                                     names(data)[(I+J+k)],1,0,1,cell.count[2,1,2])
          wyz.count[(counter-1),]<-c(names(data)[i],names(data)[(I+j)],
                                     names(data)[(I+J+k)],1,1,0,cell.count[2,2,1])
          wyz.count[(counter),]  <-c(names(data)[i],names(data)[(I+j)],
                                     names(data)[(I+J+k)],1,1,1,cell.count[2,2,2])
        }
      }
    }
    wyz.count[,1:nvars]<-lapply(wyz.count[,1:nvars], factor)
    wyz.count[,(nvars+1):(2*nvars+1)]<-lapply(wyz.count[,(nvars+1):(2*nvars+1)], 
                                              as.numeric)
    colnames(wyz.count) <- c("W", "Y", "Z", "wi", "yj", "zk", "count")
    ir.table<-vector("list", K*2)
    ir.table.names<-matrix(data=NA, nrow=1, ncol=K*2)
    counter<-0
    for (k in 1:K) {
      for (c in 0:1) {
        counter<-counter+1
        ir.table[counter]<-list(tabular(count*(mean)*W*Heading()*Factor(wi, wi, c(0,1))
                          ~Heading()*Y*Heading()*Factor(yj, yj, c(0,1)), 
                          data=wyz.count[((wyz.count[,3]==names(data)[I+J+k])&(wyz.count[,6]==c)),], 
                          suppressLabels=3))
        ir.table.names[1,counter]<-paste(names(data[(I+J+k)]), "=", c)
      }
    }
    names(ir.table)<-ir.table.names
    output<-ir.table
    #Create a summary data set
    if (create.dataframe) {
      output<-wyz.count
    }
  }
  output
}

#####################################################################################

# marginal.table() 
# A function that converts a raw data set to a marginal table

marginal.table<-function(data, I, J, K=NULL) {
  if (is.numeric(K)) {
    if (K==0) {
      K<-NULL
    }
  }
  n<-nrow(data)
  nvars<-2+is.numeric(K)
  nrows<-I*max(1,J)*max(1,K)
  if (nvars==2) {
    wy.count<-as.data.frame(matrix(data = NA, nrow = nrows, ncol = (nvars+1)))
    counter<-0
    for (i in 1:I) {
      for (j in 1:J) {
        counter<-counter+1
        cell.count<-table(data[,i],data[,(I+j)])
        wy.count[counter,]<-c(names(data)[i],names(data)[(I+j)],cell.count[2,2]) 
      }
    }
    wy.count[,1:(nvars)]<-lapply(wy.count[,1:(nvars)], factor)
    wy.count[,(nvars+1)]<-as.numeric(wy.count[,(nvars+1)])
    colnames(wy.count) <- c("W", "Y", "count")
    calc.percent <- function(x) round(100*mean(x)/n, 2)
    marginal.table<-tabular(count*W~Heading()*Y*(Heading("count")*mean + 
                            Heading("%")*calc.percent), data=wy.count, 
                            suppressLabels=2)
  }
  if (nvars==3) {
    wyz.count<-as.data.frame(matrix(data = NA, nrow = nrows, ncol = (nvars+1)))
    counter<-0
    for (i in 1:I) {
      for (j in 1:J) {
        for (k in 1:K) {
        counter<-counter+1
        cell.count<-table(data[,i],data[,(I+j)],data[,(I+J+k)])
        wyz.count[counter,]<-c(names(data)[i],names(data)[(I+j)],names(data)[(I+J+k)],
                               cell.count[2,2,2]) 
        }
      }
    }
    wyz.count[,1:(nvars)]<-lapply(wyz.count[,1:(nvars)], factor)
    wyz.count[,(nvars+1)]<-as.numeric(wyz.count[,(nvars+1)])
    colnames(wyz.count) <- c("W", "Y", "Z", "count")
    calc.percent <- function(x) round(100*mean(x)/n, 2)
    marginal.table<-vector("list", K)
    for (k in 1:K) {
      marginal.table[k]<-list(tabular(count*W~Heading()*Y*(Heading("count")*mean + 
                              Heading("%")*calc.percent), 
                              data=wyz.count[(wyz.count[,3]==names(data)[I+J+k]),], 
                              suppressLabels=2))
    }
    names(marginal.table)<-names(data[(I+J+1):(I+J+K)])
  }
  c<-NULL
  marginal.table
}

#####################################################################################

# check.zero()
# Function that adds a small constant to 0 cell counts

check.zero<-function(value, add.constant = .5) {
  if (value==0) newvalue<-add.constant else newvalue<-value
  newvalue
}

#####################################################################################

# SPMI.stat()
# A function that calculates X^2_S and X^2_S.ij

SPMI.stat<-function(data, I, J, summary.data = FALSE, add.constant = .5) {
  #For inputted data set that is a summary file
  if (summary.data) {
    data[,1:2]<-lapply(data[,1:2], factor)
    X.sq.S.ij<-matrix(data = NA, nrow = I, ncol = J)
    counter<-0
    for (i in 1:I) {
      for (j in 1:J) {
        table.11<-data[((data[,1]==levels(data[,1])[i])& 
                      (data[,2]==levels(data[,2])[j])&(data[,3]==0)&(data[,4]==0)), 5]
        table.12<-data[((data[,1]==levels(data[,1])[i])& 
                      (data[,2]==levels(data[,2])[j])&(data[,3]==0)&(data[,4]==1)), 5]
        table.21<-data[((data[,1]==levels(data[,1])[i])& 
                      (data[,2]==levels(data[,2])[j])&(data[,3]==1)&(data[,4]==0)), 5]
        table.22<-data[((data[,1]==levels(data[,1])[i])& 
                      (data[,2]==levels(data[,2])[j])&(data[,3]==1)&(data[,4]==1)), 5]
        n.table<-matrix(data = c(table.11, table.21, table.12, table.22), 
                      nrow = 2, ncol = 2)
        #Only calculate statistics for valid 2x2 tables (no 00 in row or col)
        if (min(c(rowSums(n.table), colSums(n.table))) > 0) {
          counter<-counter+1
          #Add .5 to 0 cell counts
          n.table<-apply(X = n.table, MARGIN = c(1,2), FUN = check.zero, 
                         add.constant = add.constant)
          options(warn = -1) 
          X.sq.S.ij[i,j]<-chisq.test(n.table,correct=F)$statistic
          options(warn = 0) 
        }
      }
    }
    rownames(X.sq.S.ij)<-levels(data[,1])
    colnames(X.sq.S.ij)<-levels(data[,2])
    output<-list(X.sq.S = sum(X.sq.S.ij), X.sq.S.ij = X.sq.S.ij, valid.margins = counter)
  }  
  
  #For inputted data set that is a raw file
  if (1-summary.data) {
    X.sq.S.ij<-matrix(data = NA, nrow = I, ncol = J)
    counter<-0
    for (i in 1:I) {
      for (j in 1:J) {
        #Only calculate statistics for valid tables (correct dim = 2x2)
        if (sum(dim(table(data[,i],data[,(I+j)])))==4){
          counter<-counter+1
          options(warn = -1) 
          n.table<-table(data[,i],data[,(I+j)])
          #Add .5 to 0 cell counts
          n.table<-apply(X = n.table, MARGIN = c(1,2), FUN = check.zero, 
                         add.constant = add.constant)
          X.sq.S.ij[i, j]<-chisq.test(n.table,correct=F)$statistic
          options(warn = 0) 
        }
      }
    }
    rownames(X.sq.S.ij)<-names(data)[1:I]
    colnames(X.sq.S.ij)<-names(data)[(I+1):(I+J)]
    output<-list(X.sq.S = sum(X.sq.S.ij), X.sq.S.ij = X.sq.S.ij, valid.margins = 
                counter)
  }
  c<-NULL
  output
}

#####################################################################################

# Function used in Bonferroni adjustments (p-value should not be greater than 1)

check.min<-function(value) {
  min(value, 1.00000)
}

#####################################################################################

# Function used to control display of output

print.SPMI<-function(x, ...) {
  options(scipen=10)
  type<-names(x)
  data<-x$general$data
  I<-x$general$I
  J<-x$general$J
  summary.data<-x$general$summary.data
  X.sq.S<-round(x$general$X.sq.S, 5)
  X.sq.S.ij<-round(x$general$X.sq.S.ij, 5)
  rownames(X.sq.S.ij)<-names(data)[1:I]
  colnames(X.sq.S.ij)<-names(data)[(I+1):(I+J)]
  if (summary.data) {
    rownames(X.sq.S.ij)<-levels(data[,1])
    colnames(X.sq.S.ij)<-levels(data[,2]) 
  }
  cat("Unadjusted Pearson Chi-Square Tests for Independence:", "\n")
  cat("X^2_S =", X.sq.S, "\n")
  cat("X^2_S.ij =", "\n") 
  print.default(X.sq.S.ij)
  cat("\n")
  if (any(type=="boot")){
    B.discard<-x$boot$B.discard
    B.use<-x$boot$B.use
    p.value.boot<-as.name(paste("=", round(x$boot$p.value.boot, 5))) 
    if (x$boot$p.value.boot == 0) {
      p.value.boot<-as.name(paste("<", round(1/B.use, 5)))
    }
    p.combo.prod<-as.name(paste("=", round(x$boot$p.combo.prod.boot, 5)))
    if (x$boot$p.combo.prod.boot == 0) {
      p.combo.prod<-as.name(paste("<", round(1/B.use, 5)))
    }
    p.combo.min<-as.name(paste("=", round(x$boot$p.combo.min.boot, 5)))
    if (x$boot$p.combo.min.boot ==0) {
      p.combo.min<-as.name(paste("<", round(1/B.use, 5)))
    }
    cat("Bootstrap Results:", "\n")
    if (B.discard > 0) {
      cat(B.discard, "resamples were removed from the analysis due to")  
      cat(" not having all rows or columns represented in a 2x2 table", "\n")
    }
    cat("Final results based on", B.use, "resamples", "\n")
    cat("p.boot", p.value.boot, "\n")
    cat("p.combo.prod", p.combo.prod, "\n")
    cat("p.combo.min", p.combo.min, "\n", "\n")
  }
  if (any(type=="rs2")){
    X.sq.S.rs2<-round(x$rs2$X.sq.S.rs2, 5)
    df.rs2<-round(x$rs2$df.rs2, 5)
    p.value.rs2<-as.name(paste("=", round(x$rs2$p.value.rs2, 5)))
    if (round(x$rs2$p.value.rs2, 5) < .00001) {
      p.value.rs2<-as.name(paste("<", .00001))
    }
    cat("Second-Order Rao-Scott Adjusted Results:", "\n")
    cat("X^2_S.adj =", X.sq.S.rs2, "\n")
    cat("df.adj =", df.rs2, "\n")
    cat("p.adj", p.value.rs2, "\n", "\n")
  }
  if (any(type=="bon")){
    p.value.bon<-as.name(paste("=", round(x$bon$p.value.bon, 5)))
    if (round(x$bon$p.value.bon, 5) < .00001) {
      p.value.bon<-as.name(paste("<", .00001))
    }
    X.sq.S.ij.p.bon<-format(round(x$bon$X.sq.S.ij.p.bon,5), 
                           digits = 5)
    rownames(X.sq.S.ij.p.bon)<-names(data)[1:I]
    colnames(X.sq.S.ij.p.bon)<-names(data)[(I+1):(I+J)] 
    if (summary.data) {
      rownames(X.sq.S.ij.p.bon)<-levels(data[,1])
      colnames(X.sq.S.ij.p.bon)<-levels(data[,2]) 
    }
    cat("Bonferroni Adjusted Results:", "\n")
    cat("p.adj", p.value.bon, "\n")    
    cat("p.ij.adj =", "\n")
    print.default(X.sq.S.ij.p.bon, quote=FALSE)
    cat("\n")
  }
  options(scipen=0)
  invisible(x)
}

#####################################################################################

# SPMI.test()
# A function that performs bootstrap testing, the second-order Rao-Scott approach, 
#   and Bonferroni adjustments

SPMI.test<-function(data, I, J, type = "all", B = 1999, B.max = B, 
                    summary.data = FALSE, add.constant = .5, plot.hist = FALSE, 
                    print.status = TRUE) {
  #Summary data can only be used with the Bonferroni adjustment
  if (summary.data) { 
    type<-"bon"
  }
  n<-nrow(data)
  #Observed statistics 
  observed<-SPMI.stat(data = data, I = I, J = J, 
                      summary.data = summary.data, add.constant = add.constant)
  p.value.ij<-1 - pchisq(q = observed$X.sq.S.ij, df = 1)
  p.value.min<-min(p.value.ij)
  p.value.prod<-prod(p.value.ij)    
  
  #Bootstrap
  if (type == "boot" | type == "all") {
    X.sq.S.star<-numeric(length(B.max))
    p.value.b.min<-numeric(length(B.max))
    p.value.b.prod<-numeric(length(B.max))
    X.sq.S.ij.star<-matrix(data=NA, nrow = I*J, ncol = B.max)
    discard<-0 #Count number of resamples with incorrect dimensions
    counter<-0
    b<-0
    # create progress bar
    if (print.status) {
      cat("Bootstrap Progress:", "\n")
      pb <- txtProgressBar(min = 0, max = B, style = 3)
    }
    while(((counter>=B)+(b>=B.max)) < 1) {
      b<-b+1
      # update progress bar
      if (print.status) {
        setTxtProgressBar(pb, counter+1)
      }
      #Resample W_s and Y_s independently
      W<-sample(x = 1:n, size = n, replace = TRUE) 
      Y<-sample(x = 1:n, size = n, replace = TRUE)
      data.star<-cbind(data[W, 1:I], data[Y, (I+1):(I+J)])
      stat.star<-SPMI.stat(data = data.star, I = I, J = J, 
                           summary.data = FALSE, add.constant = add.constant)
      discard<-discard+(stat.star$valid.margins<I*J) #Discard if any table < 2x2
      if (stat.star$valid.margins==I*J){ 
        counter<-counter+1
        X.sq.S.star[counter]<-stat.star$X.sq.S
        X.sq.S.ij.star[,counter]<-stat.star$X.sq.S.ij
        p.value.ij<-1 - pchisq(q = stat.star$X.sq.S.ij, df = 1)
        p.value.b.min[counter]<-min(p.value.ij)
        p.value.b.prod[counter]<-prod(p.value.ij)
      }
    }
    if (print.status) {
      close(pb)
    }
    B.use<-min(B,(B.max-discard)) #Only use desired number of resamples
    X.sq.S.star<-X.sq.S.star[1:B.use]
    X.sq.S.ij.star<-X.sq.S.ij.star[,1:B.use]
    p.value.b.min<-p.value.b.min[1:B.use]
    p.value.b.prod<-p.value.b.prod[1:B.use]
    p.value.boot<-mean(X.sq.S.star>=observed$X.sq.S)
    p.combo.min<-list(p = p.value.min, p.star = p.value.b.min, 
                      overall.p = mean(p.value.b.min<=p.value.min))
    p.combo.prod<-list(p = p.value.prod, p.star = p.value.b.prod, 
                       overall.p = mean(p.value.b.prod<=p.value.prod))
    #Histograms
    if (plot.hist) {
      par(mfrow = c(2,2))
      hist(x = X.sq.S.star, main = "Test using sum statistic", xlab = 
           expression(X[S]^{"2*"}), xlim = c(min(X.sq.S.star, observed$X.sq.S), 
           max(X.sq.S.star, observed$X.sq.S)))
      abline(v = observed$X.sq.S, col = "darkgreen", lwd = 5)
      hist(x = p.value.b.prod, main = "Test using product of p-values", xlab = 
           expression(tilde(p)[prod]^{"*"}), xlim = c(min(p.value.b.prod, 
           p.value.prod), max(p.value.b.prod, p.value.prod)))
      abline(v = p.value.prod, col = "darkgreen", lwd = 5)
      hist(x = p.value.b.min, main = "Test using minimum of p-values", xlab = 
           expression(tilde(p)[min]^{"*"}), xlim = c(min(p.value.b.min, 
           p.value.min), max(p.value.b.min, p.value.min)))
      abline(v = p.value.min, col = "darkgreen", lwd = 5) 
      par(mfrow = c(1,1))
    }
    output<-list(general = list(data = data, I = I, J = J, summary.data = 
                 summary.data, X.sq.S = observed$X.sq.S, X.sq.S.ij = 
                 observed$X.sq.S.ij), boot = list(B.use = B.use, 
                 B.discard = discard, p.value.boot = p.value.boot, 
                 p.combo.min.boot = p.combo.min$overall.p, p.combo.prod.boot = 
                 p.combo.prod$overall.p, X.sq.S.star = X.sq.S.star, X.sq.S.ij.star = 
                 X.sq.S.ij.star, p.combo.min.star = p.combo.min$p.star, 
                 p.combo.prod.star = p.combo.prod$p.star))
    output.boot<-list(B.use = B.use, B.discard = discard, 
                 p.value.boot = p.value.boot, p.combo.min.boot = 
                 p.combo.min$overall.p, p.combo.prod.boot = 
                 p.combo.prod$overall.p, X.sq.S.star = X.sq.S.star, 
                 X.sq.S.ij.star = X.sq.S.ij.star, p.combo.min.star = 
                 p.combo.min$p.star, p.combo.prod.star = p.combo.prod$p.star)
  }
  
  #Second-order Rao-Scott adjustment
  if (type == "rs2" | type == "all") {
    #Appendix A calculations
    W.counts<-as.data.frame(table(data[,1:I])) #Get all possible combos of W's
    cols <- c(1:I) #Need to order W's in ascending order starting with W1 
    W.counts.ord<-W.counts[do.call("order", as.data.frame(W.counts[,cols])), ]
    Y.counts<-as.data.frame(table(data[,(I+1):(I+J)])) #All possible Y's
    cols <- c(1:J) #Need to order Y's in ascending order starting with Y1
    Y.counts.ord<-Y.counts[do.call("order", as.data.frame(Y.counts[,cols])), ]
    n.counts<-as.data.frame(table(data)) #Get all possible combos of W's and Y's
    cols <- c(1:(I+J)) #Need to order W's and Y's in ascending order W1, W2, ...
    n.counts.ord<-n.counts[do.call("order", as.data.frame(n.counts[,cols])), ]
    G<-t(data.matrix(W.counts.ord[,1:I])-1) #rx2^r matrix 
    H<-t(data.matrix(Y.counts.ord[,1:J])-1) #cx2^c matrix 
    tau<-n.counts.ord[,(I+J+1)]/n  #Vector of multinomial probabilities
    m.row<-G%*%W.counts.ord[,(I+1)]  #Vector of W marginal counts
    m.col<-H%*%Y.counts.ord[,(J+1)]  #Vector of Y marginal counts
    GH<-kronecker(G,H)  
    m<-GH%*%n.counts.ord[,(I+J+1)]  #Marginal table counts (row 1, row 2, ...)
    pi<-m/n   #pi_ij
    pi.row<-m.row/n #pi_i.
    pi.col<-m.col/n #pi_.j  
    j.2r<-matrix(data = 1, nrow = 2^I, ncol = 1) #2^r vector of 1's
    i.2r<-diag(2^I) #(2^r)x(2^r) identity matrix
    j.2c<-matrix(data = 1, nrow = 2^J, ncol = 1) #2^c vector of 1's
    i.2c<-diag(2^J) #(2^c)x(2^c) identity matrix
    G.ij<-G%*%kronecker(i.2r,t(j.2c))
    H.ji<-H%*%kronecker(t(j.2r),i.2c)
    F<-GH - kronecker(pi.row,H.ji) - kronecker(G.ij,pi.col)
    mult.cov<-diag(tau) - tau%*%t(tau)
    sigma<-F%*%mult.cov%*%t(F)
    D<-diag(as.vector(kronecker(pi.row,pi.col)*kronecker(1-pi.row,1-pi.col)))
    Di.sigma<-solve(D)%*%sigma
    Di.sigma.eigen<-Re(eigen(Di.sigma)$values) #Only use real part of eigenvalues
    X.sq.S.rs2<-I*J*observed$X.sq.S/sum(Di.sigma.eigen^2) 
    df.rs2<-I^2*J^2/sum(Di.sigma.eigen^2) 
    X.sq.S.p.value.rs2<-1-pchisq(q = X.sq.S.rs2, df = df.rs2) 
    output<-list(general = list(data = data, I = I, J = J, summary.data = 
                 summary.data, X.sq.S = observed$X.sq.S, X.sq.S.ij = 
                 observed$X.sq.S.ij), rs2 = list(X.sq.S.rs2 = X.sq.S.rs2, df.rs2 = 
                 df.rs2, p.value.rs2 = X.sq.S.p.value.rs2))
    output.rs2<-list(X.sq.S.rs2 = X.sq.S.rs2, df.rs2 = df.rs2, 
                 p.value.rs2 = X.sq.S.p.value.rs2)
  }
  
  #Bonferroni
  if (type == "bon" | type == "all") {   
    p.value.bon<-min(p.value.min*I*J,1) #p-value should not be greater than 1
    X.sq.S.ij.p.bon<-apply(X = I*J*(1 - pchisq(q = observed$X.sq.S.ij, 
                           df = 1)), MARGIN = c(1,2), FUN = check.min)
    output<-list(general = list(data = data, I = I, J = J, summary.data = 
                 summary.data, X.sq.S = observed$X.sq.S, X.sq.S.ij = 
                 observed$X.sq.S.ij), bon = list(p.value.bon = p.value.bon, 
                 X.sq.S.ij.p.bon = X.sq.S.ij.p.bon))
    output.bon<-list(p.value.bon = p.value.bon, X.sq.S.ij.p.bon = X.sq.S.ij.p.bon) 
    if (summary.data) {
      output<-list(general = list(data = data, I = I, J = J, summary.data = 
                   summary.data, X.sq.S = observed$X.sq.S, X.sq.S.ij = 
                   observed$X.sq.S.ij), bon = list(p.value.bon = p.value.bon, 
                   X.sq.S.ij.p.bon = X.sq.S.ij.p.bon))
    }
  }
  if (type == "all") {
    output<-list(general = list(data = data, I = I, J = J, summary.data = 
                 summary.data, X.sq.S = observed$X.sq.S, X.sq.S.ij = 
                 observed$X.sq.S.ij), boot = output.boot, rs2 = output.rs2, 
                 bon = output.bon)
  }
  class(output)<-"SPMI"
  output
}

#####################################################################################