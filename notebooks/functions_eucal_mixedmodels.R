#############################################################################
## Extracting objects from a mixed model fitted by lme function #############                                               ##
#############################################################################

## Packages
{library(ggplot2);library(dplyr);library(lme4);library(nlme);
  library(lmerTest);library(cAIC4);library(splines);library(gamlss);
  library(gridExtra);library(hnp);library(varTestnlme);library(car);
  library(Rcpp);library(MASS);library(Matrix);library(tidyverse);
  library(merDeriv);library(qqplotr)}

##########################################################################

sqrt.matrix <- function(mat) 
{              
  mat <- as.matrix(mat)  
  singular_dec <- svd(mat)
  U <- singular_dec$u
  V <- singular_dec$v
  D <- diag(singular_dec$d)
  sqrtmatrix <- U %*% sqrt(D) %*% t(V)
}

###########################################################################

extract.lmeDesign2 <- function(m)
{
  start.level = 1
  data <- getData(m)
  grps <- nlme::getGroups(m)
  n <- length(grps)
  X <- list()
  grp.dims <- m$dims$ncol
  Zt <- model.matrix(m$modelStruct$reStruct, data)
  cov <- as.matrix(m$modelStruct$reStruct)
  i.col <- 1
  n.levels <- length(m$groups)
  Z <- matrix(0, n, 0)
  if (start.level <= n.levels) {
    for (i in 1:(n.levels - start.level + 1)) {
      if (length(levels(m$groups[[n.levels - i + 1]])) != 1)
      {
        X[[1]] <- model.matrix(~m$groups[[n.levels - i +
                                            1]] - 1, 
                               contrasts.arg = c("contr.treatment",
                                                 "contr.treatment"))
      }
      else X[[1]] <- matrix(1, n, 1)
      X[[2]] <- as.matrix(Zt[, i.col:(i.col + grp.dims[i] -
                                        1)])
      i.col <- i.col + grp.dims[i]
      Z <- cbind(mgcv::tensor.prod.model.matrix(X),Z)
    }
    Vr <- matrix(0, ncol(Z), ncol(Z))
    start <- 1
    for (i in 1:(n.levels - start.level + 1)) {
      k <- n.levels - i + 1
      for (j in 1:m$dims$ngrps[i]) {
        stop <- start + ncol(cov[[k]]) - 1
        Vr[ncol(Z) + 1 - (stop:start),ncol(Z) + 1 - (stop:start)] <- cov[[k]]
        start <- stop + 1
      }
    }
  }
  X <- if (class(m$call$fixed) == "name" &&  !is.null(m$data$X)) {
    m$data$X
  } else   {
    model.matrix(formula(eval(m$call$fixed)),data)
  }
  y <- as.vector(matrix(m$residuals, ncol = NCOL(m$residuals))[,NCOL(m$residuals)] + 
                   matrix(m$fitted, ncol = NCOL(m$fitted))[,NCOL(m$fitted)])
  return(list(
    Vr = Vr,                                                                 
    X = X,
    Z = Z,
    sigmasq = m$sigma ^ 2,
    lambda = unique(diag(Vr)),
    y = y,
    k = n.levels
  )
  )
}


##################################################################

res_lcr <- function(fit, graph = c("envelope", "hist", "mahalanobis distance"))
{
  { 
    data.fit <- suppressWarnings(extract.lmeDesign2(fit))
    data <-    getData(fit)
    y <- data.fit$y
    X <- data.fit$X
    N <- length(y)                                                               
    id <-  sort(as.numeric(getGroups(fit, level = 1)), index.return = TRUE)$x     
    subject <- as.numeric(unique(id))
    n <- length(as.numeric(names(table(id))))                                    
    vecni <- (table(id))                                                         
    p <- ncol(X)                                                                 
    n.levels <- length(fit$groups)                                               
    start.level <- 1
    Cgrps <- nlme::getGroups(fit, level = start.level)                           
    CCind <- levels((Cgrps))  
    sigma2 <- fit$sigma^2
  }
  obs <- numeric()
  
  for (i in 1:n)
  {
    obs <- append(obs,1:vecni[i])                                               
  }
  if (n.levels > 1) { 
    lZi <- list()
    lgi <- list()
    numrow <- numeric()
    
    mgroups <- fit$groups      
    for (n in 1:length(CCind)) {
      dgi <- data.frame(as.matrix(mgroups[mgroups == CCind[n], ]))
      nrowzi <- dim(dgi)[1]
      ncolzi <- 0
      girep <- as.numeric(length(levels(dgi[,1])))
      for (k in 2:n.levels) {
        girep <- c(girep,as.numeric(length(levels(dgi[,k]))))
      }
      for (k in 1:n.levels) {
        ncolzi <- ncolzi + as.numeric(length(levels(dgi[,k])))
      }
      auxi <- as.vector(table(dgi[,1]))
      for (i in 2:n.levels) {
        auxi <- c(auxi,as.vector(table(dgi[,i])))
      }
      l <- 1
      Zi <- matrix(0,nrowzi,ncolzi)
      for (j in 1:ncolzi) {
        Zi[l:(l + auxi[j] - 1),j] <- rep(1,auxi[j]) 
        l <- l + auxi[j]
        if (l == (nrowzi + 1)) l <- 1
      }
      
      lZi[[n]] <- Zi
      
      numrow[n] <- dim(Zi)[1]
      
      comp.var <- as.matrix(fit1$modelStruct$reStruct)
      auxg <- rep(as.numeric(comp.var[1])*sigma2,girep[1])
      for (i in 2:length(girep)) {
        auxg <- c(auxg,rep(as.numeric(comp.var[i])*sigma2,girep[i]))
      }
      lgi[[n]] <- diag(auxg)
    }
    q <- dim(lgi[[1]])[1]                     
    for (h in 2:length(CCind)) {
      q <- c(q,dim(lgi[[h]])[1])
    }
    Z <- lZi[[1]]
    for (k in 2:length(CCind)) {
      Z <- bdiag(Z,(lZi[[k]]))
    }
    Z <- as.matrix(Z)
    nrowZi <- lZi[[1]]                        
    for (h in 2:length(CCind)) {
      nrowZi <- c(nrowZi,dim(lZi[[h]])[1])
    }
    
    Gam <- lgi[[1]]
    for (k in 2:length(CCind)) {
      Gam <- bdiag(Gam,(lgi[[k]]))
    }
    Gam <- as.matrix(Gam)
  }else{
    mataux <- model.matrix(fit$modelStruct$reStruct,data)
    mataux <- as.data.frame(cbind(mataux,id))
    lZi <- list()
    lgi <- list()
    
    for (i in (as.numeric(unique(id)))) { 
      lZi[[i]] <- as.matrix((subset(split(mataux,id == i,
                                          drop = T)$`TRUE`,select = -id)))          
      lgi[[i]] <- getVarCov(fit,type = "random.effects")
    }
    Z <- as.matrix(bdiag(lZi))
    g <- getVarCov(fit,type = "random.effects")
    q <- dim(g)[1]                                                           
    Gam <- as.matrix(kronecker(diag(length(as.numeric(unique(id)))),g))
  }
  
  if (n.levels > 1) {   
    if (!inherits(fit, "lme")) 
      stop("object does not appear to be of class lme")
    grps <- nlme::getGroups(fit)
    n <- length(grps)                                                                     
    n.levels <- length(fit$groups)                                                         
    if (is.null(fit$modelStruct$corStruct)) 
      n.corlevels <- 0
    else n.corlevels <- length(all.vars(nlme::getGroupsFormula(fit$modelStruct$corStruct))) 
    if (n.levels < n.corlevels) {
      getGroupsFormula(fit$modelStruct$corStruct)
      vnames <- all.vars(nlme::getGroupsFormula(fit$modelStruct$corStruct))
      lab <- paste(eval(parse(text = vnames[1]), envir = fit$data))
      if (length(vnames) > 1) 
        for (i in 2:length(vnames)) {
          lab <- paste(lab, "/", eval(parse(text = vnames[i]), 
                                      envir = fit$data), sep = "")
        }
      grps <- factor(lab)
    }
    if (n.levels >= start.level || n.corlevels >= start.level) {
      if (n.levels >= start.level) 
        Cgrps <- nlme::getGroups(fit, level = start.level)                          
      else Cgrps <- grps
      Cind <- sort(as.numeric(Cgrps), index.return = TRUE)$ix                       
      rCind <- 1:n 
      rCind[Cind] <- 1:n
      Clevel <- levels(Cgrps)                                                      
      n.cg <- length(Clevel)                                                         
      size.cg <- array(0, n.cg)
      for (i in 1:n.cg) size.cg[i] <- sum(Cgrps == Clevel[i])  
    }
    else {
      n.cg <- 1
      Cind <- 1:n
    }
    if (is.null(fit$modelStruct$varStruct)) 
      w <- rep(fit$sigma, n)
    else {
      w <- 1/nlme::varWeights(fit$modelStruct$varStruct)
      group.name <- names(fit$groups)
      order.txt <- paste("ind<-order(data[[\"", group.name[1], 
                         "\"]]", sep = "")
      if (length(fit$groups) > 1) 
        for (i in 2:length(fit$groups)) order.txt <- paste(order.txt, 
                                                           ",data[[\"", group.name[i], "\"]]", sep = "")
      order.txt <- paste(order.txt, ")")
      eval(parse(text = order.txt))
      w[ind] <- w
      w <- w * fit$sigma
    }
    w <- w[Cind]
    if (is.null(fit$modelStruct$corStruct)) 
      lR <- array(1, n)
    else {
      c.m <- nlme::corMatrix(fit$modelStruct$corStruct)
      if (!is.list(c.m)) {
        lR <- c.m
        lR <- lR[Cind, ]
        lR <- lR[, Cind]
      }
      else {
        lR <- list()
        ind <- list()
        for (i in 1:n.cg) {
          lR[[i]] <- matrix(0, size.cg[i], size.cg[i])
          ind[[i]] <- 1:size.cg[i]
        }
        Roff <- cumsum(c(1, size.cg))
        gr.name <- names(c.m)
        n.g <- length(c.m)
        j0 <- rep(1, n.cg)
        ii <- 1:n
        for (i in 1:n.g) {
          Clev <- unique(Cgrps[grps == gr.name[i]])
          if (length(Clev) > 1) 
            stop("inner groupings not nested in outer!!")
          k <- (1:n.cg)[Clevel == Clev]
          j1 <- j0[k] + nrow(c.m[[i]]) - 1
          lR[[k]][j0[k]:j1, j0[k]:j1] <- c.m[[i]]
          ind1 <- ii[grps == gr.name[i]]
          ind2 <- rCind[ind1]
          ind[[k]][j0[k]:j1] <- ind2 - Roff[k] + 1
          j0[k] <- j1 + 1
        }
        for (k in 1:n.cg) {
          lR[[k]][ind[[k]], ] <- lR[[k]]
          lR[[k]][, ind[[k]]] <- lR[[k]]
        }
      }
    }
    if (is.list(lR)) {
      for (i in 1:n.cg) {
        wi <- w[Roff[i]:(Roff[i] + size.cg[i] - 1)]
        lR[[i]] <- as.vector(wi) * t(as.vector(wi) * lR[[i]]) 
      }
    }
    else if (is.matrix(lR)) {
      lR <- as.vector(w) * t(as.vector(w) * lR)
    }
    else {
      lR <- w^2 * lR
    }
    if (is.list(lR)) {
      R <- lR[[1]]
      for (k in 2:n.cg) {
        R <- bdiag(R,lR[[k]])
      }
      R <- as.matrix(R)
    }
    else{
      R <- diag(lR)
    }
  }else{
    R <- getVarCov(fit,type = "conditional",individual = 1)[[1]]
    for (i in 2:length(as.numeric(unique(id)))) {
      R <- as.matrix(bdiag(R,getVarCov(fit,
                                       type = "conditional",individual = i)[[1]] ) )
    }
  }
  
  {
    V <- (Z %*% Gam %*% t(Z)) + R
    iV <- solve(V)                                                
    varbeta <- solve((t(X) %*% iV %*% X))
    Q <- (iV - iV %*% X %*% (varbeta) %*% t(X) %*% iV ) 
    zq <- t(Z) %*% Q
    norm.frob.ZtQ <- sum(diag(zq %*% t(zq)))
    eblue <- as.vector(fixef(fit))
    eblup <- Gam %*% t(Z) %*% iV %*% (y - X %*% eblue)
    pre <- X %*% eblue                       
    predi <- X %*% eblue + Z %*% eblup         
    resm <- (y - pre)                        
    resc <- (y - predi)  
    var.resm <- V - X %*% solve(t(X) %*% iV %*% X) %*% t(X) 
    var.resc <- R %*% Q %*% R
    ident <- diag(N)
    auxnum <- (R %*% Q %*% Z %*% Gam %*% t(Z) %*% Q %*% R)
    auxden <- R %*% Q %*% R
    CF <- diag(auxnum)/diag(auxden)
    rescp <- resc/sqrt(diag(var.resc))
    R.half <- sqrt.matrix(R)
    auxqn <- eigen((R.half %*% Q %*% R.half), symmetric = T, only.values = FALSE) 
    lt <- sqrt(solve(diag((auxqn$values[1:(N-p)])))) %*% t(auxqn$vectors[1:N,1:(N-p)]) %*% solve(sqrt.matrix(R[1:N,1:N]))
    var.resmcp <- lt %*% var.resc[1:N,1:N] %*% t(lt)
    prive <- (lt %*% resc[1:N] )/sqrt(diag(var.resmcp))
    row_lines <- as.matrix(prive[c(1:5), ],dim = c(5, 1))
    resmcp <- rbind(prive, row_lines)
  }
  
  ####
  
  if (n.levels > 1) {
    aux <- Gam %*% t(Z) %*% Q %*% Z %*% Gam
    qm <- q - 1
    dm <- matrix(0,length(CCind),1)
    gbi <- aux[1:(q[1]),(1:q[1])]
    eblupi <- eblup[1:(q[1]),]
    dmi <- t(eblupi) %*% ginv(gbi) %*% eblupi
    dm[1] <- dmi
    for (j in 2:length(CCind)) {
      gbi <- aux[((j - 1)*q[(j - 1)] + 1 ):(q[j] + q[(j - 1)]),((j - 1)*q[(j - 1)] + 1 ):(q[j] + q[(j - 1)])]
      eblupi <- eblup[((j - 1)*q[(j - 1)] + 1 ):(q[j] + q[(j - 1)]),]
      dmi <- t(eblupi) %*% ginv(gbi) %*% eblupi
      dm[j] <- dmi
    }
  }else{
    aux <- Gam %*% t(Z) %*% Q %*% Z %*% Gam
    qm <- q - 1
    dm <- matrix(0,n,1)
    
    for (j in 1:length(CCind)) 
    {
      if (q == 1)
      {
        gbi <- aux[j,j]
        eblupi <- eblup[(q*j - qm):(q*j)]
        dmi <- t(eblupi) %*% ginv(gbi) %*% eblupi
        dm[j] <- dmi
      }
      else
      {
        gbi <- aux[(q*j - qm):(q*j),(q*j - qm):(q*j)]
        cat(gbi,'\n','\t')
        eblupi <- eblup[(q*j - qm):(q*j)]
        dmi <- t(eblupi) %*% ginv(gbi) %*% eblupi
        dm[j] <- dmi
      }
    }
    
  }
  
  ####
  
  if(graph == "envelope")
  {
    df1 <- gVals(resmcp, "norm", .95)
    
    p = ggplot(df1, aes(x = z, y = y)) + 
      geom_smooth(method = rlm, se = F, color = "red", linewidth = 1.7) +
      annotate("line", x = df1$z, y = df1$uc, color = "red", lty = 2, linewidth=1.7) +
      annotate("line", x = df1$z, y = df1$dc, color = "red", lty = 2, linewidth=1.7) +
      geom_point(size=0.8)+
      labs(x = 'N(0,1) quantiles', y = 'Least Confounded Residuals') + 
      theme(legend.position = "none",
            axis.title = element_text(size = 25),
            axis.text.x = element_text(color = "black", hjust=1),
            axis.text.y = element_text(color = "black", hjust=1),
            axis.text = element_text(size = 25),
            plot.title = element_text(size = 25))
  }
  if(graph == "hist")
  {
    Data_least_conf_res = as.data.frame(resmcp)
    
    p = ggplot(data = Data_least_conf_res, aes(x = V1)) +
        geom_histogram(color="white", bins = 20) +
        labs(x = 'Least Confounded Residuals', y = 'Density') + 
        #coord_cartesian(ylim=c(0,110))+
        theme(legend.position = "none",
              axis.title = element_text(size = 25),
              axis.text.x = element_text(color = "black", hjust=1),
              axis.text.y = element_text(color = "black", hjust=1),
              axis.text = element_text(size = 25),
              plot.title = element_text(size = 25))
  }
  if(graph == "mahalanobis distance")
  {
    QT <- data.frame(QUANTIS = dm)
    p = ggplot(QT, aes(sample = QUANTIS)) +
      stat_qq_band(distribution = "chisq", dparams = list(df = q), size=1.7) +
      stat_qq_point(distribution = "chisq", dparams = list(df = q), size = 1.5) +
      labs(x = 'Chi-squared Quantiles', y = 'Mahalanobis Distance Quantiles') + 
      stat_qq_line(distribution = "chisq", dparams = list(df = q), color = "red") +
      theme(legend.position = "none",
            axis.title = element_text(size = 25),
            axis.text.x = element_text(color = "black", hjust=1),
            axis.text.y = element_text(color = "black", hjust=1),
            axis.text = element_text(size = 25),
            plot.title = element_text(size = 25))
  }
  if(graph == "dispersion")
  {
    Data_least_conf_res = as.data.frame(resmcp)
    fitted_values <- as.data.frame(predict(fit))
    data_aux = data.frame(fitted = as.numeric(fitted_values$`predict(fit)`),
                          resid = as.numeric(Data_least_conf_res$V1))
    
    p = ggplot(aes(fitted, resid), data = data_aux) +
      geom_point(size=1.5) +
      coord_cartesian(ylim=c(-5,5)) +
      geom_hline(yintercept = 0) +
      labs(x = "Fitted Values", y = "Least Confounded Residuals") +  
      theme(axis.title = element_text(size = 25),
            axis.text.x = element_text(color = "black", hjust=1),
            axis.text.y = element_text(color = "black", hjust=1),
            axis.text = element_text(size = 25),
            plot.title = element_text(size = 25)) 
  }
  
  return(list(p = p, resmcp = resmcp))
}

##########################################################################################################
## This function constructs QQ plots for normality of random eefects                                    ##
##########################################################################################################
qqPlot2 <- function(x, distribution="norm", ..., ylab=deparse(substitute(x)),
                    xlab=paste(distribution, "quantiles"), main = NULL, 
                    las = par("las"),
                    envelope = .95,  
                    col = palette()[1], 
                    col.lines = palette()[2], lwd = 2, pch = 1, cex = par("cex"),
                    cex.lab = par("cex.lab"), cex.axis = par("cex.axis"), 
                    line = c("quartiles", "robust", "none"), 
                    labels = if (!is.null(names(x))) names(x) else seq(along = x),
                    id.method = "y", 
                    id.n = if (id.method[1] == "identify") Inf else 0,
                    id.cex = 1, id.col=palette()[1], grid = TRUE)
{
  line <- match.arg(line)
  good <- !is.na(x)
  ord <- order(x[good])
  ord.x <- x[good][ord]
  ord.lab <- labels[good][ord]
  q.function <- eval(parse(text = paste("q", distribution, sep = "")))
  d.function <- eval(parse(text = paste("d", distribution, sep = "")))
  n <- length(ord.x)
  P <- ppoints(n)
  z <- q.function(P, ...)
  plot(z, ord.x, type = "n", xlab = xlab, ylab = ylab, main = main, las = las,cex.lab = cex.lab, cex.axis = cex.axis)
  if (grid) {
    grid(lty = 1, equilogs = FALSE)
    box()}
  points(z, ord.x, col = col, pch = pch, cex = cex)
  if (line  == "quartiles" || line == "none") {
    Q.x <- quantile(ord.x, c(.25,.75))
    Q.z <- q.function(c(.25,.75), ...)
    b <- (Q.x[2] - Q.x[1])/(Q.z[2] - Q.z[1])
    a <- Q.x[1] - b*Q.z[1]
    abline(a, b, col = col.lines, lwd = lwd)
  }
  if (line == "robust") {
    coef <- coef(rlm(ord.x ~ z))
    a <- coef[1]
    b <- coef[2]
    abline(a, b)
  }
  conf <- if (envelope == FALSE) .95 else envelope
  zz <- qnorm(1 - (1 - conf)/2)
  SE <- (b/d.function(z, ...))*sqrt(P*(1 - P)/n)
  fit.value <- a + b*z
  upper <- fit.value + zz*SE
  lower <- fit.value - zz*SE
  if (envelope != FALSE) {
    lines(z, upper, lty = 2, lwd = lwd, col = col.lines)
    lines(z, lower, lty = 2, lwd = lwd, col = col.lines)
  }
}
##############################################################
############################ Plots ###########################
##############################################################
gVals <- function(y, dist, conf){ # distribution; confidence interval
  y <- sort(y)
  p <- ppoints(length(y))
  if(dist == "chisq") {
    zi <- qchisq(p, df = length(p) - 1)
    zd <- dchisq(zi, df = length(p) - 1)
    qz <- qchisq(c(.25, .75), length(p) - 1)
  } else {
    zi <- qnorm(p)
    zd <- dnorm(zi)
    qz <- qnorm(c(.25, .75))
  }
  coef <- coef(rlm(y~zi))
  a <- coef[1]
  b <- coef[2]
  z <- qnorm(1 - (1 - conf)/2)   # z = 1.96 for 95%...
  se <- (b / zd) * sqrt(p * (1 - p)/length(p))
  ft <- a + b * zi
  uc <- ft + z * se
  dc <- ft - z * se
  data.frame(z = zi, y = y, uc = uc, dc = dc)
}



