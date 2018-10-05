#' Obtains logistic fits for presence/absence and fitted probabilities of a zero occurring.
#'
#' This function is used to derive weights for feature-wise compositional estimates. Our (default)
#' intention is to derive these based on average occurrences across the dataset, as just a function
#' of sample depth, and not with particular relevance to groups.
#' @param mat count matrix
#' @param hdesign design matrix for the logistic; the default is usually sufficient.
#' @param thresh True if numerically one/zero probability occurrences must be thresholded
#' @param thresh.val if thresh is true, the numerically one/zero probability occurrences is thresholded
#'        to this value
#' @return A list with components:
#'          \itemize{
#'          \item{pi0.fit -  list with feature-wise glm.fit objects}
#'          \item{pi0 - matrix with fitted probabilities}
#'          }
#'
getHurdle <- function(mat, hdesign=model.matrix( ~-1+log(colSums(mat)) ),
                      thresh=F, thresh.val=1e-8, ... ){
  require(matrixStats)

  pi0.fit <- apply(mat, 1, function(x) glm.fit( hdesign, c(1*(x==0)), family=binomial() ) )
  pi0 <- t(sapply( pi0.fit, function(x) x$fitted.values  ))

  #tau <- colSums(mat)
  #pi0.fit <- apply(mat, 1, function(x) glm.fit( hdesign, cbind( tau-x, x ), family=binomial() ) )
  #pi0 <- t(sapply( pi0.fit, function(x) exp(tau*log(x$fitted.values)) ) )

  if(thresh){
    pi0[pi0>1-thresh.val] <- 1-thresh.val
    pi0[pi0<thresh.val] <- thresh.val
  }
  list("pi0.fit"=pi0.fit, "pi0"=pi0)
}

#' Obtain variances of logged counts.
#'
#'@param mat count matrix; rows are features and columns are samples.
#'@param design model matrix for the count matrix
#'@param plot if the mean-variance trend function (the same as that of voom) needs to be plot.
#'@param ebs2 if regularization of variances needs to be performed.
#'@param smoothed TRUE if all the variance estimates must be based on the mean-variance trend function.
#'@return a vector with variance estimates for logged feature-wise counts.
#'
gets2 <- function(mat, design=model.matrix(mat[1,]~1), plot=F, ebs2=T, smoothed=F, ...){

  require(matrixStats)
  require(locfit)
  require(limma)

  p <- nrow(mat)
  nzrows <- rowSums(mat)>0
  mat <- mat[nzrows,]

  tjs <- colSums(mat)
  tjs <- tjs/exp(mean(log(tjs)))

  design <- cbind( design, log(tjs) )
  mat <- log(mat)
  mat[!is.finite(mat)] <- NA
  fit <- lmFit( mat, design  )
  s <- fit$sigma
  s[s==0] <- NA
  mu <- rowMeans( mat, na.rm=T )

  sqrt.s <- sqrt(s)
  l <- locfit(sqrt.s~mu, family="gamma")

  if(plot){
    plot( mu, sqrt.s, pch=16, cex=.25 )
    lines(l, col='red')
  }

  k <- predict(l, mu)^4
  if(smoothed){
    s2 <- k
  } else {
    s2 <- s^2
    s2[is.na(s2)] <- k[is.na(s2)]
    if(ebs2){
      s2 <- limma::squeezeVar(s2, df = max(c(1, fit$df.residual), na.rm=T) )$var.post
    }
  }
  s2.tmp <- c(matrix(NA, nrow=p, ncol=1))
  s2.tmp[nzrows] <- s2
  s2.tmp
}

#' Marginal weight computations for wrench estimators.
#' @param res result structure of \code{wrench}
#' @param z.adj TRUE if the result structure was generated with \code{wrench} with \code{z.adj} set to TRUE.
getMargWeights <- function( res, z.adj, ...  ){
  with(res$others, {
    s2theta <- design %*% s2thetag
    tmp <- exp( sweep( replicate( nrow(design), s2 ), 2, s2theta, "+" ))-1
    if(!z.adj){
      W<-1/( (1-pi0)*( pi0 + tmp ) )
      W[!is.finite(W)] <- NA
    } else {
      W<-(1-pi0)/(pi0+ tmp)
    }
    colnames(W) <- colnames(res$others$r)
    W
  })
}

#' Postive-conditional weight computations for wrench estimators.
#' @param res result structure of \code{wrench}
getCondWeights <- function( res ) {
  with(res$others, {
    radj[radj==0] <- NA
    (1-pi0)/(( pi0 + replicate( nrow(design), exp(s2) ) -1 )*(radj^2))
  })
}

#' Log Postive-conditional weight computations for wrench estimators.
#' @param res result structure of \code{wrench}
getCondLogWeights <- function( res ) {
  with(res$others, {
    radj[radj==0] <- NA
    (1-pi0)/(( pi0 + replicate( nrow(design), exp(s2) ) -1 ))
  })
}

#' @export
getWeightedMean <- function( mat, w=rep(1, nrow(mat)) ){
  require(matrixStats)
  require(limma)
  if( is.vector(w) ){
    w <- c(w)
    res <- colWeightedMeans( mat, w )
  } else {
    res <- sapply( seq(ncol(mat)), function(j){
      yj <- mat[,j]
      wj <- w[,j]
      yj <- yj[!is.na(w)]
      wj <- wj[!is.na(w)]
      weighted.mean( yj, wj, na.rm=T )
    } )
  }
  return(res)
}

#' Obtain robust means. .
#' @param res result structure of \code{wrench}
#' @param estim.type estimator type
estimSummary <- function( res, estim.type="s2.w.mean", ...  ){
  require(matrixStats)

  if(estim.type=="s2.w.mean"){ #weights based on s2
    with(res$others, colWeightedMeans( radj, 1/s2 ))
  } else if( estim.type=="gs.mean" ){ #this performs well for real-life data
    with(res$others,{
      y <- apply( radj, 2, function(x) exp((1/length(x))*sum( log( x[x>0] ) )) )
      y/exp(mean(log(y)))
    })
  } else if(estim.type=="mean"){
    with(res$others,{
      colMeans(radj)/exp(mean(log(colMeans(radj))))
    })
  } else if( estim.type=="median" ) {
    with(res$others,  colMedians(radj)/exp(mean(log(colMedians(radj))))
    )
  } else if( estim.type=="w.marg.median" ){
    W <- getMargWeights( res, ... )
    with(res$others, {
      sapply( seq(ncol(radj)), function(j){
        y <- radj[,j]
        w <- W[,j]
        y <- y[!is.na(w)]
        w <- w[!is.na(w)]
        weighted.median( y, w )
      } )
    })
  } else if (estim.type=="w.marg.mean"){
    W <- getMargWeights( res, ... )
    with(res$others, {
      sapply( seq(ncol(radj)), function(j){
        y <- radj[,j]
        w <- W[,j]
        y <- y[!is.na(w)]
        w <- w[!is.na(w)]
        weighted.mean( y, w )
      } )
    })
  }else if( estim.type=="w.cond.median" ){
    W <- getCondWeights( res )
    with(res$others, {
      sapply( seq(ncol(radj)), function(j){
        y <- radj[,j]
        w <- W[,j]
        y <- y[!is.na(w)]
        w <- w[!is.na(w)]
        weighted.median( y, w )
      } )
    })
  } else if( estim.type=="w.cond.mean"){
    W <- getCondWeights( res )
    with(res$others, {
      sapply( seq(ncol(radj)), function(j){
        y <- radj[,j]
        w <- W[,j]
        y <- y[!is.na(w)]
        w <- w[!is.na(w)]
        weighted.mean( y, w )
      } )
    })
  } else if( estim.type=="w.cond.log.median" ){
    W <- getCondLogWeights( res )
    with(res$others, {
      sapply( seq(ncol(radj)), function(j){
        y <- log(radj[,j])
        y[!is.finite(y)] <- NA
        w <- W[,j]
        y <- y[!is.na(w)]
        w <- w[!is.na(w)]
        exp(weighted.median( y, w, na.rm=T ))
      } )
    })
  } else if( estim.type=="w.cond.log.mean"){
    W <- getCondLogWeights( res )
    with(res$others, {
      sapply( seq(ncol(radj)), function(j){
        y <- log(radj[,j])
        y[!is.finite(y)] <- NA
        w <- W[,j]
        y <- y[!is.na(w)]
        w <- w[!is.na(w)]
        exp(weighted.mean( y, w, na.rm=T ))
      } )
    })
  } else if( estim.type=="hurdle.w.mean"  ){
    #hurdle weights applied to regularized ratios
    W <- 1/(1-res$others$pi0) #hurdle weights
    with(res$others, {
      sapply( seq(ncol(r)), function(j){
        y <- r[,j]
        w <- W[,j]
        weighted.mean( r[,j], W[,j])
      } )
    })
  }
}

#'
getReference <- function( mat, ref.est="sw.means", ... ){
  require(matrixStats)
  tots <- colSums(mat)
  if(ref.est=="logistic"){
    qref <- 1-plogis(
      apply( mat, 1, function(x){
        glm( cbind(tots-x, x)~1, family=binomial() )$coefficients
      }  )
    )
  } else if(ref.est == "sw.means"){ #sample-wise means
    qmat <- sweep(mat, 2, colSums(mat), "/")
    qref <- rowMeans(qmat)
  } else {
    stop("Unknown reference type.")
  }
  qref
}

#'
detrend.ccf <- function(ccf, tau, grp, plt.detrends.all=F){
  logccf <- log(ccf)
  logtau <- log(tau)
  grp <- as.factor(grp)
  df <- data.frame("logccf"=logccf, "grp"=grp, "logtau"=logtau)

  for( g in levels(df$grp) ){ #fit lm and remove trend
    g.indcs <- which(df$grp==g)
    g.sub <- df[g.indcs, ]
    g.fit <- lm( logccf ~ logtau, data=g.sub )
    pred.values <- predict( g.fit )
    if(plt.detrends.all){
      plot( g.sub$logccf ~ g.sub$logtau, col = "black" )
      points( pred.values ~ g.sub$logtau, col = "red", pch=19 )
    }

    max.tau.indx <- which.max( g.sub$logtau )
    g.sub$logccf <- g.sub$logccf + ( pred.values[max.tau.indx] - pred.values )
    df$logccf[g.indcs] <- g.sub$logccf
  }

  ccfs.detr <- exp(df$logccf)
  res <- list()
  res$ccf.detr.un <- ccfs.detr
  res$ccf.detr <- ccfs.detr / exp( mean(log(ccfs.detr)) )
  res
}

#' @title Normalization for sparse, under-sampled count data.
#'
#' @description Obtain normalization factors for sparse, under-sampled count data that often arise with
#' metagenomic count data.
#'
#' @param mat count matrix; rows are features and columns are samples
#' @param condition a vector with group information on the samples
#' @param etype weighting strategy with the following options:
#'        \itemize{
#'        \item{ hurdle.w.mean, the W1 estimator in manuscript.}
#'        \item{ w.marg.mean, the W2 estimator in manuscript. These are appropriately computed depending on
#'                             whether \code{z.adj}=T (see below)
#'        }
#'        \item{s2.w.mean, weight by inverse of feature-variances of logged count data. }
#'        }
#' @param ebcf TRUE if empirical bayes regularization of ratios needs to be performed. Default recommended.
#' @param z.adj TRUE if the feature-wise ratios need to be adjusted
#'              by hurdle probabilities (arises when taking marginal expectation). Default recommended.
#' @param phi.adj TRUE if estimates need to be adjusted for variance terms
#'                (arises when considering positive-part expectations). Default recommended.
#' @param detrend FALSE if any linear dependence between sample-depth and compositional factors needs to be removed.
#'                (setting this to TRUE reduces variation in compositional factors and can improve accuracy, but requires an extra assumption that no linear dependence between compositional factors and sample depth is present in samples).
#' @return a \code{list} with components:
#'         \itemize{
#'         \item{ \code{nf}, \emph{normalization factors} for samples passed.
#'                Samples with zero total counts are removed from output. }
#'         \item{ \code{ccf}, \emph{compositional correction factors}.
#'                Samples with zero total counts are removed from output.
#'           }
#'         \item{ \code{others},  a \code{list} with results from intermediate computations. }
#'         \itemize{
#'         \item{ \code{qref},  reference chosen. }
#'         \item{ \code{design},  design matrix used for computation of positive-part parameters. }
#'         \item{ \code{s2},   feature-wise variances of logged count data. }
#'         \item{ \code{r},  (regularized) ratios of feature-wise proportions. }
#'         \item{ \code{radj},  adjustments made to the regularized ratios based
#'                              on z.adj and phi.adj settings.}
#'         }
#'         }
#' @examples
#' #Obtain counts matrix and some group information
#' require(metagenomeSeq)
#' data(mouseData)
#' cntsMatrix <- MRcounts(mouseData)
#' group <- pData(mouseData)$diet

#' #Running wrench with defaults
#' W <- wrench( cntsMatrix, condition=group  )
#' compositionalFactors <- W$ccf
#' normalizationFactors <- W$nf
#'
#'#Introducing the above normalization factors for the most
#'# commonly used tools is shown below.
#'

#' #If using metagenomeSeq
#' normalizedObject <- mouseData
#' normFactors(normalizedObject) <- normalizationFactors
#'

#' #If using edgeR, we must pass in the compositional factors
#'require(edgeR)
#'edgerobj <- DGEList( counts=cntsMatrix,
#'                      group = as.matrix(group),
#'                      norm.factors=compositionalFactors )
#'

#' #If using DESeq/DESeq2
#'require(DESeq2)
#'deseq.obj <- DESeqDataSetFromMatrix(countData = cntsMatrix,
#'                                    DataFrame(group),
#'                                    ~ group )
#'DESeq2::sizeFactors(deseq.obj) <- normalizationFactors

#' @author M. Senthil Kumar
#' @export
wrench <- function( mat, condition, etype="w.marg.mean",
                    ebcf=T, z.adj=F, phi.adj=T, detrend=F, ... ){

  require(matrixStats)

  #trim
  mat <- mat[rowSums(mat)>0,]
  nzcols <- colSums(mat)>0
  mat <- mat[,nzcols]
  condition <- condition[nzcols]


  #feature-wise parameters: hurdle, variance, reference and raw ratios
  n <- ncol(mat)
  tots <- colSums(mat)
  compute.pi0 <- !((etype %in% c("mean", "median", "s2.w.mean")) & !z.adj )
  if(compute.pi0){
    pi0 <- getHurdle( mat, ... )$pi0
  }
  group <- as.character(condition)
  if(length(unique(group)) == 1){
    design <- model.matrix(mat[1,]~1)
  } else {
    design <- model.matrix(~-1+group)
  }
  s2 <- gets2(mat,design,...)

  #refernece
  qref <- getReference( mat, ... )

  #sample-wise ratios
  qmat <- sweep(mat, 2, colSums(mat), "/")
  r <- qmat/qref

  if(ebcf){
    #group-wise ratios
    Yg <- sapply( unique(group), function(g){
      if(sum(group==g)>1){
        rowSums(mat[,group==g])
      } else {
        mat[,group==g]
      }
    }  )
    qg <- sweep(Yg, 2, colSums(Yg), "/") #weighted estimator
    rg <- qg/qref
    lrg <- log(rg)
    lrg[!is.finite(lrg)] <- NA
    s2thetag <- colVars(lrg, na.rm=T)
    s2thetag_rep <- design %*% s2thetag

    thetag <- colMeans(rg)
    thetag_rep <- c(design %*% thetag)

    #regularized estimation of positive means.
    r <- sweep(r, 2, thetag_rep, "/")
    thetagj <- exp( sapply(seq(n), function(j){
      x <- log(r[,j])
      x[!is.finite(x)] <- NA
      weighted.mean( x, w=1/(s2+s2thetag_rep[j]), na.rm=T )
    }))

    p <- nrow(r)
    thetagi <- t( sapply(seq(p), function(i){
      exp(
        (s2thetag_rep/(s2[i]+s2thetag_rep))*( log(r[i,]) - log(thetagj) )
      )
    }) )

    r <- sweep( thetagi, 2, thetagj*thetag_rep, "*")
  }

  #adjustments for marginal, and truncated means.
  phi2 <- exp(s2)
  radj <- r
  if(z.adj){
    radj <- radj/(1-pi0)
  }
  if(phi.adj){
    radj <- sweep( radj, 1, sqrt(phi2), "/" )
  }


  #return result structure
  res <- list()
  res$others <- list()
  if(ebcf){
    res$others <- list( "rg" = rg,
                        "thetag" = thetag,
                        "thetagi" = thetagi,
                        "s2thetag" = s2thetag
    )
  }
  res$others <- c(
    res$others,
    list(
      "qref"=qref,
      "design" = design,
      "s2" = s2,
      "r" = r,
      "radj" = radj)
  )
  if(compute.pi0){
    res$others <- c(res$others, list("pi0"=pi0))
  }

  res$ccf <- estimSummary( res, estim.type = etype, z.adj=z.adj, ... )
  res$ccf <- with(res, ccf/exp(mean(log(ccf))))

  if( detrend ){
    res$others$ccf0 <- res$ccf
    detrended <- detrend.ccf( res$ccf, tots, condition )
    res$others$ccf.detr.un <- detrended$ccf.detr.un
    res$ccf <- detrended$ccf.detr
  }

  tjs <- colSums(mat)/exp(mean(log(colSums(mat))))
  res$nf <- res$ccf * tjs


  res
}
