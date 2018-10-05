## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

require(metagenomeSeq)
require(edgeR)
require(Wrench)
require(DESeq2)


## ---- warning=FALSE------------------------------------------------------
#extract count and group information for from the mouse microbiome data in the metagenomeSeq package
data(mouseData)
counts <- MRcounts( mouseData, norm=F ) 
group <- pData(mouseData)$diet

#Running wrench with defaults
W <- wrench( counts, condition=group  )
compositionalFactors <- W$ccf
normalizationFactors <- W$nf


## ---- warning=FALSE------------------------------------------------------
#Introducing the above normalization factors for the most
# commonly used tools is shown below.

# -- If using metagenomeSeq
normalizedObject <- mouseData
normFactors(normalizedObject) <- normalizationFactors

# -- If using edgeR, we must pass in the compositional factors
edgerobj <- DGEList( counts=counts,
                     group = as.matrix(group),
                     norm.factors=compositionalFactors )

# -- If using DESeq/DESeq2
deseq.obj <- DESeqDataSetFromMatrix(countData = counts,
                                   DataFrame(group),
                                   ~ group )
sizeFactors(deseq.obj) <- normalizationFactors


## ---- warning=FALSE------------------------------------------------------
time <- as.numeric(as.character(pData(mouseData)$relativeTime))
time.levs <- cut( time, breaks = c(0, 6, 28, 42, 56, 70) )
overall_group <- paste( group, time.levs ) #merge the time information and the group information together
W <- wrench( counts, condition = overall_group )

