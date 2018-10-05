pkgname <- "Wrench"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('Wrench')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("wrench")
### * wrench

flush(stderr()); flush(stdout())

### Name: wrench
### Title: Normalization for sparse, under-sampled count data.
### Aliases: wrench

### ** Examples

#Obtain counts matrix and some group information
require(metagenomeSeq)
data(mouseData)
cntsMatrix <- MRcounts(mouseData)
group <- pData(mouseData)$diet
#Running wrench with defaults
W <- wrench( cntsMatrix, condition=group  )
compositionalFactors <- W$ccf
normalizationFactors <- W$nf

#Introducing the above normalization factors for the most
# commonly used tools is shown below.

#If using edgeR, we must pass in the compositional factors
require(edgeR)
edgerobj <- DGEList( counts=cntsMatrix,
                     group = as.matrix(group),
                     norm.factors=compositionalFactors )

#If using DESeq/DESeq2
require(DESeq2)
deseq.obj <- DESeqDataSetFromMatrix(countData = cntsMatrix,
                                   DataFrame(group),
                                   ~ group )
DESeq2::sizeFactors(deseq.obj) <- normalizationFactors
#If using metagenomeSeq
normalizedObject <- mouseData
pData(normalizedObject@expSummary$expSummary)$normFactors <- normalizationFactors




### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
