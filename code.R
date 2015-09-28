install.packages('gplots')
biocLite('edgeR')
biocLite('DESeq2')

source("http://bioconductor.org/biocLite.R")
library(edgeR)
library(DESeq2)
library(limma)
library(gplots)
library(Biobase)

# ReCount data: http://bowtie-bio.sourceforge.net/recount/
# ReCount paper: http://www.biomedcentral.com/1471-2105/12/449
# Bottomly et al. data set: http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0017820
# (mouse)

# Load data from eset from ReCount
load('bottomly_eset.RData')

#Filtering out non-expressors
isexpr <- rowSums(cpm(exprs(bottomly.eset))>1) >= 3
sum(isexpr)
bottomly.eset <- bottomly.eset[isexpr,]

# Convert metadata to factors
pData(bottomly.eset)$experiment.number <- factor(pData(bottomly.eset)$experiment.number)
pData(bottomly.eset)$lane.number <- factor(pData(bottomly.eset)$lane.number)

# Get indices for each sample group 
C57 <- which(pData(bottomly.eset)$strain == 'C57BL/6J')
DBA <- which(pData(bottomly.eset)$strain == 'DBA/2J')

# Randomly sample 5 reps from each
reps <- c(sample(C57, 5), sample(DBA, 5))
bottomly.5reps <- bottomly.eset[ ,reps]

# Randomly sample 2 reps from each
reps <- c(sample(C57, 2), sample(DBA, 2))
bottomly.2reps <- bottomly.eset[ ,reps]

# Clean up
rm(reps)
rm(C57)
rm(DBA)

# Plot 2 replicate dataset
eset <- bottomly.2reps
cpm.mat <- log(cpm(exprs(eset)))
mean.vec <- apply(cpm.mat, 1, mean)
sdvec <- apply(cpm.mat, 1, sd)
plot(mean.vec, sdvec, pch=".", main="2 replicates", ylab="sd", xlab="Average logCPM")

# Plot 5 replicate dataset
eset <- bottomly.5reps
cpm.mat <- log(cpm(exprs(eset)))
mean.vec <- apply(cpm.mat, 1, mean)
sdvec <- apply(cpm.mat, 1, sd)
plot(mean.vec, sdvec, pch=".", main="5 replicates", ylab="sd", xlab="Average logCPM")

# Plot 10 replicate dataset
eset <- bottomly.eset
cpm.mat <- log(cpm(exprs(eset)))
mean.vec <- apply(cpm.mat, 1, mean)
sdvec <- apply(cpm.mat, 1, sd)
plot(mean.vec, sdvec, pch=".", main="10 replicates", ylab="sd", xlab="Average logCPM")

#
# DEseq2
#

# Create DESeq2 datasets
dds <- DESeqDataSetFromMatrix(countData = exprs(bottomly.eset), colData = pData(bottomly.eset), design = ~ strain )
dds <- DESeq(dds)

dds.5rep <- DESeqDataSetFromMatrix(countData = exprs(bottomly.5reps), colData = pData(bottomly.5reps), design = ~ strain )
dds.5rep <- DESeq(dds.5rep)

dds.2rep <- DESeqDataSetFromMatrix(countData = exprs(bottomly.2reps), colData = pData(bottomly.2reps), design = ~ strain )
dds.2rep <- DESeq(dds.2rep)


# Plot dispersion estimates
plotDispEsts(dds.2rep)
plotDispEsts(dds.5rep)
plotDispEsts(dds)

#
# edgeR
#

dge <- DGEList(counts=exprs(bottomly.eset), group=pData(bottomly.eset)$strain)
# Normalize by total count
dge <- calcNormFactors(dge)

# Create the contrast matrix
design.mat <- model.matrix(~ 0 + dge$samples$group)
colnames(design.mat) <- levels(dge$samples$group)

# Estimate dispersion parameter for GLM
dge <- estimateGLMCommonDisp(dge, design.mat)
dge <- estimateGLMTrendedDisp(dge, design.mat, method="power")
dge<- estimateGLMTagwiseDisp(dge,design.mat)

# Do it all over again for 5 replicates
dge.5reps <- DGEList(counts=exprs(bottomly.5reps), group=pData(bottomly.5reps)$strain)
dge.5reps <- calcNormFactors(dge.5reps)
design.mat <- model.matrix(~ 0 + dge.5reps$samples$group)
colnames(design.mat) <- levels(dge.5reps$samples$group)

dge.5reps <- estimateGLMCommonDisp(dge.5reps, design.mat)
dge.5reps <- estimateGLMTrendedDisp(dge.5reps, design.mat, method="power")
dge.5reps<- estimateGLMTagwiseDisp(dge.5reps,design.mat)

# Do it all over again for 2 replicates
dge.2reps <- DGEList(counts=exprs(bottomly.2reps), group=pData(bottomly.2reps)$strain)
dge.2reps <- calcNormFactors(dge.2reps)
design.mat <- model.matrix(~ 0 + dge.2reps$samples$group)
colnames(design.mat) <- levels(dge.2reps$samples$group)

dge.2reps <- estimateGLMCommonDisp(dge.2reps, design.mat)
dge.2reps <- estimateGLMTrendedDisp(dge.2reps, design.mat, method="power")
dge.2reps<- estimateGLMTagwiseDisp(dge.2reps,design.mat)

# Plot mean-variance
plotBCV(dge.2reps)
plotBCV(dge.5reps)
plotBCV(dge)

#
# limma/voom
#

# Create design matrix for 2 replicates data set
design <- model.matrix(~ pData(bottomly.2reps)$strain)

# apply voom transformation & plot
nf <- calcNormFactors(bottomly.2reps)
v.2reps <- voom(exprs(bottomly.2reps), design, lib.size=colSums(exprs(bottomly.2reps))*nf, 
                normalize.method="quantile", plot=TRUE)

# Do same for 5 replicate dataset
design <- model.matrix(~ pData(bottomly.5reps)$strain)
nf <- calcNormFactors(bottomly.5reps)
v.5reps <- voom(exprs(bottomly.5reps), design, lib.size=colSums(exprs(bottomly.5reps))*nf, 
                normalize.method="quantile", plot=TRUE)

# Same, for all data.
design <- model.matrix(~ pData(bottomly.eset)$strain)
nf <- calcNormFactors(bottomly.eset)
v <- voom(exprs(bottomly.eset), design, lib.size=colSums(exprs(bottomly.eset))*nf, normalize.method="quantile", plot=TRUE)

##### Now, let's calculate differentially expressed genes for 2 repl at p 0.05

p.threshold <- 0.05

## edgeR ##
# Design matrix
design.mat <- model.matrix(~ 0 + dge.2reps$samples$group)
colnames(design.mat) <- c("C57BL", "DBA")

# Model fitting
fit.edgeR <- glmFit(dge.2reps, design.mat)

# Differential expression
contrasts.edgeR <- makeContrasts(C57BL - DBA, levels=design.mat)
lrt.edgeR <- glmLRT(fit.edgeR, contrast=contrasts.edgeR)

# Access results tables
edgeR_results <- lrt.edgeR$table
sig.edgeR <- decideTestsDGE(lrt.edgeR, adjust.method="BH", p.value = p.threshold)
genes.edgeR <- row.names(edgeR_results)[which(sig.edgeR == 1)]

## DESeq2 ##
contrast.deseq2 <- list("strainC57BL.6J", "strainDBA.2J")
deseq2_results <- results(dds, contrast=contrast.deseq2)
deseq2_results$threshold <- as.logical(deseq2_results$padj < p.threshold)
genes.deseq <- row.names(deseq2_results)[which(deseq2_results$threshold)]

## voom-limma ##
design <- model.matrix(~ pData(bottomly.2reps)$strain)
fit.voom <- lmFit(v.2reps, design)
fit.voom <- eBayes(fit.voom)

voom_results <- topTable(fit.voom, coef=2,  adjust="BH", number = nrow(exprs(bottomly.eset)))
voom_results$threshold <- as.logical(voom_results$adj.P.Val < p.threshold)
genes.voom <- row.names(voom_results)[which(voom_results$threshold)]

# now, look at overlap

venn(list(edgeR = genes.edgeR, DESeq2 = genes.deseq, voom = genes.voom))

##### Now, let's calculate differentially expressed genes for 5 repl at p 0.05

p.threshold <- 0.05

## edgeR ##
# Design matrix
design.mat <- model.matrix(~ 0 + dge.5reps$samples$group)
colnames(design.mat) <- c("C57BL", "DBA")

# Model fitting
fit.edgeR <- glmFit(dge.5reps, design.mat)

# Differential expression
contrasts.edgeR <- makeContrasts(C57BL - DBA, levels=design.mat)
lrt.edgeR <- glmLRT(fit.edgeR, contrast=contrasts.edgeR)

# Access results tables
edgeR_results <- lrt.edgeR$table
sig.edgeR <- decideTestsDGE(lrt.edgeR, adjust.method="BH", p.value = p.threshold)
genes.edgeR <- row.names(edgeR_results)[which(sig.edgeR == 1)]

## DESeq2 ##
contrast.deseq2 <- list("strainC57BL.6J", "strainDBA.2J")
deseq2_results <- results(dds.5rep, contrast=contrast.deseq2)
deseq2_results$threshold <- as.logical(deseq2_results$padj < p.threshold)
genes.deseq <- row.names(deseq2_results)[which(deseq2_results$threshold)]

## voom-limma ##
design <- model.matrix(~ pData(bottomly.5reps)$strain)
fit.voom <- lmFit(v.5reps, design)
fit.voom <- eBayes(fit.voom)

voom_results <- topTable(fit.voom, coef=2,  adjust="BH", number = nrow(exprs(bottomly.eset)))
voom_results$threshold <- as.logical(voom_results$adj.P.Val < p.threshold)
genes.voom <- row.names(voom_results)[which(voom_results$threshold)]

# now, look at overlap

venn(list(edgeR = genes.edgeR, DESeq2 = genes.deseq, voom = genes.voom))

##### Finally, let's calculate differentially expressed genes for all repl at p 0.05

p.threshold <- 0.05

## edgeR ##
# Design matrix
design.mat <- model.matrix(~ 0 + dge$samples$group)
colnames(design.mat) <- c("C57BL", "DBA")

# Model fitting
fit.edgeR <- glmFit(dge, design.mat)

# Differential expression
contrasts.edgeR <- makeContrasts(C57BL - DBA, levels=design.mat)
lrt.edgeR <- glmLRT(fit.edgeR, contrast=contrasts.edgeR)

# Access results tables
edgeR_results <- lrt.edgeR$table
sig.edgeR <- decideTestsDGE(lrt.edgeR, adjust.method="BH", p.value = p.threshold)
genes.edgeR <- row.names(edgeR_results)[which(sig.edgeR == 1)]

## DESeq2 ##
contrast.deseq2 <- list("strainC57BL.6J", "strainDBA.2J")
deseq2_results <- results(dds, contrast=contrast.deseq2)
deseq2_results$threshold <- as.logical(deseq2_results$padj < p.threshold)
genes.deseq <- row.names(deseq2_results)[which(deseq2_results$threshold)]

## voom-limma ##
design <- model.matrix(~ pData(bottomly.eset)$strain)
fit.voom <- lmFit(v, design)
fit.voom <- eBayes(fit.voom)

voom_results <- topTable(fit.voom, coef=2,  adjust="BH", number = nrow(exprs(bottomly.eset)))
voom_results$threshold <- as.logical(voom_results$adj.P.Val < p.threshold)
genes.voom <- row.names(voom_results)[which(voom_results$threshold)]

# now, look at overlap

venn(list(edgeR = genes.edgeR, DESeq2 = genes.deseq, voom = genes.voom))

