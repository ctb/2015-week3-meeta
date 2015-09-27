source("http://bioconductor.org/biocLite.R")
biocLite('edgeR')
biocLite('DESeq2')

library(edgeR)
library(DESeq2)
library(limma)

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

# Create DESeq2 datasets
dds <- DESeqDataSetFromMatrix(countData = exprs(bottomly.eset), colData = pData(bottomly.eset), design = ~ strain )
dds <- DESeq(dds)

dds.5rep <- DESeqDataSetFromMatrix(countData = exprs(bottomly.5reps), colData = pData(bottomly.5reps), design = ~ strain )
dds.5rep <- DESeq(dds.5rep)

dds.2rep <- DESeqDataSetFromMatrix(countData = exprs(bottomly.2reps), colData = pData(bottomly.2reps), design = ~ strain )
dds.2rep <- DESeq(dds.2rep)


# Plot dispersion estimates
plotDispEsts(dds)
plotDispEsts(dds.5rep)
plotDispEsts(dds.2rep)


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
plotBCV(dge)
plotBCV(dge.5reps)
plotBCV(dge.2reps)


# Create design matrix
design <- model.matrix(~ pData(bottomly.eset)$strain)

# Apply voom transformation
nf <- calcNormFactors(bottomly.eset)
v <- voom(exprs(bottomly.eset), design, lib.size=colSums(exprs(bottomly.eset))*nf, normalize.method="quantile", plot=TRUE)

# Do same for 5 replicate dataset
design <- model.matrix(~ pData(bottomly.5reps)$strain)
nf <- calcNormFactors(bottomly.5reps)
v.5reps <- voom(exprs(bottomly.5reps), design, lib.size=colSums(exprs(bottomly.5reps))*nf, 
                normalize.method="quantile", plot=TRUE)

# Do same for 2 replicates dataset
design <- model.matrix(~ pData(bottomly.2reps)$strain)
nf <- calcNormFactors(bottomly.2reps)
v.2reps <- voom(exprs(bottomly.2reps), design, lib.size=colSums(exprs(bottomly.2reps))*nf, 
                normalize.method="quantile", plot=TRUE)

#####

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
# Create design matrix
design <- model.matrix(~ pData(bottomly.eset)$strain)

# Usual limma pipeline
fit.voom <- lmFit(v, design)
fit.voom <- eBayes(fit.voom)

voom_results <- topTable(fit.voom, coef=2,  adjust="BH", number = nrow(exprs(bottomly.eset)))
voom_results$threshold <- as.logical(voom_results$adj.P.Val < p.threshold)
genes.voom <- row.names(voom_results)[which(voom_results$threshold)]

install.packages('gplots')

library(gplots)
venn(list(edgeR = genes.edgeR, DESeq2 = genes.deseq, voom = genes.voom))

