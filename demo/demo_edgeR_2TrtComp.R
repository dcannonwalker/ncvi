library(edgeR)
edgeRUsersGuide()

# simualted data for CRD with 2 trt groups and 3 replicates per group
# For how this dataset was simulated, check Simulation 3 of Kvam , Liu, and Si (2012).
data <- read.csv("Simu3_NB_rep3.csv", header = TRUE)
data <- data[apply(data[, 1:6], 1, sum) > 0, ]
head(data)

trt <- as.factor(c(1, 1, 1, 2, 2, 2))

dge <- DGEList(counts = data[, 1:6], group = trt)
dge <- calcNormFactors(dge)
names(dge)
dge$samples


# estimatation of phi_g in two-group comparison
# need to estimate common dispersion before estimating the tagwise dispersion
dge <- estimateCommonDisp(dge, verbose = TRUE)
sqrt(0.03028)

dge <- estimateTagwiseDisp(dge) 
# dge <- estimateDisp(dge)
hist(dge$tagwise.dispersion, breaks = 1000)

# detecting DE genes
et <- exactTest(dge)
topTags(et)
names(et)

# p-values adjusted for multiple testing using Benjamini and Hochberg method
p_bh <- p.adjust(et$table[, 3], method = "BH")

# number of rejected null hypotheses  
sum(p_bh <= 0.05)

## Because this is a simulated data, we know which 
## genes are truly DE and which genes are not.

# number of true positives
sum(p_bh[data[, 7] == 1] <= 0.05)

# We can also check the performance 
# of the test and FDR control. 
# True positive rate
sum(p_bh[data[, 7] == 1] <= 0.05) / sum(data[, 7] == 1)

# False positive rate
sum(p_bh[data[, 7] == 0] <= 0.05) / sum(data[, 7] == 0)

# V/R
sum(p_bh[data[, 7] == 0] <= 0.05) / sum(p_bh <= 0.05)



##########################
# Androgen data example 
# two treatment comparison
##########################

targets <- readTargets("pnas_target.txt", row.names = "Name")
targets
# read in counts 
# the pnas_expression.txt dataset can be downloaded at
# https://sites.google.com/site/davismcc/useful-documents
x <- read.delim("pnas_expression.txt", row.names = 1, stringsAsFactors = FALSE)
head(x)

# Put the counts and other information into a DGEList object:
y <- DGEList(counts = x[, 1:7], group = targets$Treatment, 
             genes = data.frame(Length = x[, 8]))
colnames(y) <- targets$Label

head(y$counts)

#  sum(apply(y, 1, sum) == 0)
#  [1] 15558

# apply(y, 2, sum)

# keep genes that has at least one count per million (cpm) in at least three samples:
keep <- rowSums(cpm(y) > 1) >= 3
y <- y[keep, ]
dim(y)

# Re-compute the library sizes:
y$samples$lib.size <- colSums(y$counts)

# Compute effective library sizes using TMM normalization:
y <- calcNormFactors(y)
y$samples


# An MDS plots shows distances, in terms of leading log-FC between pair of samples:
plotMDS(y)
x11()
# An MDS plots shows distances, in terms of BCV between pair of samples:
plotMDS(y, method = "bcv")

# The common dispersion estimates the overall BCV of the dataset, averaged over all genes:
y <- estimateCommonDisp(y, verbose = TRUE)
# estimate gene-specific dispersions:
y <- estimateTagwiseDisp(y)
# Plot the estimated dispersions:
plotBCV(y)

# Perform the exact tests for DE between treatment and control
et <- exactTest(y)
top <- topTags(et)
top

# The total number of DE genes at 5% FDR is given by
summary(de <- decideTestsDGE(et))

# Plot the log-fold-changes, highlighting the DE genes:
detags <- rownames(y)[as.logical(de)]
plotSmear(et, de.tags = detags)
abline(h = c(-1, 1), col = "blue")



##################################################
# edgeR GLM analysis for DE genes
# suppose there are 3 treatment groups, each with 2 replicates
##################################################

group <- factor(c(1, 1, 2, 2, 3, 3))

# 1. Treatment-contrast parameterization
design <- model.matrix(~ group)
design

# if you have estimated dispersion and normalization factor, 
# you could use the following line to fit the GLM model
# suppose that y is a DGEList that contains the counts 
# and normalization factors 
fit <- glmFit(y, design)

# To compare group 2 vs 1:
lrt.2vs1 <- glmLRT(fit, coef = 2)
topTags(lrt.2vs1)
# To compare group 3 vs 1:
lrt.3vs1 <- glmLRT(fit, coef = 3)
topTags(lrt.3vs1)
# To compare 3 vs 2:
lrt.3vs2 <- glmLRT(fit, contrast = c(0, -1, 1))
topTags(lrt.3vs2)
# To find genes different between any of the three groups:
lrt <- glmLRT(fit, coef = 2:3)
topTags(lrt)

# 2. Group-means parameterization
design <- model.matrix(~ 0 + group)
design
fit <- glmFit(y, design)

# To compare group 2 vs 1:
lrt.2vs1 <- glmLRT(fit, contrast = c(-1, 1, 0))
topTags(lrt.2vs1)
# To compare group 3 vs 1:
lrt.3vs1 <- glmLRT(fit, contrast = c(-1, 0, 1))
topTags(lrt.3vs1)
# To compare 3 vs 2:
lrt.3vs2 <- glmLRT(fit, contrast = c(0, -1, 1))
topTags(lrt.3vs2)
# To find genes different between any of the three groups:
my.contrasts <- makeContrasts(group2vs1 = group2 - group1, 
                              group3vs1 = group3 - group2, 
                              levels = design)
lrt <- glmLRT(fit, contrast = my.contrasts)
topTags(lrt)

