## simulate data for comparison to edgeR HT method
library(edgeR)
library(tidyr)
devtools::load_all()
set.seed(193)
## reject if abs(beta1) > cut_off
cut_off <- 5.85

## back of the napkin calculations
qnorm(0.1, mean = 0, sd = 1, lower.tail = F)
pnorm(1, mean = 1.85, lower.tail = T)

## set up
P = 2
G = 100
U = 6
N = 12
mu0 <- c(0, cut_off - 0.85)
precision_beta = 1
precision_u = 5
beta <- mvtnorm::rmvnorm(G, mean = mu0,
                         sigma = diag(1/precision_beta, P))
mu0
coin <- as.logical(rbinom(G, 1, 0.5))

for (i in seq(1,G)) {
  if (beta[i, 2] < cut_off) beta[i, 2] <- 0
  # if (coin[i]) beta[i, 2] = -beta[i, 2]
}

beta

u <- mvtnorm::rmvnorm(G, mean = rep(0, U),
                      sigma = diag(1/precision_u, U))
data("data_REH")
# X <- data_REH$C[,1:P]
# Z <- data_REH$C[,(P+1):(P+U)]
C <- data_REH$C
phi <- cbind(beta, u)
eta <- c(apply(phi, 1, function(p) C %*% p))
y_vector = rpois(N * G, lambda = exp(eta))
y <- vector2list(y_vector, G, N)
y

data_for_fit_ncvi <- list(y = y,
                 C = C,
                 G = G,
                 P = P,
                 U = U,
                 y_vector = y_vector)
init_for_fit_ncvi <- generate_inits(G, N, P, U,
                                    precision_mu0 = 1,
                                    precision_beta = precision_beta,
                                    precision_u = precision_u,
                                    mu_adjust = 3)
fit_vi <- fit_ncvi(data = data_for_fit_ncvi,
         init = init_for_fit_ncvi,
         elbo = elbo_REH,
         differentials = differentials_new,
         update_pars = update_pars_REH)

t(sapply(fit_vi$pars$phi, function(x) x$mu[1:2]))
beta

dfcounts <- data.frame(t(data.frame(y)))
rownames(dfcounts) <- paste0("symbol", 1:G)
colnames(dfcounts) <- c(paste0("cont", 1:U),
                        paste0("trtm", 1:U))
group <- factor(c(rep(1, N / P), rep(2, N / P)))


# data("arab")
# head(arab)
Treat <- factor(substring(colnames(dfcounts), 1, 4))
Treat <- relevel(Treat, ref = "cont")
Subj <- factor(substring(colnames(dfcounts), 5, 5))

data_edgeR <- DGEList(counts = dfcounts, group = Treat)
design = model.matrix(~ Subj + Treat)
rownames(design) <- colnames(data_edgeR)

data_edgeR <- estimateDisp(data_edgeR,
                           design = design,
                           robust = T)
data_edgeR$common.dispersion
plotBCV(data_edgeR)
fit <- glmQLFit(data_edgeR, design = design, robust = T)
plotQLDisp(fit)
qlf <- glmQLFTest(fit)
top <- rownames(topTags(qlf))
cpm(data_edgeR)[top,]

summary(decideTests(qlf))
test <- decideTests(qlf)


ci <- make_ci(fit_vi, par_names = c("beta1"))
ci

ci[beta[, 2] > 0, ]
View(ci[beta[, 2] == 0, ])

test <- make_test(fit_vi, cut_off, par_names = "beta1")
test
test %>% dplyr::arrange(lowertail)

