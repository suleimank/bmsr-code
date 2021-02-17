suppressPackageStartupMessages(library(pacman))
suppressPackageStartupMessages(p_load(Biobase))
suppressPackageStartupMessages(p_load(SummarizedExperiment))
suppressPackageStartupMessages(p_load(S4Vectors))
suppressPackageStartupMessages(p_load(PharmacoGx))
suppressPackageStartupMessages(p_load(rstan))
suppressPackageStartupMessages(p_load(plyr))
suppressPackageStartupMessages(p_load(dplyr))
suppressPackageStartupMessages(p_load(mRMRe))
suppressPackageStartupMessages(p_load(foreach))
suppressPackageStartupMessages(p_load(parallel))
suppressPackageStartupMessages(p_load(ggplot2))

# Integrative modeling of Lapatinib across GDSC and CCLE breast cancer cell lines
# using BMSR for joint analysis over GDSC and CCLE, mRMRe for initial feature selection, and
# PharmacoGx for GDSC and CCLE data access.

num.cores <- detectCores()
if(!is.na(num.cores) && (num.cores > 1)) {
  suppressPackageStartupMessages(p_load("doMC"))
  cat(paste("Registering ", num.cores-1, " cores.\n", sep=""))
  registerDoMC(cores=(num.cores-1))
}
num.processes <- num.cores - 1
## Following is for rstan
options(mc.cores = num.processes)

source("dataset_utils.R")
source("bmsr.R")
source("bmsr_wrappers.R")

model.file = "bmsr.stan";
mssr.model = rstan::stan_model(file=model.file)

# Download the GDSC and CCLE datasets (using PharmacoGx)
datasets <- download.and.prepare.datasets()

# Post-process the expression and drug response data to: (1) exclude all-NA samples, (2) ensure expression and drug response
# data have the same samples, and (3) standardize both expression and drug response

expression.matrices <- llply(datasets, .fun = function(dataset) dataset[["expr"]])
drug.response.matrices <- llply(datasets, .fun = function(dataset) dataset[["response"]])

# Exclude all-NA samples
expression.matrices <- llply(expression.matrices, .fun = function(mat) na.omit(mat[, !unlist(lapply(1:ncol(mat), function(i) all(is.na(mat[,i]))))]))
drug.response.matrices <- llply(drug.response.matrices, .fun = function(mat) mat[, !unlist(lapply(1:ncol(mat), function(i) all(is.na(mat[,i]))))])

# Ensure samples are consistent between expression and drug response data
for(ds in names(expression.matrices)) {
  common.samples <- intersect(colnames(expression.matrices[[ds]]), colnames(drug.response.matrices[[ds]]))
  expression.matrices[[ds]] <- expression.matrices[[ds]][, common.samples, drop=F]
  drug.response.matrices[[ds]] <- drug.response.matrices[[ds]][, common.samples, drop=F]
}

# Standardize expression and drug response data
expression.matrices <- llply(expression.matrices, .fun = function(mat) t(scale(t(mat), center = TRUE, scale = TRUE)))
drug.response.matrices <- llply(drug.response.matrices, .fun = function(mat) t(scale(t(mat), center = TRUE, scale = TRUE)))

# Extract response for drug to model
drug <- "Lapatinib"
drug.response.vectors <- llply(drug.response.matrices, .fun = function(mat) mat[drug, ,drop=FALSE])

# Perform feature/gene selection using mRMRe (ensemble-based minimum redundancy, maximum relevance)
# to improve computational efficiency of BMSR
nms <- names(datasets)
names(nms) <- nms

features <-
  llply(nms, .parallel = TRUE,
        .fun = function(ds) {
                 new.df <- data.frame(target = as.numeric(drug.response.vectors[[ds]]), t(expression.matrices[[ds]]))
                 new.df <- new.df[!is.na(new.df$target), ]
                 data <- mRMR.data(data = new.df)
                 res <- mRMR.ensemble(data = data, target_indices = 1, feature_count = 50, solution_count = 20)
		 genes <- lapply(1:ncol(res@filters[[1]]), function(i) res@feature_names[res@filters[[1]][,i]])
		 sort(unique(as.character(unlist(genes))))
               })

features <- unique(unname(unlist(features)))

# Limit expression matrices to genes selected by mRMRe
expression.matrices <- llply(expression.matrices, .fun = function(mat) mat[features, ])

# Set a random seed
set.seed(1234)

# Train BMSR
fit <- train.mssr(expression.matrices, drug.response.vectors)
out = fit$out

# Predict drug response (here using test data = training data)
ypred <- predict.mssr(expression.matrices, out)

# Plot predicted versus observed drug response
prep <- prepare.bmsr.data(expression.matrices, drug.response.vectors)
yobs <- as.numeric(prep$data.y)
names(yobs) <- colnames(prep$data.y)

common.samples <- intersect(names(ypred), names(yobs))
ypred <- as.numeric(ypred[common.samples])
yobs <- as.numeric(yobs[common.samples])

df <- data.frame(ypred = ypred, yobs = yobs)
g <- ggplot(data = df, aes(x = ypred, y = yobs)) + geom_point() + xlab("Predicted Response") + ylab("Observed Response")

# Get the posterior means for the beta vectors across both GDSC and CCLE
posteriorFunction = paste0("posterior.",model.file);
post = getPosterior(posteriorFunction,out)

betaShared = as.numeric(post$betaShared)
names(betaShared) <- features

# Plot the posterior betas
df <- data.frame(gene = features, dataset1 = post$beta[1,], dataset2 = post$beta[2,])
g <- ggplot(data = df, aes(x = dataset1, y = dataset2)) + geom_point()
g <- g + xlab("GDSC Posterior Beta") + ylab("CCLE Posterior Beta")


