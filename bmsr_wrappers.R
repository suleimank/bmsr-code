## gene.expression.matrices is a list of expression matrices, one for each dataset.
## drug.response.matrices is a list of drug response matrices (e.g., aucs or ic50s), one
## for each dataset.
## both lists are assumed named with the corresponding dataset
prepare.bmsr.data <- function(gene.expression.matrices, drug.response.vectors = NULL) {

  data.y <- NULL
  
  if(!is.null(drug.response.vectors)) {
    if(!(all(sort(names(gene.expression.matrices)) == sort(names(drug.response.vectors))))) {
      stop("Expected gene expression and drug response matrices to be named according to dataset\n")
    }
  }

  ## Merge the drugs across data sets, ensuring that we only keep those present in all of them.
  if(!is.null(drug.response.vectors)) {
    common.drugs <- Reduce("intersect", lapply(drug.response.vectors, rownames))
    data.y <- llply(drug.response.vectors, .fun = function(df) df[common.drugs, , drop=F])
  }
  
  ## Restrict to common features / genes
  common.x <- Reduce("intersect", lapply(gene.expression.matrices, rownames))
  data.x <- llply(gene.expression.matrices, .fun = function(df) df[common.x, , drop=F])

  for(ds in names(gene.expression.matrices)) {
    flag <- !is.na(colSums(data.x[[ds]]))
    data.x[[ds]] <- data.x[[ds]][, flag]
  }

  if(!is.null(drug.response.vectors)) {  
    for(ds in names(drug.response.vectors)) {
      flag <- !is.na(data.y[[ds]])
      data.y[[ds]] <- data.y[[ds]][, flag, drop = F]
    }
  }

  ## Restrict to common samples across features and responses within each data set
  if(!is.null(drug.response.vectors)) {  
    for(ds in names(gene.expression.matrices)) {
      common.samples <- intersect(colnames(data.x[[ds]]), colnames(data.y[[ds]]))
      data.x[[ds]] <- data.x[[ds]][, common.samples, drop=F]
      data.y[[ds]] <- data.y[[ds]][, common.samples, drop=F]
    }
  }

  for(ds in names(gene.expression.matrices)) {
    colnames(data.x[[ds]]) <- paste0(colnames(data.x[[ds]]), "_", ds)
  }

  if(!is.null(drug.response.vectors)) {
    for(ds in names(drug.response.vectors)) {
      colnames(data.y[[ds]]) <- paste0(colnames(data.y[[ds]]), "_", ds)
    }
  }
      
  ## Label the samples according to dataset
  sample.labels <- unlist(llply(1:length(data.x), .fun = function(i) rep(names(data.x)[i], ncol(data.x[[i]]))))

  num.samples <- unlist(llply(data.x, ncol))
  names(num.samples) <- names(data.x)

  data.x <- as.matrix(Reduce("cbind", data.x))
  if(!is.null(drug.response.vectors)) {  
    data.y <- as.matrix(Reduce("cbind", data.y))
  }
  
  return(list("data.x" = data.x, "data.y" = data.y, "sample.labels" = sample.labels, "num.samples" = num.samples))
}

train.mssr_ <- function(data.x, data.y, num.samples) {

  ## Set up options
  model.file <- "bmsr.stan"

  ## Set model parameters
  p0 = 100 #model parameter p0
  opts = list(iter=2000,seeds=c(12,345,6789),inference="Sampling")

  ## Input data to MMSR is:
  ## S:  number of data sets
  ## nY: a named vector holding the number of samples in each data set (i.e., num rows in X and Y)
  ## dX: number of columns in X matrix
  ## dY: number of columns in Y matrix
  ## Y: matrix of responses (rows are samples; columns are drugs)
  ## X: matrix of features (rows are samples; columns are features)
  S <- length(num.samples)
  nY <- num.samples
  X <- t(data.x)
  Y <- t(data.y)
  dX <- ncol(X)
  dY <- ncol(Y)
  data <- list("S" = S, "nY" = nY, "dX" = dX, "dY" = dY, "X"= X, "Y" = Y, "p0" = p0)
  res = runSTAN(model.file, data, opts)
  res$features <- colnames(X)
  res
}

train.mssr <- function(expression.matrices, drug.response.vectors) {

  prep <- prepare.bmsr.data(expression.matrices, drug.response.vectors)
  data.x <- prep$data.x
  data.y <- prep$data.y
  sample.labels <- prep$sample.labels
  num.samples <- prep$num.samples
  
  y <- data.y

if(FALSE) {
  flag <- !is.na(y) & !is.na(colSums(data.x))
  y <- y[,flag, drop=F]
  data.x <- data.x[, flag]
  sample.labels <- sample.labels[flag]

  ## Calculate the number of samples in each dataset
  tmp <- as.data.frame(table(sample.labels))
  num.samples <- tmp$Freq
  names(num.samples) <- tmp[,1]
  num.samples <- num.samples[unique(sample.labels)]
  }
  
  cat(paste0("Applying MSSR with ", nrow(data.x), " features\n"))

  fit <- train.mssr_(data.x = data.x, data.y = y, num.samples = num.samples)
  fit

}

predict.mssr_ <- function(data.x, num.samples, out) {

  model.file = "bmsr.stan";
  predFunction = paste0("predict.",model.file);

  ## Input data to MMSR is:
  ## S:  number of data sets
  ## nY: a named vector holding the number of samples in each data set (i.e., num rows in X and Y)
  ## dX: number of columns in X matrix
  ## dY: number of columns in Y matrix
  ## Y: matrix of responses (rows are samples; columns are drugs)
  ## X: matrix of features (rows are samples; columns are features)
  S <- length(num.samples)
  nY <- num.samples
  X <- t(data.x)
  dX <- ncol(X)
  dataTest <- list("S" = S, "nY" = nY, "dX" = dX, "X"= X)
  yPred = predictSTAN(predFunction,out,xTest=dataTest$X,nTest=dataTest$nY,yN= NULL) # get prediction on test data
  vec <- as.numeric(yPred)
  names(vec) <- rownames(dataTest$X)
  vec
}

predict.mssr <- function(expression.matrices, out) {

  prep <- prepare.bmsr.data(expression.matrices, drug.response.vectors = NULL)
  data.x <- prep$data.x
  sample.labels <- prep$sample.labels
  num.samples <- prep$num.samples

  if(FALSE) {
  flag <- !is.na(colSums(data.x))
  data.x <- data.x[, flag]
  sample.labels <- sample.labels[flag]

  ## Calculate the number of samples in each dataset
  tmp <- as.data.frame(table(sample.labels))
  num.samples <- tmp$Freq
  names(num.samples) <- tmp[,1]
  num.samples <- num.samples[unique(sample.labels)]
  }
  
  predict.mssr_(data.x, num.samples, out)
}
