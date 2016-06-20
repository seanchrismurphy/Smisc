# Datamerge: A function I wrote to merge together two data frames where one or both have variables the other doesn't
# Note that this function assumes variables of the same name in different dataframes are the same variable
Datamerge <- function(data1, data2) {
  add1 <- setdiff(colnames(data2), colnames(data1))
  add2 <- setdiff(colnames(data1), colnames(data2))
  hold1 <- data.frame(matrix(data = NA, nrow = nrow(data1), ncol=length(add1)))
  hold2 <- data.frame(matrix(data = NA, nrow = nrow(data2), ncol=length(add2)))
  colnames(hold1) <- add1
  colnames(hold2) <- add2
  end1 <- cbind(data1, hold1)
  end2 <- cbind(data2, hold2)
  final <- rbind(end1, end2)
  final
}

# All credit to the original source of this, which I think is here: https://stat.ethz.ch/pipermail/r-help/2008-March/156583.html
# I have just adapted it slightly and load it here for my convenience.
corstars <- function(x){ 
  require(Hmisc) 
  x <- as.matrix(x) 
  R <- rcorr(x)$r 
  p <- rcorr(x)$P 
  
  ## define notions for significance levels; spacing is important.
  mystars <- ifelse(p < .001, "***", ifelse(p < .01, "** ", ifelse(p < .05, "* ", " ")))
  
  ## trunctuate the matrix that holds the correlations to two decimal
  R <- format(round(cbind(rep(-1.11, ncol(x)), R), 2))[,-1] 
  
  ## build a new matrix that includes the correlations with their apropriate stars 
  Rnew <- matrix(paste(R, mystars, sep=""), ncol=ncol(x)) 
  diag(Rnew) <- paste(diag(R), " ", sep="") 
  rownames(Rnew) <- colnames(x) 
  colnames(Rnew) <- paste(colnames(x), "", sep="") 
  
  ## remove upper triangle
  Rnew <- as.matrix(Rnew)
  Rnew[upper.tri(Rnew, diag = TRUE)] <- ""
  diag(Rnew) <- 1
  Rnew <- as.data.frame(Rnew) 
  
  ## remove last column and return the matrix (which is now a data frame)
  Rnew <- cbind(Rnew[1:length(Rnew)])
  colnames(Rnew) <- 1:ncol(Rnew)
  return(Rnew) 
}

# Just a function to save a line of code when reading in spss files using haven
spss.load <- function(x) {
  require(haven)
  data <- read_sav(x)
  data <- as.data.frame(data, stringsAsFactors = FALSE)
  data
}

# Function to quickly pull the descriptions from an spss object loaded through haven
getdesc <- function(x) {
  sapply(x, function(y) attributes(y)$label)
}

# Qualtrics doubles up on column names by having the question names and descriptions each with their own column, which
# causes everything to be loaded as a factor/character. This function loads the first line, then everything after the
# second line, and puts them together so everything is unaltered.

# It also adds the description attribute to the data to keep those Qualtrics questions in there, and can be accessed 
# with the description function below in the same way as colnames, though the attributes aren't per-column like in 
# haven loadings. 

### Just realised I had missed the header = FALSE argument here which meant I was missing the first row of every dataset I've used
### this function on. Time to go back over stuff!
qual.load <- function(x) {
  data <- read.csv(x, stringsAsFactors = FALSE, skip = 2, header = FALSE)
  names <- read.csv(x, stringsAsFactors = FALSE, nrow = 2)
  descriptions <- names[1,]
  colnames(data) <- colnames(names)
  attributes(data)$description <- as.character(descriptions[1,])
  rm(names, descriptions)
  data
}

description <- function(x) { 
  return(attributes(x)$description)
}
