# Datamerge: A function I wrote to merge together two data frames where one or both have variables the other doesn't
# Note that this function assumes variables of the same name in different dataframes are the same variable

# Note that this is now (and possible was when I wrote it, though it has been a while) implemented in plyr
# by rbind.fill. I will maintain it here so that legacy code doesn't break, but will try to remember to use
# rbind.fill in future.
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

# Just a function to save a line of code when reading in spss files using haven. Though I never remember
# to use it.
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

# I believe that with the label = TRUE argument, this applies labels from Qualtrics to the data in exactly
# the same way as haven applies descriptions from spss, meaning that they are available through getdesc. 
# this seems to be the case in my experience, as it seemed to work in my dating profile examples. Thus,
# the description argument (and the description function) are likely deprecated.

# Updated this function to use the new readr version of read_csv. 
qual.load <- function(x, label = TRUE) {
  require(Hmisc)
  require(readr)
  data <- read_csv(x, skip = 2, col_names = FALSE)
  labs <- read_csv(x, n_max = 1)
  colnames(data) <- colnames(labs)
  labs <- as.character(labs[1,]); names(labs) <- 1:length(labs)
  if (label == TRUE) {
    label(data) <- lapply(labs, function(x) { label(data[, as.numeric(names(x))]) = x })
    
  }
  rm(labs)
  data
}



# Just a method to print the output from FA nicely without having to specify 20 settings, because it's a garbled
# mess of a function that spits out way too much output and is way too difficult to interpret.
nice_fa <- function(data, nfactors, rotate = 'varimax', loadings.only = TRUE, cut = .3, digits = 2, sort = TRUE, ...) {
  args <- list(...)
  if (loadings.only) {
    print(loadings(fa(data, nfactors = nfactors, rotate = rotate, args)), cut = cut, digits = digits, sort = sort)
  } else {
    print((fa(data, nfactors = nfactors, rotate = rotate, args)), cut = cut, digits = digits, sort = sort)
  }
  
}



# Leveraging the code that tells the alpha function in psych when things are likely reverse coded, to 
# automagically reverse code entire scales. This returns the entire dataset so needs to be used assignment
# wise. Update
auto.reverse.code <- function(data, max = NULL, min = NULL) {
  if (is.null(max)) {
    max <- max(data, na.rm = TRUE)
  }
  
  if (is.null(min)) {
    min <- min(data, na.rm = TRUE)
  }
  
  require(psych)
  p1 <- principal(data)
  if (any(p1$loadings < 0)) {
    # Converting to logical because tibbles dont like being subset by single dimension matrices. 
    reversed <- as.logical(p1$loadings < 0)
    # This uses the empirical max of the entire scale, whereas I think the reverse.code function in psych
    # uses only the item by default, though there you can specify your max and min
    data[, reversed] <- (max + min - data[, reversed])
    colnames(data)[reversed] <- paste0(colnames(data)[reversed], '-')
    print(paste('items ', paste(colnames(data)[reversed], collapse = ', '), ' were reverse coded.'))
  }
  data
}

# I very often have to create scales from a grep of scale names. So just making a shortcut to save 
# some time. 
create_scale <- function(data, scale) {
  data[[scale]] <- rowMeans(data[, grep(scale, colnames(data))], na.rm = TRUE)
  data
}

# A function to check and then load packages (though only on CRAN, not from github)

safe_load <- function(packages) {
  
  if (length(setdiff(packages, installed.packages()[,1])) > 0) {
    install.packages(setdiff(packages, installed.packages()))
  }
  
  sapply(packages, require, character.only = TRUE)
}




### Taking inspiration from the APA package, this is a simpler function that uses the broom package to
### print regression output easily in APA format. 

regression_apa <- function(model, variable) {
  require(broom)
  # This lets us enter variable names without strings, but doens't work with spaces (not that we 
  # should have any)
  variable <- deparse(substitute(variable))
  tab <- tidy(model, conf.int = TRUE)
  # Changed this so it finds the variable or the variable enclosed by brackets so we can capture 'scale' 
  # use. Have made as restrictive as possible so shouldn't ever find two variables, except in weird edge
  # cases with nonstandard bracket use in variable names. 
  tab <- tab[grepl(paste0('^', variable, '$'), tab$term)|grepl(paste0('\\(', variable, '\\)'), tab$term),]
  beta <- tab['estimate']
  p <- tab['p.value']
  t <- tab['statistic']
  df <- glance(model)['df.residual']
  
  p <- round(p, 3)
  p <- max(p, .001)
  
  if (p == .001) {
    text = paste0('$\\beta$ = ', gsub('0\\.', '\\.', sprintf("%.2f", round(beta, 2))), 
                  ', *t*(', df, ') = ', gsub('0\\.', '\\.', (sprintf("%.2f", round(t, 2)))), 
                                             ', *p* < .001')
  } else {
    text = paste0('$\\beta$ = ', gsub('0\\.', '\\.', sprintf("%.2f", round(beta, 2))), 
                  ', *t*(', df, ') = ', gsub('0\\.', '\\.',(sprintf("%.2f", round(t, 2)))), 
                                             ', *p* = ', gsub('0\\.', '\\.', sprintf("%.3f", p)))
  }
  
  text
}

alpha_apa <- function(dataframe) {
  require(psych)
  text = paste0('$\\alpha$ = ', gsub('0\\.', '\\.',(sprintf("%.2f", psych::alpha(dataframe)$total$raw_alpha, 2))))
  text
}

fit_many_regressions <- function(dependent, predictors, data) {
  models <- vector(mode = 'list', length = length(dependent))
  names(models) <- dependent
  
  # So basically what I'm doing here is just fitting all the regressions that I'm going to need to report, using dependent as a list of DVs and predictors 
  # as a list of IVs. 
  for (i in 1:length(dependent)) {
    # Getting a bit clever with vectorized paste here, perhaps at risk of making this hard to read/change later. Just pasting scale around both
    # predictors before collapsing them. 
    
    # Get the predictor classes, and use a fairly hacky approach to set the function for factors to as.factor. 
    # Not that it's needed but it matches the closing bracket and makes things a bit easier for me. 
    
    classes <- sapply(predictors, function(x) class(data[[x]]))
    fun_list <- rep('scale(', length = length(classes))
    fun_list[classes == 'factor'] <- 'as.factor('
    # However by trying to scale everything we make it impossible to have factors as controls. 
    models[[i]] <- as.formula(paste0('scale(', dependent[i], ') ~ ', paste(paste0(fun_list, predictors, ')'), collapse = ' + ')))
  }
  
  fitted <- sapply(models, function(x) do.call('lm', list(x, data = data)))
  fitted
}

# Add a report-sig function to make things a little easier to update in knitr. Basically, if one of my 
# tests (regression or correlation, but it should work on anything that spits out a p value) is significant,
# it adds nothing. But if it's not significant, it adds 'not '. Hopefully I can use this flexibly so that
# my sentences will update correctly, at least at the rough level. 
flip_sig <- function(text, word = 'not') {
  require(stringr)
  outstring_neg <- c('not ', 'did not predict')
  outstring_pos <- c('', 'predicted')
  names(outstring_neg) <- c('not', 'predicted')
  names(outstring_pos) <- c('not', 'predicted')
  pval <- as.numeric(str_match(text, 'p.*(\\.[0-9]*)')[,2])
  if (pval < .05) {
    out = outstring_pos[word]
  } else {
    out = outstring_neg[word]
  }
  out
}

# Just a trimmed down version of cor_apa from the apa package with sensible defaults to save space when calling
# it in-text
cor_apa_slim <- function(x, y) {
  require(apa)
  cor_apa(cor.test(x, y), print = FALSE, format = 'markdown')
}

# This function hijacks apa.cor.table from the apaTables package, but alters it to better suit my liking
# (no random spaces between rows, removes a lot of the extra output like table notes and such for simplicities
# sake for now). For best use, wrap the resulting table in a pander() call in a knitr chunk with results
# = 'asis', and it should output in a fairly reasonable format to word. 
apa_cor_table <- function(input, double.stars = FALSE) {
  require(apaTables); require(pander)
  panderOptions('table.alignment.rownames', 'left')
  mytab <- apa.cor.table(input, show.conf.interval = FALSE)
  mytab <- as.data.frame(mytab$table.body[seq(1, nrow(mytab$table.body), 2), 2:ncol(mytab$table.body)], row.names = mytab$table.body[seq(1, nrow(mytab$table.body), 2), 'Variable'])
  # Replace the double significance stars with single stars for convenience
  if (!double.stars) {
    mytab[] <- sapply(mytab, function(x) gsub('\\*\\*', '*', x))
  }
  
  mytab
}





# A little convenience function to figure out what a number was divided by. I use it when I encounter
# data that someone else has re-scaled and need to figure out the original scale. target = 0 means we're
# looking for whole numbers, .5 means half numbers just in case. Can use with sapply over a vector of
# names in a table to test a whole vector of divided numbers. 
undivide <- function(n, target = 0) {
  found = FALSE
  while (!found & i < 100) {
    # Deal with floating point problems. Round to five decimal places to stop tiny errors being a big deal.
    placeholder = round(n*i, digits = 5)
    found <- (((placeholder) %% 1) < (target + .01)) & (((placeholder) %% 1) > (target - .01))
    i <- i + 1
  }
  
  if (i < 100) {
    print(paste0('Lowest whole number ',n,' might have been divided from is ', i-1))
    print(paste0('The original number before division could be ', (i-1)*n))
  } else {
    print(paste0('Could not find a divisor for ', n))
  }
  
}