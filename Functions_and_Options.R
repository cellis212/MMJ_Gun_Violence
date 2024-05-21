
# Set table and figure paths ----------------------------------------------
if(grepl(pattern = "celli", getwd())){
  tabPath = "C:\Users\celli\Dropbox\Apps\Overleaf\Medical Cannabis and Gun Violence\Tables"
  figPath = "C:\Users\celli\Dropbox\Apps\Overleaf\Medical Cannabis and Gun Violence\Figures"
  mapPath = "C:\Users\celli\Dropbox\Apps\Overleaf\Medical Cannabis and Gun Violence\Maps"
  datPath = "C:\Users\celli\Dropbox\School\Papers\MMJ Gun Violence\Data"
}


#### Options I Like ####
options(dplyr.width = Inf) # Makes dplyr show all columns in commands like head
options(stringsAsFactors = FALSE) # loads string variables as strings instead of factors
options(scipen = 999) # turns off scientific notation
options(error = function() {
  stop("An error occurred. Check the console for more information.")
}) # Stop all execution on an error


# Function to transform suitable columns to numeric
transform_numeric <- function(data_frame, exclude_cols = NULL) {
  # Get the columns to be checked (excluding the ones provided)
  columns_to_check <- setdiff(names(data_frame), exclude_cols)
  
  for(col_name in columns_to_check) {
    # Check if the column is not already numeric
    if(!is.numeric(data_frame[[col_name]])) {
      # Identify which values cannot be converted to numeric
      problematic_entries <- is.na(as.numeric(data_frame[[col_name]]))
      
      # Convert the column to numeric if not all values are problematic
      if(any(problematic_entries) && !all(problematic_entries)) {
        data_frame[[col_name]] <- as.numeric(data_frame[[col_name]])
      }
    }
  }
  return(data_frame)
}

# Function to report problematic values
report_problematic <- function(data_frame, exclude_cols = NULL) {
  problematic_values <- list()
  
  # Get the columns to be checked (excluding the ones provided)
  columns_to_check <- setdiff(names(data_frame), exclude_cols)
  
  for(col_name in columns_to_check) {
    # Check if the column is not already numeric
    if(!is.numeric(data_frame[[col_name]])) {
      # Identify which values cannot be converted to numeric
      problematic_entries <- is.na(as.numeric(data_frame[[col_name]]))
      if(any(problematic_entries) && !all(problematic_entries)) {
        problematic_values[[col_name]] <- unique(data_frame[[col_name]][problematic_entries])
      }
    }
  }
  
  # Convert the list to a data frame for reporting
  df_problematic <- if(length(problematic_values) > 0) {
    data.frame(column = names(problematic_values), 
               non_numeric_values = sapply(problematic_values, paste, collapse=", "), 
               stringsAsFactors = FALSE)
  } else {
    NULL
  }
  return(df_problematic)
}

# Usage:
# df <- data.frame(a = c("1", "2", "hello"), b = c("3", "apple", "4"), c = c("banana", "5", "6"))
# transformed_df <- transform_numeric(df, exclude_cols = c("b"))
# problematic_values_df <- report_problematic(df, exclude_cols = c("b"))


# Usage:
# df <- data.frame(a = c("1", "2", "hello"), b = c("3", "apple", "4"), c = c("banana", "5", "6"))
# result <- transform_and_report(df, exclude_cols = c("b"))
# transformed_df <- result$transformed_data
# problematic_values_df <- result$report


# Usage:
# df <- data.frame(a = c("1", "2", "hello"), b = c("3", "apple", "4"), c = c("banana", "5", "6"))
# result <- transform_and_report(df)
# transformed_df <- result$transformed_data
# problematic_values_df <- result$report

####### Rounding for table creation #######
# Mean
meanPr <- function(data, digits = 3, comma = TRUE)
{
  m <- mean(data, na.rm = TRUE)
  if(abs(m) < 1000 |
     comma == FALSE){
    return(formatC(m, digits = digits, format = "f"))
  } else {
    return(paste0(paste0("\\multicolumn{2}{", formatC(m, digits = 0, format = "f")), "}"))
  }
}

# Standard Deviation
sdPr <- function(data, digits = 3)
{
  m <- sd(data, na.rm = TRUE)
  if(m < 1000){
    return(formatC(m, digits = digits, format = "f"))
  } else {
    return(paste0(paste0("\\multicolumn{2}{", formatC(m, digits = 0, format = "f")), "}"))
  }
}

# Andrews and Armstrong (2017) unbiased IV for exactly identified models=
aaniv <- function(form1,
                  form2,
                  data,
                  cvar,
                  nboot = 100)
{
  # Required Packages
  require(lfe)
  
  ests <- matrix(0,
                 nrow = nboot,
                 ncol = 2)
  
  for(b in 1:nboot){
    # Block Bootstrapping
    if(b > 1){
      cInd <- sample(unique(cvar), length(unique(cvar)), replace = TRUE)
      tInd <- c()
      for(t in 1:length(cInd)){
        tInd <- c(tInd, which(cvar == cInd[t]))
      }
      dataB <- data[tInd,]
    } else {
      dataB <- data
    }
    
    mod <- felm(form1, dataB)
    ests[b,1] <- mod$coefficients[1]
    
    mod <- felm(form2, dataB)
    ests[b,2] <- mod$coefficients[1]
  }
  
  vc <- matrix(c(var(ests[2:nboot,1]), 
                 cov(ests[2:nboot,1], ests[2:nboot,2]),
                 cov(ests[2:nboot,1], ests[2:nboot,2]),
                 var(ests[2:nboot,2])),
               nrow = 2,
               ncol = 2)
  
  bu <- 1/(vc[2,2]^.5) * 
    (1 - pnorm(ests[1,2]/(vc[2,2]^.5))) / dnorm(ests[1,2]/(vc[2,2]^.5)) * 
    (ests[1,1] - vc[1,2]/vc[2,2] * ests[1,2]) +
    vc[1,2]/vc[2,2]
  
  return(list("bu" = bu,
              "ests" = ests,
              "vc" = vc))
}

# Bootstrap diff-in-coefficients test
cdifboot <- function(g, # Function that inputs data and returns a coefficient
                     data,
                     cutVar, # character variable that defines subset
                     cvar, # cluster id variable
                     nboot = 100,
                     seed = 42)
{
  set.seed(seed)
  ests <- matrix(0,
                 nrow = nboot,
                 ncol = 2) 
  cutCol <- which(colnames(data) == cutVar) 
  
  for(b in 1:nboot){
    if(b > 1){
      cInd <- sample(unique(cvar), length(unique(cvar)), replace = TRUE)
      tInd <- c()
      for(t in 1:length(cInd)){
        tInd <- c(tInd, which(cvar == cInd[t]))
      }
      dataB <- data[tInd,]
    } else {
      dataB <- data
    }
    
    ests[b,1] <- g(dataB[dataB[,cutCol] == 0,])
    ests[b,2] <- g(dataB[dataB[,cutCol] == 1,])
  }
  
  dif <- ests[1,1] - ests[1,2]
  varDif <- var(ests[2:nboot, 1]) + var(ests[2:nboot, 2]) - 2 * cov(ests[2:nboot, 1], ests[2:nboot, 2])
  
  return(2*pnorm(-abs(dif/(varDif)^.5)))
}

# T-test for difference in means
tdif <- function(var1, var2, digits = 3)
{
  # Come back and add stars
  m <- t.test(var1, var2)
  return(formatC(m$statistic, digits = digits, format = "f"))
}

# Difference in means
difPr <- function(var1, var2, digits = 3)
{
  m <- mean(var1, na.rm = TRUE) - mean(var2, na.rm = TRUE)
  if(abs(m) < 1000){
    return(formatC(m, digits = digits, format = "f"))
  } else {
    return(paste0(paste0("\\multicolumn{2}{", formatC(m, digits = 0, format = "f")), "}"))
  }
}

# Within R-squared for FELM
WR2 <- function(model, digits = 3){
  return(round(summary(model)$P.r.squared, digits = digits))
}

# Clean up RD coefficients (from rdrobust) for printing 
# To do - 
# Add other rd packages
RDcoef <- function(model, digits = 3){
  if(model$pv[1] < .01){
    return(paste0(round(model$Estimate[1], digits = 3),"^{***}"))
  } else if(model$pv[1] < .05) {
    return(paste0(round(model$Estimate[1], digits = 3),"^{**}"))
  } else if(model$pv[1] < .1) {
    return(paste0(round(model$Estimate[1], digits = 3),"^{*}"))
  } else {
    return(paste0(round(model$Estimate[1], digits = 3),"^{}"))
  }
}

# Clean up RD ses (from rdrobust) for printing 
# To do - 
# Add other rd packages
RDse <- function(model, digits = 3){
  if(round(model$se[1],digits = 3) > 0) {
    return(paste0("(", round(model$se[1],digits = 3), ")"))
  } else {
    return(paste0("(<", 1/10^digits, ")"))
  }
}

# Clean up RD p-values (from rdrobust) for printing 
# To do - 
# Add other rd packages
RDp <- function(model, digits = 3){
  if(round(model$se[1],digits = 3) > 0) {
    return(round(model$pv[1],digits = 3))
  } else {
    return(paste0("<", 1/10^digits))
  }
}

# Clean up RD Effective Sample Size (from rdrobust) for printing 
# To do - 
# Add other rd packages
RDnobs <- function(model){
  nobs <- model$N_b[1] + model$N_h[1]
  return(paste0("\\multicolumn{1}{c}{", format(nobs, big.mark=","), "}"))
}

RDbw <- function(model){
  return(round(model$bws[1],digits = 3))
}

##### Get the mode of a data set #######
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

###### Better config for DataExplorer #######
# create_report(data, y = "", config = dconfig)
dconfig <- list(
  "introduce" = list(),
  "plot_str" = list(
    "type" = "diagonal",
    "fontSize" = 35,
    "width" = 1000,
    "margin" = list("left" = 350, "right" = 250)
  ),
  "plot_missing" = list(),
  "plot_histogram" = list(),
  "plot_bar" = list(),
  "plot_correlation" = list("use" = "pairwise.complete.obs"),
  "plot_boxplot" = list(),
  "plot_scatterplot" = list()
)

####### Tic/toc functionality ########
tic <- function(gcFirst = TRUE, type=c("elapsed", "user.self", "sys.self"))
{
  type <- match.arg(type)
  assign(".type", type, envir=baseenv())
  if(gcFirst) gc(FALSE)
  tic = 0
  tic <- proc.time()[type]
  assign(".tic", tic, envir=baseenv())
  invisible(tic)
}

toc <- function()
{
  type <- get(".type", envir=baseenv())
  toc = 0
  toc <- proc.time()[type]
  tic <- get(".tic", envir=baseenv())
  print(toc - tic)
  invisible(toc)
}

####### Add the following code to TexStudio options/build for .rnw file compiling  ########
# user commands to allow for .Rnw files
# "Rscript.exe" -e  "knitr::knit2pdf('%.Rnw')" | pdflatex -synctex=1 -interaction=nonstopmode %.tex | txs:///view-pdf


####### Add scalebox, fix notes, and have variable of interest in stargazer ########
fixStargazer <- function(tab, nmodels, scalesize = .8, nVarInt = 0, file = 0, label = 0)
{
  # tab should be a character vector of captured output from stargazer
  # this version requires you to use booktabs in your preamble of your tex doc
  # adding a scalebox
  scaleChar <- paste("\\scalebox{", scalesize, "}{", sep = "")
  a <- sub("\\begin{tabular}", paste(scaleChar, "\\begin{tabular}", sep = ""), tab, fixed = TRUE)
  b <- sub("\\end{tabular}", "\\end{tabular}}", a, fixed = TRUE)
  
  # making the notes go across the table
  noteChar <- paste("\\multicolumn{", nmodels + 1, "}{l}{\\parbox", sep = "")
  patChar <- paste(" & \\multicolumn{", nmodels, "}{l}{\\parbox", sep = "")
  d <- sub(patChar,  noteChar, b, fixed = TRUE)
  
  # dropping the comment at the top
  d[2:4] <- ""
  
  # adding Var of Interest and control designation
  if(nVarInt == 1 ){
    d[16] <- sub("\\\\[-2.1ex]", "\\\\[-2.1ex] \\textbf{Variable of Interest:} \\\\", d[16], fixed = TRUE)
    d[16 + nVarInt*2] <- paste("\\textbf{Control Variables:} \\\\", d[16 + nVarInt*2])
  }
  
  if(nVarInt > 1 ){
    d[16] <- sub("\\\\[-2.1ex]", "\\\\[-2.1ex] \\textbf{Variables of Interest:} \\\\", d[16], fixed = TRUE)
    d[16 + nVarInt*2] <- paste("\\textbf{Control Variables:} \\\\", d[16 + nVarInt*2])
  }
  
  # getting line endings
  a <- gsub("\\\\", "\\\\ \n", d, fixed = TRUE)
  
  # printing output (if trunctated, go to RStudiio Options, Code, Display, and set the number at the bottom higher.
  # or printing to a file for auto loading
  if(file == 0 & label == 0) cat(a)
  if(file != 0 & label == 0) cat(a, file = file)
  if(file !=0 & label != 0) cat(a, file = paste0(file, label, ".tex"))
}


#### Fix the weird unlabeled dataset problem ####
# This is also fixed by loading the packages in the order from the master file
# The problem occurs when you load HMisc and tidyverse in the wrong order
labelDataset <- function(data) {
  correctLabel <- function(x) {
    
    if(!is.null(attributes(x)$labels)) {
      class(attributes(x)$labels) <- typeof(x)
    }
    return(x)
  }
  for(i in colnames(data)) {
    data[, i] <- correctLabel(data[, i])
  }
  return(data)
}


#### Produce a well-labeled table for sum stats of a treated and non-treated sample ####
twoSampleSumStats <- function(data, 
                              treatvar, 
                              varnames, 
                              fancyvarnames, 
                              basicPlot = FALSE,
                              caption = "", 
                              label = "", 
                              note = "", 
                              treatcolname = "Treated",
                              nontreatcolname = "Not Treated", 
                              scale = 1,
                              parbox = "21cm",
                              file = 0)
{
  # need to add a way to order the variables like varname orders them.
  # need to add a version that works for beamer (dropping cmidrule).
  # need to add a way to save it to a file instead of printing.
  # should add an option for digits
  
  # pull variables from data
  l <- which(names(data) %in% varnames)
  data <- data[,l]
  
  # separete treat vs nontreat samples
  treatdata    <- data[treatvar == 1, ]
  nontreatdata <- data[treatvar == 0, ]
  vars <- names(data)
  
  tab <- vector()
  # top of the table
  tab[1] <- paste0("\\begin{table}[h!] \\centering \n")
  tab[2] <- paste0("\\caption{",caption,"} \n")
  tab[3] <- paste0("\\label{",label,"} \n")
  tab[4] <- paste0("\\scalebox{", scale, "}{\\begin{tabular}{@{\\extracolsep{2pt}}lD{.}{.}{-3} D{.}{.}{-3} D{.}{.}{-3} D{.}{.}{-3} D{.}{.}{-3} D{.}{.}{-3} D{.}{.}{-3} D{.}{.}{-3} D{.}{.}{-3}} \n")
  tab[5] <- paste0("\\hline \n")
  tab[6] <- paste0("\\\\[0.5mm] \n")
  tab[7] <- paste0("Variable & \\multicolumn{2}{c}{Full Sample:} & \\multicolumn{2}{c}{",treatcolname,":} & \\multicolumn{2}{c}{",nontreatcolname,":}  & \\multicolumn{2}{c}{Diff. Means:}\\\\ \n")
  tab[8] <- paste0("\\cmidrule{2-3} \\cmidrule{4-5} \\cmidrule{6-7} \\cmidrule{8-9} \n")
  tab[9] <- paste0("\\\\[-1.8ex] \n")
  tab[10]<- paste0("& \\multicolumn{1}{c}{Mean} & \\multicolumn{1}{c}{St. Dev.} & \\multicolumn{1}{c}{Mean} & \\multicolumn{1}{c}{St. Dev.} & \\multicolumn{1}{c}{Mean} & \\multicolumn{1}{c}{St. Dev.} & \\multicolumn{1}{c}{Diff.} &\\multicolumn{1}{c}{t}\\\\ \n")
  tab[11]<- paste0("\\hline \\\\[-1.8ex] \n")
  
  row = length(tab)
  # gen sum stats
  for(i in 1:length(vars)){
    fm  <- meanPr(data[[i]])
    fs  <- sdPr(data[[i]])
    tm  <- meanPr(treatdata[[i]])
    ts  <- sdPr(treatdata[[i]])
    ntm <- meanPr(nontreatdata[[i]])
    nts <- sdPr(nontreatdata[[i]])
    d   <- difPr(treatdata[[i]], nontreatdata[[i]])
    dm  <- tdif(treatdata[[i]], nontreatdata[[i]])
    
    row <- row + 1
    tab[row] <- paste(fancyvarnames[i],"&",fm,"&",fs,"&",tm,"&",ts,"&",ntm,"&",nts,"&",d,"&", dm,"\\\\ \n")
  }
  
  # bottom
  tab[row+1] <- paste0("N & \\multicolumn{2}{c}{",nrow(data),"} & \\multicolumn{2}{c}{",nrow(treatdata),"} & \\multicolumn{2}{c}{",nrow(nontreatdata),"} & \\multicolumn{2}{c}{-} \\\\ \n")
  tab[row+2] <- paste0("\\hline \\\\[-1.8ex] \n")
  if (note != "") {
    tab[row+3] <- paste0("\\multicolumn{9}{l}{\\parbox[t]{",parbox,"}{\\textbf{Note:}",note,"}} \\\\ \n")
  }  else {
    tab[row+3] <- paste0(" \n")
  }
  tab[row+4] <- paste0("\\end{tabular}} \n")
  tab[row+5] <- paste0("\\end{table} \n")
  
  
  if(basicPlot == TRUE) {
    tab <- tab[10:(length(tab)-3)]
  }
  
  # print or save table
  if(file == 0) cat(tab)
  if(file !=0) cat(tab, file = paste0(file, label, ".tex"))
}


#### Function to make FELM formula from strings ####
felmForm <- function(y, treat, x, fe = 0, inst = 0, clus = 0)
{
  
  
  if(length(x) > 0){
    form <- paste(y, "~", paste(treat, collapse = "+"), "+", paste(x, collapse = " + "), "|", paste(fe, collapse = " + "), "|", inst, "|", paste(clus, collapse = " + "))
  } else {
    form <- paste(y, "~", paste(treat, collapse = "+"), "|", paste(fe, collapse = " + "), "|", inst, "|", paste(clus, collapse = " + "))
  }
  form <- formula(form, env = globalenv())
  return(form)
}














