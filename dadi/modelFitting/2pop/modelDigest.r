## a script to digest my model-fitting results
## look through each file and spit out the best-fit parameters for each file


datasets <- c("cp1-cp2", "cp1-cp4", "cp2-cp1", "cp2-cp4")

for (set in datasets) {

  print(set) 
    # to track progress
  
  write.table(paste("\n\n\n", set), file = "model_fits.txt", append = TRUE, quote = FALSE, sep = "")
  # write the dataset name to the output file
  
  models <- dir(set, pattern =".fit") # nsl data
    # get the names of all the files in this dataset
  
  for (model in models) {
    
    print(model)
      # to track progress
    
    data <- read.table(paste(set, "/" , model, sep = ""), header = TRUE)
      # read in the model
    
    data <-  data[!grepl("ll", data[,1]),]
      # scrub extra text out of the dataframe
    
    best <- tail(data[order(as.numeric(as.character(data[,1]))),], 1) 
      # assuming the loglikelihood (ll) value is first col
      # get the row with the highest (least neg) ll value
    
    write.table(best, file = "model_fits.txt", append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE)
      # write this model's best results to the output file
    
  } 
  
}