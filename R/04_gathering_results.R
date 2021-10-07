## This file is used to gather the results from simulations

#setwd("/Users/gr8lawrence/Documents/Dissertation/conv_opt_code/cleaned_algorithm/PFSM/R")

round_number <- 1

files <- dir("./output/")

res = readRDS(paste("./output/", files[1], sep = ''))
for (i in 2:length(files)) {
  res1 = readRDS(paste("./output/", files[i], sep = '')) # load results files
  res = rbind(res, res1)
}

save(res, file=paste0("./rdata/sim_round_", round_number, ".RData"))

