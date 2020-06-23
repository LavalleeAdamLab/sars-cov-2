setwd("/Users/Rachel/eclipse-files/network2020/mcl_analysis/ontologizer_files/")


go_data <- read.csv("2020.04.06.kroganMCLi2_allGOs.txt", sep = "\t", header = F)
go_data$p_adj <- p.adjust(go_data$V2, method = "fdr")

for(i in (c(0.01,0.001))){
  significantGO <- go_data[go_data$p_adj < i, c(1,3)]
  sortedGO <- significantGO[order(significantGO$V1, abs(significantGO$p_adj)),]
  finalGO <- sortedGO[ !duplicated(sortedGO$V1), ]
  
  output_file <- paste("2020.04.06.kroganMCLi2_SignificantCorrectedGO_p", i, ".txt", sep="")
  readr::write_tsv(finalGO, output_file, quote=F)  
}

