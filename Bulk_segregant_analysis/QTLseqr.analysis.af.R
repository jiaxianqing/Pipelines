#
# Xianqing Jia
# 2019-01-07
#
# https://github.com/bmansfeld/QTLseqr
# https://dl.sciencesocieties.org/publications/tpg/abstracts/11/2/180006
# 
# devtools::install_github("bmansfeld/QTLseqr")

#########################################################################
filename       <- commandArgs(TRUE)[1]
HighBulk       <- commandArgs(TRUE)[2]
LowBulk        <- commandArgs(TRUE)[3]
winsize        <- commandArgs(TRUE)[4]
outfile       <- commandArgs(TRUE)[5]
#########################################################################

library("QTLseqr")
library("ggplot2")

outfile1 <- paste(outfile,".QTLs.deltaSNP.csv", sep="")
outfile2 <- paste(outfile,".QTLs.Gprime.csv", sep="")
outfile3 <- paste(outfile,".sites_info.csv", sep="")
outfile4 <- paste(outfile,".QTLs.log", sep="")

winsize <- as.numeric(winsize)
Chroms <- c("chr01", "chr02", "chr03", "chr04", "chr05", "chr06", "chr07", "chr08", "chr09", "chr10", "chr11", "chr12")
sink(outfile4)

df <- importFromTable(
  file = filename,
  highBulk = HighBulk,
  lowBulk = LowBulk,
  chromList = Chroms,
  sep = "\t"
)

df_filt <- filterSNPs(
  SNPset = df,
  refAlleleFreq = 0.20,
  minTotalDepth = 100,
  maxTotalDepth = 400,
  minSampleDepth = 40,
  minGQ = 99
)

#Run QTLseq analysis
df_filt <- runQTLseqAnalysis(
  SNPset = df_filt,
  windowSize = winsize,
  popStruc = "F2",
  bulkSize = c(14, 15),
  replications = 10000,
  intervals = c(95, 99),
  maxk = 10000
)
#Run G' analysis
df_filt <- runGprimeAnalysis(
  SNPset = df_filt,
  windowSize = winsize,
  outlierFilter = "deltaSNP",
  filterThreshold = 0.1,
  maxk = 10000
)


getQTLTable(SNPset = df_filt, method = "deltaSNP", alpha = 0.01, export = TRUE, fileName = outfile1)
getQTLTable(SNPset = df_filt, method = "Gprime", alpha = 0.01, export = TRUE, fileName = outfile2)
write.table(df[], file = outfile3, sep = "\t",quote = F)

sink()