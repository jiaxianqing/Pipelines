##############################################################
# R code for QTL mapping
#
# http://www.rqtl.org
#
# 2018-1-25
##############################################################
#########################################################################
wkdir       <- commandArgs(TRUE)[1]
in_file     <- commandArgs(TRUE)[2]
outhk       <- commandArgs(TRUE)[3]
outem       <- commandArgs(TRUE)[4]
outimp      <- commandArgs(TRUE)[5]
#########################################################################

library("qtl")
setwd(wkdir)
#read data
hd <- read.cross("csvr", genotypes=c("AA","AB","BB"), alleles=c("A", "B"), dir = wkdir, in_file)

###################### summary of raw data ######################
# stats
sink("samples.summary_raw_data.txt")
print("#summary of raw data")
summary(hd)
sink()

hd <- calc.genoprob(hd, step=1)

##count CO number
# nxo <- countXO(hd)
# pdf(file = "crossover.count.pdf")
# plot(nxo, ylab="No. crossovers")
# dev.off()

######################## single-QTL ########################
out.em <- scanone(hd, method="em")
out.hk <- scanone(hd, method="hk")
out.imp <- scanone(hd, method="imp")

write.table(out.hk[], file = outhk, sep = "\t",quote = F)
write.table(out.hk[], file = outem, sep = "\t",quote = F)
write.table(out.hk[], file = outimp, sep = "\t",quote = F)
# stats
sink("samples.single.qtl.txt")
print("#summary of out.em")
summary(out.em)
print("#summary of out.hk")
summary(out.hk)
print("#summary of out.imp")
summary(out.imp)
sink()
#plot
pdf(file = "out.em.pdf")
plot(out.em)
dev.off()
pdf(file = "out.hk.pdf")
plot(out.hk)
dev.off()
pdf(file = "out.imp.pdf")
plot(out.imp)
dev.off()

#pdf(file = "out.hk.chr01.pdf")
#plot(out.hk, chr="01")
#dev.off()

###################### Permutation tests ######################
operm <- scanone(hd, method="hk", n.perm=1000)
#stats
sink("Permutation.tests.txt")
print("#summary of operm")
summary(operm)
print("#summary of operm, 0.01 and 0.05")
summary(operm, alpha=c(0.01, 0.05))
print("#summary of operm, 1% significance level")
summary(out.hk, perms=operm, alpha=0.01, pvalues=TRUE)
print("#summary of operm, 5% significance level")
summary(out.hk, perms=operm, alpha=0.05, pvalues=TRUE)
sink()

#plot
pdf(file = "operm.pdf")
plot(operm)
dev.off()
###################### Interval estimates of QTL location ######################
sink("Interval.estimates.hk.1.8_LOD.ex2marker.txt")
print("#LOD support intervals, 1.8-LOD, expand to marker")
lodint(out.hk, chr="01", drop=1.8, expandtomarkers=TRUE)
lodint(out.hk, chr="02", drop=1.8, expandtomarkers=TRUE)
lodint(out.hk, chr="03", drop=1.8, expandtomarkers=TRUE)
lodint(out.hk, chr="04", drop=1.8, expandtomarkers=TRUE)
lodint(out.hk, chr="05", drop=1.8, expandtomarkers=TRUE)
lodint(out.hk, chr="06", drop=1.8, expandtomarkers=TRUE)
lodint(out.hk, chr="07", drop=1.8, expandtomarkers=TRUE)
lodint(out.hk, chr="08", drop=1.8, expandtomarkers=TRUE)
lodint(out.hk, chr="09", drop=1.8, expandtomarkers=TRUE)
lodint(out.hk, chr="10", drop=1.8, expandtomarkers=TRUE)
lodint(out.hk, chr="11", drop=1.8, expandtomarkers=TRUE)
lodint(out.hk, chr="12", drop=1.8, expandtomarkers=TRUE)
sink()
sink("Interval.estimates.hk.Bayes_0.95.ex2marker.txt")
print("#Bayes credible intervals, 95%, expand to marker")
bayesint(out.hk, chr="01", prob=0.95, expandtomarkers=TRUE)
bayesint(out.hk, chr="02", prob=0.95, expandtomarkers=TRUE)
bayesint(out.hk, chr="03", prob=0.95, expandtomarkers=TRUE)
bayesint(out.hk, chr="04", prob=0.95, expandtomarkers=TRUE)
bayesint(out.hk, chr="05", prob=0.95, expandtomarkers=TRUE)
bayesint(out.hk, chr="06", prob=0.95, expandtomarkers=TRUE)
bayesint(out.hk, chr="07", prob=0.95, expandtomarkers=TRUE)
bayesint(out.hk, chr="08", prob=0.95, expandtomarkers=TRUE)
bayesint(out.hk, chr="09", prob=0.95, expandtomarkers=TRUE)
bayesint(out.hk, chr="10", prob=0.95, expandtomarkers=TRUE)
bayesint(out.hk, chr="11", prob=0.95, expandtomarkers=TRUE)
bayesint(out.hk, chr="12", prob=0.95, expandtomarkers=TRUE)
sink()
