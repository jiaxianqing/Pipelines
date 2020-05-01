#########################################################################
input          <- commandArgs(TRUE)[1]
centrome_file  <- commandArgs(TRUE)[2]
genome_file    <- commandArgs(TRUE)[3]
v005           <- commandArgs(TRUE)[4]
v001           <- commandArgs(TRUE)[5]
ylabname       <- commandArgs(TRUE)[6]
titlename      <- commandArgs(TRUE)[7]
output         <- commandArgs(TRUE)[8]
#########################################################################

library(ggplot2)
library(ggthemes)

v005 <- as.numeric(v005)
v001 <- as.numeric(v001)

df <- read.table(input, sep = "\t", header = T)
centrome <- read.table(centrome_file, sep = "\t", col.names = c("Chrom", "start", "end"))
genome <- read.table(genome_file, sep = "\t", col.names = c("Chrom", "length"))
df.flt <- subset(df, Site_count >0)

ggplot(df.flt) +
  geom_rect(data = centrome, aes(xmin=start/1000000, xmax=end/1000000, ymax=max(df$Stats), ymin=0), fill="grey")+
  geom_line(aes(x = (Start + End)/2/1000000, y = Stats), size = 0.5, colour="black") + 
  #geom_line(aes(x = Bin, y = Stats), size = 0.5, colour="black") + 
  facet_grid(Chrom~.)+
  theme_bw() + 
  scale_x_continuous(expand=c(0,0), breaks = seq(0, max(genome$length)/1000000, 5)) +
  #scale_x_continuous(expand=c(0,0)) +
  scale_fill_discrete(guide=FALSE)+
  labs(x="Chromosome lLength (Mbp)",
       y=ylabname,
       title=titlename) +
  geom_hline(aes(yintercept=v005), colour="#990000", linetype="dashed", size = 0.5) + 
  geom_hline(aes(yintercept=v001), colour="darkblue", linetype="dashed", size = 0.5) + 
  theme(legend.position="none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.text=element_text(face="bold",size=rel(1)))
ggsave(filename = output, heigh = 12, width =8, device = "pdf")
