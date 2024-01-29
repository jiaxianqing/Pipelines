library(ggplot2)
library(reshape2)
library(ggthemes)
library(ggsci)

# tdata <- read.table("read_count.FC2.SA.fpkm_avr1.c06.txt", header = T, row.names= 1, sep = "\t")
# aa <- scale(tdata)
# 
# write.table(aa,"read_count.FC2.SA.fpkm_avr1.c06.z_score.txt",quote=F, row.names=T, col.names=T, sep="\t")

rm(list=ls())
adata <- read.table("all_deg.fc1_5.T1_T2.fpkm.avr.zscore.txt", header = T, sep = "\t")
fdata <- melt(adata, id.vars = c('ID','Cluster'), variable.name = 'Treat', value.name = 'Value')

write.table(fdata,"all_deg.fc1_5.T1_T2.fpkm.avr.zscore.melt.txt",
            quote=F,row.names=F,col.names=T,sep="\t")

aa <- read.table("clipboard", header = T, sep = "\t")
faa <- melt(aa, id.vars = c('Cluster'), variable.name = 'Treat', value.name = 'Value')
write.table(faa,"aa.mod.txt", quote=F,row.names=F,col.names=T,sep="\t")


ffdata <- read.table("all_deg.fc1_5.T1_T2.fpkm.avr.zscore.melt.mod.txt", header = T, sep = "\t")


ggplot(ffdata) +
  geom_line(aes(x = Treat, y = Value, group =ID), color= "grey") +
  facet_grid(Cluster~Group, scales = "free", space="free")+
  geom_line(data=aa, aes(x=Treat, y=Value, group = Cluster, color = Group), size = 0.8)+
  theme_few()+
  scale_color_lancet()+
  theme(legend.position="none",strip.background = element_rect(fill = "transparent"),
        axis.text.x = element_text(angle = 45, hjust = 0.4, color = "black", vjust = 0.5),
        axis.text.y = element_text(color = "black"))

ggsave(filename = "all_deg.fc1_5.T1_T2.fpkm.avr.zscore.melt.mod.pdf", heigh = 15, width =3.5, device = "pdf")

