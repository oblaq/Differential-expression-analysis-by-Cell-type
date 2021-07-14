library(ggplot2)
library(dplyr)


AB_MTG_exon_DNMgene_normalized_data<-read.delim("AB_MTG_exon_normalized_data.txt")
DNMgene_entrezID = read.delim(file="entrez_convert.txt",header=T)[,c(1,2)]
colnames(DNMgene_entrezID)=c("entrezID","gene_name")
AB_MTG_sample_columns =read.csv("human_MTG_2018-06-14_samples-columns.csv")



AB_MTG_exon_DNMgene_normalized_matrix = as.matrix(AB_MTG_exon_DNMgene_normalized_data)

AB_MTG_exon_DNMgene_t = t(AB_MTG_exon_DNMgene_normalized_matrix)
colnames(AB_MTG_exon_DNMgene_t) = DNMgene_entrezID$gene_name[match(colnames(AB_MTG_exon_DNMgene_t),DNMgene_entrezID$entrezID)]

AB_MTG_sample_columns$cell_type=gsub( "\\s.*", "",AB_MTG_sample_columns$cluster)
AB_MTG_sample_columns_ = AB_MTG_sample_columns[,colnames(AB_MTG_sample_columns) %in% c("sex","age_days","brain_subregion","class","cell_type")]

AB_MTG_exon_group= cbind(AB_MTG_sample_columns_,AB_MTG_exon_DNMgene_t)
AB_MTG_exon_group=AB_MTG_exon_group[!is.na(colnames(AB_MTG_exon_group))]

ten_gene=c("ZNF595","SDK1","ZNF718","DPYSL2","KCNQ1","gcnt2","SNX9","MSI2","NBPF1","AACS")

# for (i in 1:10){
#   target_gene=ten_gene[i]
#   print(target_gene)
#   groupup=AB_MTG_exon_group %>% group_by(class) %>% summarise(gene_mean = mean(target_gene),gene_sd = sd(target_gene))
#   group_celltype=AB_MTG_exon_group %>% group_by(cell_type) %>% summarise(gene_mean = mean(target_gene),gene_sd = sd(target_gene))
#   
#   p1.title=paste(target_gene,"(by cell-type)",sep=" ")
#   p1=ggplot(group_celltype, aes(x=cell_type, y=gene_mean)) + 
#     geom_bar(stat="identity") +
#     geom_errorbar(aes(ymin=gene_mean-gene_sd,ymax=gene_mean+gene_sd),
#                   width=.2) +
#     ggtitle(p1.title)
#   png(filename=paste(p1.title,".png",sep=""),width=1000,height=1000)
#   p1
#   dev.off()
#   
#   
#   
#   p2.title=paste(target_gene,"(by class)",sep=" ")
#   #by class
#   p2=ggplot(groupup, aes(x=class, y=gene_mean)) + 
#     geom_bar(stat="identity") +
#     geom_errorbar(aes(ymin=gene_mean-gene_sd,ymax=gene_mean+gene_sd),
#                   width=.2)+
#     ggtitle(p2.title)
#   png(filename=paste(p2.title,".png",sep=""),width=1000,height=1000)
#   p2
#   dev.off()
# }


#ZNF595 SDK1 ZNF718 KCNQ1 gcnt2 SNX9 MSI2 NBPF1 AACS


#dpysl2
groupup=AB_MTG_exon_group %>% group_by(class) %>% summarise(gene_mean = mean(DPYSL2),gene_sd = sd(DPYSL2))
group_celltype=AB_MTG_exon_group %>% group_by(cell_type) %>% summarise(gene_mean = mean(DPYSL2),gene_sd = sd(DPYSL2))
p1=ggplot(group_celltype, aes(x=cell_type, y=gene_mean)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=gene_mean-gene_sd,ymax=gene_mean+gene_sd),
                width=.2) +
  ggtitle("DPYSL2 - cell type")
png(filename="DPYSL2_celltype.png",width=1000,height=1000)
p1
dev.off()

#by class
p2=ggplot(groupup, aes(x=class, y=gene_mean)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=gene_mean-gene_sd,ymax=gene_mean+gene_sd),
                width=.2) +
  ggtitle("DPYSL2 - class")
png(filename="DPYSL2_class.png",width=1000,height=1000)
p2
dev.off()



#ZNF595
groupup=AB_MTG_exon_group %>% group_by(class) %>% summarise(gene_mean = mean(ZNF595),gene_sd = sd(ZNF595))
group_celltype=AB_MTG_exon_group %>% group_by(cell_type) %>% summarise(gene_mean = mean(ZNF595),gene_sd = sd(ZNF595))

p1=ggplot(group_celltype, aes(x=cell_type, y=gene_mean)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=gene_mean-gene_sd,ymax=gene_mean+gene_sd),
                  width=.2) +
  ggtitle("ZNF595 - cell type")
png(filename="ZNF595_celltype.png",width=1000,height=1000)
p1
dev.off()

#by class
p2=ggplot(groupup, aes(x=class, y=gene_mean)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=gene_mean-gene_sd,ymax=gene_mean+gene_sd),
                  width=.2) +
  ggtitle("ZNF595 - class")
png(filename="ZNF595_class.png",width=1000,height=1000)
p2
dev.off()

#SDK1
groupup=AB_MTG_exon_group %>% group_by(class) %>% summarise(gene_mean = mean(SDK1),gene_sd = sd(SDK1))
group_celltype=AB_MTG_exon_group %>% group_by(cell_type) %>% summarise(gene_mean = mean(SDK1),gene_sd = sd(SDK1))

p1=ggplot(group_celltype, aes(x=cell_type, y=gene_mean)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=gene_mean-gene_sd,ymax=gene_mean+gene_sd),
                width=.2) +
  ggtitle("SDK1- cell type")
png(filename="SDK1_celltype.png",width=1000,height=1000)
p1
dev.off()

#by class
p2=ggplot(groupup, aes(x=class, y=gene_mean)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=gene_mean-gene_sd,ymax=gene_mean+gene_sd),
                width=.2) +
  ggtitle("SDK1 - class")
png(filename="SDK1_class.png",width=1000,height=1000)
p2
dev.off()

#ZNF718
groupup=AB_MTG_exon_group %>% group_by(class) %>% summarise(gene_mean = mean(ZNF718),gene_sd = sd(ZNF718))
group_celltype=AB_MTG_exon_group %>% group_by(cell_type) %>% summarise(gene_mean = mean(ZNF718),gene_sd = sd(ZNF718))

p1=ggplot(group_celltype, aes(x=cell_type, y=gene_mean)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=gene_mean-gene_sd,ymax=gene_mean+gene_sd),
                width=.2) +
  ggtitle("ZNF718 - by cell type")
png(filename="ZNF718_celltype.png",width=1000,height=1000)
p1
dev.off()

#by class
p2=ggplot(groupup, aes(x=class, y=gene_mean)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=gene_mean-gene_sd,ymax=gene_mean+gene_sd),
                width=.2) +
  ggtitle("ZNF718 - by class")
png(filename="ZNF718_class.png",width=1000,height=1000)
p2
dev.off()

#KCNQ1
groupup=AB_MTG_exon_group %>% group_by(class) %>% summarise(gene_mean = mean(KCNQ1),gene_sd = sd(KCNQ1))
group_celltype=AB_MTG_exon_group %>% group_by(cell_type) %>% summarise(gene_mean = mean(KCNQ1),gene_sd = sd(KCNQ1))

p1=ggplot(group_celltype, aes(x=cell_type, y=gene_mean)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=gene_mean-gene_sd,ymax=gene_mean+gene_sd),
                width=.2) +
  ggtitle("KCNQ1 - by cell type")
png(filename="KCNQ1_celltype.png",width=1000,height=1000)
p1
dev.off()

#by class
p2=ggplot(groupup, aes(x=class, y=gene_mean)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=gene_mean-gene_sd,ymax=gene_mean+gene_sd),
                width=.2) +
  ggtitle("KCNQ1 - by class")
png(filename="KCNQ1_class.png",width=1000,height=1000)
p2
dev.off()

#gcnt2
groupup=AB_MTG_exon_group %>% group_by(class) %>% summarise(gene_mean = mean(gcnt2),gene_sd = sd(gcnt2))
group_celltype=AB_MTG_exon_group %>% group_by(cell_type) %>% summarise(gene_mean = mean(gcnt2),gene_sd = sd(gcnt2))

p1=ggplot(group_celltype, aes(x=cell_type, y=gene_mean)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=gene_mean-gene_sd,ymax=gene_mean+gene_sd),
                width=.2) +
  ggtitle("gcnt2 - by cell type")
png(filename="gcnt2_celltype.png",width=1000,height=1000)
p1
dev.off()

#by class
p2=ggplot(groupup, aes(x=class, y=gene_mean)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=gene_mean-gene_sd,ymax=gene_mean+gene_sd),
                width=.2) +
  ggtitle("gcnt2 - by class")
png(filename="gcnt2_class.png",width=1000,height=1000)
p2
dev.off()





#SNX9
groupup=AB_MTG_exon_group %>% group_by(class) %>% summarise(gene_mean = mean(SNX9),gene_sd = sd(SNX9))
group_celltype=AB_MTG_exon_group %>% group_by(cell_type) %>% summarise(gene_mean = mean(SNX9),gene_sd = sd(SNX9))

p1=ggplot(group_celltype, aes(x=cell_type, y=gene_mean)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=gene_mean-gene_sd,ymax=gene_mean+gene_sd),
                width=.2) +
  ggtitle("SNX9 - by cell type")
png(filename="SNX9_celltype.png",width=1000,height=1000)
p1
dev.off()

#by class
p2=ggplot(groupup, aes(x=class, y=gene_mean)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=gene_mean-gene_sd,ymax=gene_mean+gene_sd),
                width=.2) +
  ggtitle("SNX9 - by class")
png(filename="SNX9_class.png",width=1000,height=1000)
p2
dev.off()

#MSI2
groupup=AB_MTG_exon_group %>% group_by(class) %>% summarise(gene_mean = mean(MSI2),gene_sd = sd(MSI2))
group_celltype=AB_MTG_exon_group %>% group_by(cell_type) %>% summarise(gene_mean = mean(MSI2),gene_sd = sd(MSI2))

p1=ggplot(group_celltype, aes(x=cell_type, y=gene_mean)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=gene_mean-gene_sd,ymax=gene_mean+gene_sd),
                width=.2) +
  ggtitle("MSI2 - by cell type")
png(filename="MSI2_celltype.png",width=1000,height=1000)
p1
dev.off()

#by class
p2=ggplot(groupup, aes(x=class, y=gene_mean)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=gene_mean-gene_sd,ymax=gene_mean+gene_sd),
                width=.2) +
  ggtitle("MSI2 - by class")
png(filename="MSI2_class.png",width=1000,height=1000)
p2
dev.off()


#NBPF1
groupup=AB_MTG_exon_group %>% group_by(class) %>% summarise(gene_mean = mean(NBPF1),gene_sd = sd(NBPF1))
group_celltype=AB_MTG_exon_group %>% group_by(cell_type) %>% summarise(gene_mean = mean(NBPF1),gene_sd = sd(NBPF1))

p1=ggplot(group_celltype, aes(x=cell_type, y=gene_mean)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=gene_mean-gene_sd,ymax=gene_mean+gene_sd),
                width=.2) +
  ggtitle("NBPF1 - by cell type")
png(filename="NBPF1_celltype.png",width=1000,height=1000)
p1
dev.off()

#by class
p2=ggplot(groupup, aes(x=class, y=gene_mean)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=gene_mean-gene_sd,ymax=gene_mean+gene_sd),
                width=.2) +
  ggtitle("NBPF1 - by class")
png(filename="NBPF1_class.png",width=1000,height=1000)
p2
dev.off()




#violin plot
AB_MTG_exon_DNMgene_normalized_matrix_with_pheno = cbind(AB_MTG_exon_DNMgene_normalized_matrix_,AB_MTG_sample_annot)
AB_MTG_exon_DNMgene_normalized_matrix_with_pheno$cluster = factor(AB_MTG_exon_DNMgene_normalized_matrix_with_pheno$cluster,levels = c("GABAergic_Inh_L1","GABAergic_Inh_L1-2","GABAergic_Inh_L1-3","GABAergic_Inh_L1-4","GABAergic_Inh_L2-3","GABAergic_Inh_L2-4","GABAergic_Inh_L2-5","GABAergic_Inh_L2-6","GABAergic_Inh_L3-5","GABAergic_Inh_L3-6","GABAergic_Inh_L4-5","GABAergic_Inh_L4-6","GABAergic_Inh_L5-6","Glutamatergic_Exc_L2","Glutamatergic_Exc_L2-3","Glutamatergic_Exc_L2-4","Glutamatergic_Exc_L3-4","Glutamatergic_Exc_L3-5","Glutamatergic_Exc_L4-5","Glutamatergic_Exc_L4-6","Glutamatergic_Exc_L5-6","Glutamatergic_Exc_L6","Non-neuronal_Astro_L1-2","Non-neuronal_Astro_L1-6","Non-neuronal_Endo_L2-6","Non-neuronal_Micro_L1-3","Non-neuronal_Oligo_L1-6","Non-neuronal_OPC_L1-6","no class"))

p_exon = list()
p_exon[[1]] = ggplot(AB_MTG_exon_DNMgene_normalized_matrix_with_pheno, aes(x=cluster,y=GJC1)) + geom_violin() + geom_jitter(shape=16, position = position_jitter(0.2)) + coord_flip()


groupup=AB_MTG_exon_group %>% group_by(class) %>% summarise(gene_mean = mean(AACS),gene_sd = sd(AACS))
group_celltype=AB_MTG_exon_group %>% group_by(cell_type) %>% summarise(gene_mean = mean(AACS),gene_sd = sd(AACS))

p1=ggplot(group_celltype, aes(x=cell_type, y=gene_mean)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=gene_mean-gene_sd,ymax=gene_mean+gene_sd),
                width=.2) +
  ggtitle("AACS - by cell type")
png(filename="AACS_celltype.png",width=1000,height=1000)
p1
dev.off()

#by class
p2=ggplot(groupup, aes(x=class, y=gene_mean)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=gene_mean-gene_sd,ymax=gene_mean+gene_sd),
                width=.2) +
  ggtitle("AACS - by class")
png(filename="AACS_class.png",width=1000,height=1000)
p2
dev.off()