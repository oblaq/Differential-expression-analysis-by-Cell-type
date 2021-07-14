library(aws.s3)
library(dplyr)
Sys.setenv("AWS_ACCESS_KEY_ID" = "AKIAQQFS7Y5XOBCCACW2",
           "AWS_SECRET_ACCESS_KEY" = "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",
           "AWS_DEFAULT_REGION" = "us-west-2")
get_bucket("brainscrna")
AB_MTG_exon_DNMgene=s3read_using(FUN = read.csv, bucket = "brainscrna", object = "human_MTG_2018-06-14_exon-matrix.csv")
AB_MTG_intron_DNMgene=s3read_using(FUN = read.csv, bucket = "brainscrna", object = "human_MTG_2018-06-14_intron-matrix.csv")

AB_MTG_exon_DNMgene_t = t(AB_MTG_exon_DNMgene)
colnames(AB_MTG_exon_DNMgene_t)=AB_MTG_exon_DNMgene_t[1,]
AB_MTG_exon_DNMgene_t=AB_MTG_exon_DNMgene_t[-1,]

AB_MTG_intron_DNMgene_t = t(AB_MTG_intron_DNMgene)
colnames(AB_MTG_intron_DNMgene_t)=AB_MTG_intron_DNMgene_t[1,]
AB_MTG_intron_DNMgene_t=AB_MTG_intron_DNMgene_t[-1,]

AB_MTG_sample_columns =s3read_using(FUN = read.csv, bucket = "brainscrna", object = "human_MTG_2018-06-14_samples-columns.csv")
DNMgene_entrezID = read.delim(file="entrez_convert.txt",header=T)[,c(1,2)]
colnames(DNMgene_entrezID)=c("entrezID","gene_name")


##rename colnames based on DNMgene_entrezID
colnames(AB_MTG_exon_DNMgene_t) = DNMgene_entrezID$gene_name[match(colnames(AB_MTG_exon_DNMgene_t),DNMgene_entrezID$entrezID)]
colnames(AB_MTG_intron_DNMgene_t) = DNMgene_entrezID$gene_name[match(colnames(AB_MTG_intron_DNMgene_t),DNMgene_entrezID$entrezID)]

AB_MTG_sample_columns_ = AB_MTG_sample_columns[,colnames(AB_MTG_sample_columns) %in% c("sex","age_days","brain_subregion","class")]

AB_MTG_exon_DNMgene= cbind(AB_MTG_sample_columns_,AB_MTG_exon_DNMgene_t)
AB_MTG_intron_DNMgene = cbind(AB_MTG_sample_columns_,AB_MTG_intron_DNMgene_t)

AB_MTG_exon_DNMgene=AB_MTG_exon_DNMgene[!is.na(names(AB_MTG_exon_DNMgene))]
AB_MTG_intron_DNMgene=AB_MTG_intron_DNMgene[!is.na(names(AB_MTG_intron_DNMgene))]


AB_MTG_exon_DNMgene_split = split(AB_MTG_exon_DNMgene,AB_MTG_exon_DNMgene$class)
AB_MTG_intron_DNMgene_split = split(AB_MTG_intron_DNMgene,AB_MTG_intron_DNMgene$class)

groupup=AB_MTG_exon_DNMgene %>% group_by(class) %>% summarise(DPYSL2_mean = mean(DPYSL2),DPYSL2_sd = sd(DPYSL2))

library(ggplot2)
ggplot(groupup, aes(x=class, y=DPYSL2_mean)) + 
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=DPYSL2_mean-DPYSL2_sd,ymax=DPYSL2_mean+DPYSL2_sd),
                width=.2)