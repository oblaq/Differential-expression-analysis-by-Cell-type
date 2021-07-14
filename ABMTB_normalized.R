library(aws.s3)
library(dplyr)

Sys.setenv("AWS_ACCESS_KEY_ID" = "AKIAQQFS7Y5XOBCCACW2",
           "AWS_SECRET_ACCESS_KEY" = "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",
           "AWS_DEFAULT_REGION" = "us-west-2")
get_bucket("brainscrna")
AB_MTG_exon_DNMgene=s3read_using(FUN = read.csv, bucket = "brainscrna", object = "human_MTG_2018-06-14_exon-matrix.csv")
AB_MTG_sample_columns =s3read_using(FUN = read.csv, bucket = "brainscrna", object = "human_MTG_2018-06-14_samples-columns.csv")
DNMgene_entrezID = read.delim(file="entrez_convert.txt",header=T)[,c(1,2)]
colnames(DNMgene_entrezID)=c("entrezID","gene_name")

rownames(AB_MTG_exon_DNMgene)=AB_MTG_exon_DNMgene[,1]
AB_MTG_exon_DNMgene=AB_MTG_exon_DNMgene[,-1]


library(Seurat)
AB_MTG_exon = CreateSeuratObject(counts = AB_MTG_exon_DNMgene,project = "AB_MTG_exon_snRNAseq", min.cells = 0, min.features = 0)
AB_MTG_exon = NormalizeData(AB_MTG_exon,normalization.method = "LogNormalize", scale.factor = 10000)
AB_MTG_exon_normalized_data = AB_MTG_exon[["RNA"]]@data


AB_MTG_exon_DNMgene_normalized_data = AB_MTG_exon_normalized_data[rownames(AB_MTG_exon_normalized_data) %in% DNMgene_entrezID$entrezID,]
AB_MTG_exon_DNMgene_normalized_matrix = as.matrix(AB_MTG_exon_DNMgene_normalized_data)

AB_MTG_exon_DNMgene_t = t(AB_MTG_exon_DNMgene_normalized_matrix)
colnames(AB_MTG_exon_DNMgene_t) = DNMgene_entrezID$gene_name[match(colnames(AB_MTG_exon_DNMgene_t),DNMgene_entrezID$entrezID)]

AB_MTG_sample_columns$cell_type=gsub( "\\s.*", "",AB_MTG_sample_columns$cluster)
AB_MTG_sample_columns_ = AB_MTG_sample_columns[,colnames(AB_MTG_sample_columns) %in% c("sex","age_days","brain_subregion","class","cell_type")]

AB_MTG_exon_group= cbind(AB_MTG_sample_columns_,AB_MTG_exon_DNMgene_t)
AB_MTG_exon_group=AB_MTG_exon_group[!is.na(colnames(AB_MTG_exon_group))]

groupup=AB_MTG_exon_group %>% group_by(class) %>% summarise(DPYSL2_mean = mean(DPYSL2),DPYSL2_sd = sd(DPYSL2))
group_celltype=AB_MTG_exon_group %>% group_by(cell_type) %>% summarise(DPYSL2_mean = mean(DPYSL2),DPYSL2_sd = sd(DPYSL2))

#write to s3
# get_bucket(celltyperesult)
s3write_using(groupup, FUN = write.table,quote=FALSE,sep="\t",row.names=FALSE,col.name=TRUE,object = "class_summarise.txt", bucket = "celltyperesult")
s3write_using(group_celltype, FUN = write.table,quote=FALSE,sep="\t",row.names=FALSE,col.name=TRUE,object = "celltype_summarise.txt", bucket = "celltyperesult")

# s3write_using(AB_MTG_exon_normalized_data, FUN = write.table,quote=FALSE,sep="\t",row.names=TRUE,col.name=TRUE,object = "AB_MTG_exon_normalized_data.txt", bucket = "celltyperesult")

write.table(AB_MTG_exon_normalized_data, file="AB_MTG_exon_normalized_data_v2.txt",quote=FALSE,sep="\t",row.names=TRUE,col.name=TRUE)


put_object("/home/rstudio/AB_MTG_exon_normalized_data_v2.txt", object = "AB_MTG_exon_normalized_data_v2.txt", bucket = "celltyperesult",show_progress = TRUE, multipart = TRUE)




library(ggplot2)
p1=ggplot(group_celltype, aes(x=cell_type, y=DPYSL2_mean)) + 
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=DPYSL2_mean-DPYSL2_sd,ymax=DPYSL2_mean+DPYSL2_sd),
                width=.2)



#no error bar
p2=ggplot(group_celltype, aes(x=cell_type, y=DPYSL2_mean)) + 
  geom_bar(stat="identity")

#by class
p3=ggplot(groupup, aes(x=class, y=DPYSL2_mean)) + 
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=DPYSL2_mean-DPYSL2_sd,ymax=DPYSL2_mean+DPYSL2_sd),
                width=.2)
