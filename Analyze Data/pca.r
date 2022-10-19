library(devtools)
dev_mode(on=T)
library('SNPRelate')
dir<-"/project/arpawong_181/HRS_AsMa/keith/phased2/"
dirout<-"/project/burghard_687/genetic_data/phased2/"
gds_files = c()
for (ch in c(1:22)){
        print(ch)
        ped_f = paste0(paste(dir,"chr",sep=""),toString(ch),".ped")
        map_f = paste0(paste(dir,"chr",sep=""),toString(ch),"_imputed.map")
        gds_f = paste0(paste(dirout,"chr",sep=""),toString(ch),".gds")
        #gds_files[[ch]] <- gds_f
        gds_files <- append(gds_files,gds_f)

        #if (!file.exists(gds_f)){
        snpgdsPED2GDS( ped_f, map_f, gds_f)
        #}
}


print(gds_files)
full_gds = paste(dirout,"all_chr.gds",sep="")
if (!file.exists(full_gds)){
    snpgdsCombineGeno(gds_files, full_gds)
}


gds = snpgdsOpen(full_gds)
pca = snpgdsPCA(gds)
out_file=paste(dirout,"pca_phased",sep="")
save(pca, file = out_file)
ids = pca$sample.id
evals = pca$eigenval
evecs = pca$eigenvect
var_expl = pca$varprop
print(ids)
print(var_expl)
evals_f = paste(dir,"all_chr_gds_evals.csv",sep="")
write.table(evals, file = evals_f,row.names = FALSE,col.names=FALSE,quote = FALSE,sep=",")

ids_f = paste(dir,"all_ids.csv",sep="")
write.table(ids, file = ids_f,row.names = FALSE,col.names=FALSE,quote = FALSE,sep=",")

var_expl_f = paste(dir,"all_chr_gds_var_expl.csv",sep="")
write.table(var_expl, file = var_expl_f,row.names = FALSE,col.names=FALSE,quote = FALSE,sep=",")

evecs_f = paste(dir,"all_chr_gds_evecs.csv",sep="")
write.table(evecs, file = evecs_f,row.names = FALSE,col.names=FALSE,quote = FALSE,sep=",")
