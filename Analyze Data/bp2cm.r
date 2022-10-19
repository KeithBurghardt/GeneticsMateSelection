# interpolation inspired by https://www.biostars.org/p/16236/
# Note we need a marey map to convert from bp to cm
# MM uses ~1000 known relations (note these differ from men and women; we average the differences for now)
# interpolate with a number of methods, I use spline, use "plot(map)" to show fit

library(MareyMap)
library(stats)
data(Homo_sapiens_mean)


for (ch in c(1:22)){
        padded_chr <-  sprintf("Chromosome %02d", ch)
        in_file <- sprintf("/project/arpawong_181/HRS_AsMa/keith/phased/chr%d.map",ch)
        out_file <- sprintf("/project/arpawong_181/HRS_AsMa/keith/phased/interp_chr%d.map",ch)
	chr_data <- read.table(in_file,header=FALSE,sep="\t")
	human_chr <- Homo_sapiens_mean[[padded_chr]]
	itr <- MMSpline3()
	human_chr <- human_chr + itr
	interpolation <- (interpolations(human_chr)$spline)@model
        #print(predict(interpolation,chr_data[,4]))
	chr_data[,3] <- predict(interpolation,chr_data[,4])$y
        write.table(chr_data, file = out_file,row.names = FALSE,col.names=FALSE,quote = FALSE,sep="\t")
}
