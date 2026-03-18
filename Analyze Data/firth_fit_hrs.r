library(dplyr)
library(base)
library(fastDummies)
library(logistf)
library(data.table)
library(caret)

args <- commandArgs(trailingOnly = TRUE)
boot="False"
if (length(args) > 0){
	select_slice = args[2]
	train80 = args[4]
        one_ethnic = "False"
	fold=-1
        if(args[5]=='--boot'){boot = args[6]}
	else{if (length(args) > 4){
	    fold = args[6]
            one_ethnic = args[8]
	}}
}else{
	select_slice = "None"
	train80 = "True"
        one_ethnic = "False"
	fold = 0
}

one_eth <- ""
if (one_ethnic == "True"){
one_eth <- "_one_ethnic"
}

# directory to collect/store data
dir <- "/project/arpawong_181/burghard_687/genetic_data/"
train_outcome <- read.csv(paste(dir,"train_outcome_all_data_fold=-1_cog27.csv",sep=""))
train_X <- read.csv(paste(dir,"train_X_non_slice_features_all_data_fold=-1_cog27",".csv",sep=""))
num_marriages = length(subset(train_outcome, marriage == 1)$marriage)

#slice_features = c(select_slice)

library(stringr)

# Get list of matching files
file_list <- list.files("/project/arpawong_181/burghard_687/genetic_data/slices_combined/",
                        pattern = "^hrs_slice_.*\\.csv$", full.names = TRUE)

# Extract and transform filenames
slice_features <- sapply(file_list, function(sf) {
  filename <- basename(sf)
  chr_id <- str_replace_all(filename, c("^hrs_slice_" = "", "\\.csv$" = ""))
  paste0("chr", chr_id)
})
slice_features <- sample(slice_features)

print(slice_features)
for (select_slice in slice_features){
  for (abs_features in c(T)){
    for (s in c("full")){ # ,"random","nom")){
        print(select_slice)
        ######## CLUSTER ##########
        cluster_ids <-read.csv(paste("/project/arpawong_181/burghard_687/genetic_data/slices_combined/hrs_slice_",gsub("chr", "",select_slice),".csv",sep=""))
        if(max(cluster_ids) == -1){next}
        #train_slice <- train_slice[hrs_indices,]
        # cut off dummy variables
        #train_slice[train_slice>30] <- 31
        cluster_ids <- data.frame(cluster_ids)
	colnames(cluster_ids) <- c("newclust")
        # dummify the data
        #dmy <- dummyVars(" ~ .", data = train_slice)
        #train_slice <- data.frame(predict(dmy, newdata = train_slice))
	# make dummy/1-hot variables
        cluster_ids <- dummy_cols(cluster_ids)
        dmy_names <- colnames(cluster_ids)
	# only select dummy names
        dmy_names <- dmy_names[2:length(dmy_names)]
        cluster_ids <- select(cluster_ids, dmy_names)
        colnames(cluster_ids) <- gsub("-","minus",gsub("_","",colnames(cluster_ids)))
	# if (length(colnames(cluster_ids))==1 & colnames(cluster_ids)[1] == "newclust-1"){next}
        ###############################
	# combine all training features
	df <- cbind(cbind(train_X,cluster_ids),train_outcome)
        # abs
        if (abs_features){
		df <- abs(df)
	}
	# outcome variable
	y <- "marriage"
	# create name for slice
	slices <- colnames(df)[grepl("newclust", colnames(df), fixed = TRUE)]
        print(slices)
	# name of non-slice features
	covars_list <- colnames(df)[!grepl("newclust", colnames(df), fixed = TRUE)]
	all_covars <- c("place","ethnicity","relig","edu","height","cog27")#"bmi"
        covar_combos <- c("none",paste(all_covars,collapse="+"))
		print(covar_combos)
	for (covar in covar_combos){
            # file names
	    folder = ""
	    if (abs_features){
		folder = "_abs"
	    }else{
		folder = "_nonabs"
	    }
            coef_file = paste(c( dir,"slices_all",folder,one_eth,"_cog27_",fold,"_new2/","firth_",covar,"_coefs_",select_slice,"_",s,folder,"_no_cog-bmi_cog27_fold=",fold,"_boot=",boot,"_newclusters.csv"),collapse="")
            prob_file = paste(c( dir,"slices_all",folder,one_eth,"_cog27_",fold,"_new2/","firth_",covar,"_prob_",select_slice,"_",s,folder,"_no_cog-bmi_cog27_fold=",fold,"_boot=",boot,"_newclusters.csv"),collapse="")#"-EV+TestSep.csv"),collapse="")
            old_coef_file = paste(c( dir,"slices_all",folder,one_eth,"_cog27_",fold,"_new2/","firth_",covar,"_coefs_",select_slice,"_",s,folder,"_no_cog-bmi_cog27_fold=",fold,".csv"),collapse="")#,"-EV+TestSep.csv"),collapse="")
            old_prob_file = paste(c( dir,"slices_all",folder,one_eth,"_cog27_",fold,"_new2/","firth_",covar,"_prob_",select_slice,"_",s,folder,"_no_cog-bmi_cog27_fold=",fold,".csv"),collapse="")#"-EV+TestSep.csv"),collapse="")

            # ignore analysis if we have already calculated values
            if (file.exists(prob_file) || file.exists(coef_file)){#(!file.exists(old_prob_file) || file.exists(prob_file)){
                print('ignore file')
            	next
            } 
            print(select_slice)
            #print(prob_file)
            # collect covariate features
	    covar_features <- c()
	    if (grepl("place",covar)){
	        covar_features <- c(covar_features,c('RABPLACE_1.0','RABPLACE_2.0','RABPLACE_3.0','RABPLACE_4.0','RABPLACE_5.0','RABPLACE_6.0','RABPLACE_7.0','RABPLACE_8.0','RABPLACE_9.0','RABPLACE_10.0','RABPLACE_11.0'))
	    }
            if (grepl("ethnicity_only",covar)){
                covar_features <- c(covar_features,c('ethnicity_1.0','ethnicity_2.0','ethnicity_3.0','ethnicity_4.0'))
                #for (ev in 1:10){ # first 10 PCA components
                #       covar_features <- c(covar_features,paste("EV",as.character(ev),sep=""))
                #}
            }
            if (grepl("pca_distance",covar)){
                covar_features <- c(covar_features,c('pca_distance'))
                #for (ev in 1:10){ # first 10 PCA components
                #       covar_features <- c(covar_features,paste("EV",as.character(ev),sep=""))
                #}
            }
	    if (grepl("edu_only",covar)){
	        covar_features <- c(covar_features,c('RAEDYRS'))
	    }
	    if (grepl("cog27",covar)){ 
	        covar_features <- c(covar_features,c('cog27'))
	    }
	    if (grepl("income",covar)){ 
	        covar_features <- c(covar_features,c('iearn'))
	    }
	    if (grepl("num_marriage",covar)){ 
	        covar_features <- c(covar_features,c('mrct'))
	    }
	    if (grepl("num_div",covar)){ 
	        covar_features <- c(covar_features,c('mdiv'))
	    }

	    if (grepl("relig",covar)){
	        covar_features <- c(covar_features,c('RARELIG_1.0','RARELIG_2.0','RARELIG_3.0','RARELIG_4.0','RARELIG_5.0'))
	    }
	    #if (grepl("bmi",covar)){
	    #    covar_features <- c(covar_features,c('bmi'))
	    #}
	    if (grepl("height",covar)){
	        covar_features <- c(covar_features,c('height'))
	    }
            if (grepl("ethnicity",covar) & covar != "ethnicity_only"){
		if(grepl("nopca",covar)){
   		    covar_features <- c(covar_features,c('ethnicity_1.0','ethnicity_2.0','ethnicity_3.0','ethnicity_4.0'))
		}else{
	            covar_features <- c(covar_features,c('ethnicity_1.0','ethnicity_2.0','ethnicity_3.0','ethnicity_4.0','pca_distance'))
		}
                #for (ev in 1:10){ # first 10 PCA components
		#	covar_features <- c(covar_features,paste("EV",as.character(ev),sep=""))
		#}
	    }

	    if (grepl("edu",covar) & covar != "edu_only"){
	        covar_features <- c(covar_features,c('RAEDYRS','iearn','mdiv','mrct'))#'cogtot'
	    }
            if(one_ethnic == "True"){
                # if european ethnicity not in covar_features...
		if(!("ethnicity_1.0" %in% covar_features)){
		    covar_features <- c(covar_features,c('ethnicity_1.0'))
		}
	    }
            # create data frame
	    df_coef = data.frame(rbind(c(c("intercept"),covar_features,slices)))
	    colnames(df_coef) <- c(c("intercept"),covar_features,slices)
	    df_prob = data.frame(rbind(c(c("intercept"),covar_features,slices)))
	    colnames(df_prob) <- c(c("intercept"),covar_features,slices)
	    df_lowhigh = data.frame(rbind(c(c("intercept"),covar_features,slices,c("intercept"),covar_features,slices)))
	    colnames(df_lowhigh) <- c(c("intercept"),covar_features,slices,c("intercept"),covar_features,slices)
            #slice_ind <- 0
	    # for (slice in slices) {
            
            if(one_ethnic == "True"){
                    df <- df[df$ethnicity_1.0 == 1,]
	    }
            cleaned_covar_features <- c()

            for (col in c(covar_features)){
		if(sd(df[,col])> 0.000001){
			cleaned_covar_features <- c(cleaned_covar_features,col)
		}
	    }
	    cleaned_slices <- c()
            for (col in c(slices)){
                if( (sd(df[,col])> 0.00001) & (sum(df[,col])>3)){
                print(sum(df[,col]))
		cleaned_slices <- c(cleaned_slices,col)
		}
	    }
            # if there are no viable slices, ignore
	    if (length(cleaned_slices)==0){
		print('NO SLICES with variance')
		next
	    }	
            # formulate model
            formula_str <- paste(y, paste(c(cleaned_covar_features,cleaned_slices), collapse=" + "), sep=" ~ ")
    	    mymodel <- as.formula(formula_str)
            print(formula_str)
      	    df_coef = c()
   	    df_prob = c()
	    df_lowhigh = c()
            # run logistic regression
            if(FALSE){#boot=="True"){
	    	for (b in 1:100){
	            df_boot <- df[sample(1:length(df))]
                    lf <- logistf(formula = mymodel,data=df_boot,pl=FALSE)
		    betas <- coef(lf)
	            # collect data
		    df_coef <- rbind(df_coef,lf$coefficients)#[nrow(df_coef)+1,] = lf$coefficients
          	    rownames(df_coef)[length(rownames(df_coef))] <- select_slice
        	    df_prob<- rbind(df_prob,lf$prob)
        	    rownames(df_prob)[length(rownames(df_prob))] <- select_slice
        	    df_lowhigh <- rbind(df_lowhigh,c(lf$ci.lower,lf$ci.upper))#[nrow(df_lowhigh)+1,] = c(lf$ci.lower,lf$ci.upper)
        	    rownames(df_lowhigh)[length(rownames(df_lowhigh))] <- select_slice       
		}
	    }else{
                lf <- logistf(formula = mymodel,data=df,pl=FALSE)
		betas <- coef(lf)
		# collect data
		df_coef <- rbind(lf$coefficients)#[nrow(df_coef)+1,] = lf$coefficients
         	rownames(df_coef)[length(rownames(df_coef))] <- select_slice
        	df_prob<- rbind(lf$prob)
        	rownames(df_prob)[length(rownames(df_prob))] <- select_slice
        	df_lowhigh <- rbind(c(lf$ci.lower,lf$ci.upper))#[nrow(df_lowhigh)+1,] = c(lf$ci.lower,lf$ci.upper)
        	rownames(df_lowhigh)[length(rownames(df_lowhigh))] <- select_slice       
	    }
	    # save file    	    
    	    write.csv(df_coef, file =  coef_file)
            print(prob_file)
            write.csv(df_prob, file = prob_file)
        }
   }
 }
}
