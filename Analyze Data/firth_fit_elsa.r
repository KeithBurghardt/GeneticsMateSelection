library(dplyr)
library(base)
library(fastDummies)
library(logistf)
library(data.table)
library(caret)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0){
	select_slice = args[2]
	train80 = args[4]
        one_ethnic = "False"
	fold=-1
	if (length(args) > 4){
	    fold = args[6]
            one_ethnic = args[8]
	}
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
dir <- "/project/arpawong_181/HRS_AsMa/keith/ELSA/"
print("loading data")
train_outcome <- read.csv(paste(dir,"train_outcome_all_data_fold=-1_heldout.csv",sep=""))
train_X <- read.csv(paste(dir,"train_X_non_slice_features_all_data_fold=-1_heldout.csv",sep=""))
print("data loaded")

library(stringr)

# Get list of matching files
file_list <- list.files("/project/arpawong_181/burghard_687/genetic_data/slices_combined/",
                        pattern = "^elsa_slice_.*\\.csv$", full.names = TRUE)

# Extract and transform filenames
slice_features <- sapply(file_list, function(sf) {
  filename <- basename(sf)
  chr_id <- str_replace_all(filename, c("^elsa_slice_" = "", "\\.csv$" = ""))
  paste0("chr", chr_id)
})
slice_features <- sample(slice_features)

for (select_slice in slice_features){
for (abs_features in c(T)){
    for (s in c("full")){ # ,"random","nom")){
        all_files_found <- TRUE
	all_covars <- c("place","ethnicity","relig","edu","height","cog27")
        covar_combos <- c("none",paste(all_covars,collapse="+"))
        out_dir <- "/project/arpawong_181/burghard_687/genetic_data/"
        for (covar in covar_combos){
            print(covar)
            # file names
            folder = ""
            if (abs_features){
                folder = "_abs"
            }else{
                folder = "_nonabs"
            }
            
            #coef_file = paste(c( out_dir,"slices_all",folder,one_eth,"_heldout_",fold,"_new2/","firth_",covar,"_coefs_",select_slice,"_",s,folder,"_no_cog-bmi_cog27_fold=",fold,"_heldout_wald_newslice.csv"),collapse="")
            coef_file = paste(c( out_dir,"slices_all",folder,one_eth,"_heldout_",fold,"_new2/","firth_",covar,"_coefs_",select_slice,"_",s,folder,"_no_cog-bmi_cog27_fold=",fold,"_heldout_wald-multicluster.csv"),collapse="")#,"-EV+TestSep.csv"),collapse="")
            print(coef_file)
            #prob_file = paste(c( out_dir,"slices_all",folder,one_eth,"_heldout_",fold,"_new2/","firth_",covar,"_prob_",select_slice,"_",s,folder,"_no_cog-bmi_cog27_fold=",fold,"_heldout_wald_newslice.csv"),collapse="")
            prob_file = paste(c( out_dir,"slices_all",folder,one_eth,"_heldout_",fold,"_new2/","firth_",covar,"_prob_",select_slice,"_",s,folder,"_no_cog-bmi_cog27_fold=",fold,"_heldout_wald-multicluster.csv"),collapse="")#"-EV+TestSep.csv"),collapse="")

            coef_file_alt = paste(c( out_dir,"slices_all",folder,one_eth,"_heldout_",2,"_new2/","firth_",covar,"_coefs_",select_slice,"_",s,folder,"_no_cog-bmi_cog27_fold=",fold,"_heldout_wald.csv"),collapse="")
            prob_file_alt = paste(c( out_dir,"slices_all",folder,one_eth,"_heldout_",2,"_new2/","firth_",covar,"_prob_",select_slice,"_",s,folder,"_no_cog-bmi_cog27_fold=",fold,"_heldout_wald.csv"),collapse="")
            # ignore analysis if we have already calculated values
            if ((file.exists(coef_file) && file.exists(prob_file)) || (file.exists(coef_file_alt) && file.exists(prob_file_alt))){# && file.exists(lowhigh_file) && file.exists(ll_file) && file.exists(r2_file)){
                print('old file')
                next
            }else{
		# if we have not found all files, continue (skipping now is faster than waiting to load data first)
		all_files_found <- FALSE
		}
        }
        #print(all_files_found)
	if (all_files_found){
	    next
	}

        print('loading train slice')
        ######## CLUSTER ##########
        slice_file <- paste("/project/arpawong_181/burghard_687/genetic_data/slices_combined/elsa_slice_",gsub("chr", "",select_slice),".csv",sep="")
        #print(slice_file)
        cluster_ids <-read.csv(slice_file)
        print(max(cluster_ids))
	if(max(cluster_ids) == -1){next}
        # cut off dummy variables
        cluster_ids <- data.frame(cluster_ids)
        colnames(cluster_ids) <- c("newclust")
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
        ###############################

        # this ensures all rows are valid

        df <-df[complete.cases(df), ]

        # abs
        if (abs_features){
		df <- abs(df)
	}
        #print('df finalized')

	# outcome variable
	y <- "marriage"
	# create name for slice
        slices <- colnames(df)[grepl("newclust", colnames(df), fixed = TRUE)]
	# name of non-slice features
        covars_list <- colnames(df)[!grepl("newclust", colnames(df), fixed = TRUE)]
	# all relevant feature combinations
        out_dir <- '/project/arpawong_181/burghard_687/genetic_data/'
	for (covar in covar_combos){
            print(covar)
            # file names
	    folder = ""
	    if (abs_features){
		folder = "_abs"
	    }else{
		folder = "_nonabs"
	    }
            coef_file = paste(c( out_dir,"slices_all",folder,one_eth,"_heldout_",fold,"_new2/","firth_",covar,"_coefs_",select_slice,"_",s,folder,"_no_cog-bmi_cog27_fold=",fold,"_heldout_wald-multicluster.csv"),collapse="")#,"-EV+TestSep.csv"),collapse="")
            ll_file = paste(c( out_dir,"slices_all",folder,one_eth,"_heldout_",fold,"_new2/","firth_",covar,"_loglik_",select_slice,"_",s,folder,"_no_cog-bmi_cog27_fold=",fold,"_heldout_wald-multicluster.csv"),collapse="")#,"-EV+TestSep.csv"),collapse="")
            r2_file = paste(c( out_dir,"slices_all",folder,one_eth,"_heldout_",fold,"_new2/","firth_",covar,"_r2_",select_slice,"_",s,folder,"_no_cog-bmi_cog27_fold=",fold,"_heldout_wald-multicluster.csv"),collapse="")#,"-EV+TestSep.csv"),collapse="")
            prob_file = paste(c( out_dir,"slices_all",folder,one_eth,"_heldout_",fold,"_new2/","firth_",covar,"_prob_",select_slice,"_",s,folder,"_no_cog-bmi_cog27_fold=",fold,"_heldout_wald-multicluster.csv"),collapse="")#"-EV+TestSep.csv"),collapse="")
            lowhigh_file = paste(c( out_dir,"slices_all",folder,one_eth,"_heldout_",fold,"_new2/","firth_",covar,"_lowhigh_",select_slice,"_",s,folder,"_no_cog-bmi_cog27_fold=",fold,"_heldout_wald-multicluster.csv"),collapse="")#,"-EV+TestSep.csv"),collapse="")
            

            coef_file_alt = paste(c( out_dir,"slices_all",folder,one_eth,"_heldout_",2,"_new2/","firth_",covar,"_coefs_",select_slice,"_",s,folder,"_no_cog-bmi_cog27_fold=",fold,"_heldout.csv"),collapse="")#,"-EV+TestSep.csv"),collapse="")
            prob_file_alt = paste(c( out_dir,"slices_all",folder,one_eth,"_heldout_",2,"_new2/","firth_",covar,"_prob_",select_slice,"_",s,folder,"_no_cog-bmi_cog27_fold=",fold,"_heldout.csv"),collapse="")#"-EV+TestSep.csv"),collapse="")
            # prob_file = dump_dir+'slices_all'+folder+u_one_eth+'_cog27_'+str(fold)+'_new/firth_'+cov+'_prob_'+slice+'_'+sample+folder+'_no_cog-bmi_cog27_fold='+str(fold)+'_heldout.csv'
	    if (train80 == "True"){
                coef_file = paste(c( dir,"slices",folder,one_eth,"_heldout_",fold,"/","firth_",covar,"_coefs_",select_slice,"_",s,"_train80_fold=",fold,folder,"_one_ethnic=",one_ethnic,"_heldout.csv"),collapse="")#,"-EV+TestSep.csv"),collapse="")
                prob_file = paste(c( dir,"slices",folder,one_eth,"_heldout_",fold,"/","firth_",covar,"_prob_",select_slice,"_",s,"_train80_fold=",fold,folder,"_one_ethnic=",one_ethnic,"_heldout.csv"),collapse="")#,"-EV+TestSep.csv"),collapse="")
                lowhigh_file = paste(c( dir,"slices",folder,one_eth,"_heldout_",fold,"/","firth_",covar,"_lowhigh_",select_slice,"_",s,"_train80_fold=",fold,folder,"_one_ethnic=",one_ethnic,"_heldout.csv"),collapse="")#,"-EV+TestSep.csv"),collapse="")
                ll_file = paste(c( dir,"slices",folder,one_eth,"_heldout_",fold,"/","firth_",covar,"_loglik_",select_slice,"_",s,"_train80_fold=",fold,folder,"_one_ethnic=",one_ethnic,"_heldout.csv"),collapse="")#,"-EV+TestSep.csv"),collapse="")
                r2_file = paste(c( dir,"slices",folder,one_eth,"_heldout_",fold,"/","firth_",covar,"_r2_",select_slice,"_",s,"_train80_fold=",fold,folder,"_one_ethnic=",one_ethnic,"_heldout.csv"),collapse="")#,"-EV+TestSep.csv"),collapse="")

                coef_file_alt = paste(c( dir,"slices",folder,one_eth,"_heldout_2/","firth_",covar,"_coefs_",select_slice,"_",s,"_train80_fold=",fold,folder,"_one_ethnic=",one_ethnic,"_heldout.csv"),collapse="")#,"-EV+TestSep.csv"),collapse="")
                prob_file_alt = paste(c( dir,"slices",folder,one_eth,"_heldout_2/","firth_",covar,"_prob_",select_slice,"_",s,"_train80_fold=",fold,folder,"_one_ethnic=",one_ethnic,"_heldout.csv"),collapse="")#,"-EV+TestSep.csv"),collapse="")
            }
            print(prob_file)
            # ignore analysis if we have already calculated values
            if ((file.exists(coef_file) && file.exists(prob_file))){ # || (file.exists(coef_file_alt) && file.exists(prob_file_alt))){# && file.exists(lowhigh_file) && file.exists(ll_file) && file.exists(r2_file)){
                print('old file')
            	next
            } 
            # collect covariate features
	    covar_features <- c()
	    if (grepl("place",covar)){
	        covar_features <- c(covar_features,c())
	    }
            if (grepl("ethnicity_only",covar)){
                covar_features <- c(covar_features,c())
            }
            if (grepl("pca_distance",covar)){
                covar_features <- c(covar_features,c('pca_distance'))
            }
	    if (grepl("edu_only",covar)){
	        covar_features <- c(covar_features,c('edyrs'))
	    }
	    if (grepl("cog27",covar)){ 
	        covar_features <- c(covar_features,c('cognition'))
	    }
	    if (grepl("income",covar)){ 
	        covar_features <- c(covar_features,c('htot'))
	    }
	    if (grepl("height",covar)){
	        covar_features <- c(covar_features,c('mheight'))
	    }
            if (grepl("ethnicity",covar) & covar != "ethnicity_only"){
		if(grepl("nopca",covar)){
   		    covar_features <- c()
		}else{
	            covar_features <- c(covar_features,c('pca_distance'))
		}
	    }

	    if (grepl("edu",covar) & covar != "edu_only"){
	        covar_features <- c(covar_features,c('edyrs','htot'))#'cogtot'
	    }
            print(covar_features)
            # create data frame
            df_coef = data.frame(rbind(c(c("intercept"),covar_features,slices)))
            colnames(df_coef) <- c(c("intercept"),covar_features,slices)
            df_prob = data.frame(rbind(c(c("intercept"),covar_features,slices)))
            colnames(df_prob) <- c(c("intercept"),covar_features,slices)
            df_lowhigh = data.frame(rbind(c(c("intercept"),covar_features,slices,c("intercept"),covar_features,slices)))
            colnames(df_lowhigh) <- c(c("intercept"),covar_features,slices,c("intercept"),covar_features,slices)


	    for (slice in slices) {
                cleaned_covar_features <- c()
                for(col in c(covar_features)){
			if(sd(df[,col])> 0.000001){
                        print(col)
                        print(sd(df[,col]))
		        cleaned_covar_features <- c(cleaned_covar_features,col)
			}
		}
            cleaned_slices <- c()
            for (col in c(slices)){
                if( (sd(df[,col])> 0.00001) & (sum(df[,col])>3)){
                print(col)
                print(sum(df[,col]))
                print(sd(df[,col]))
                cleaned_slices <- c(cleaned_slices,col)
                }
            }
            # if there are no viable slices, ignore
            if (length(cleaned_slices)==0){
                print('NO SLICES with variance')
                next
            }
            print(cleaned_covar_features)
            print(cleaned_slices)
            # formulate model
            formula_str <- paste(y, paste(c(cleaned_covar_features,cleaned_slices), collapse=" + "), sep=" ~ ")
                mymodel <- as.formula(formula_str)
                print(formula_str)
                # run logistic regression
                lf <- logistf(formula = mymodel,data=df,pl=FALSE)

		betas <- coef(lf)
                
		X <- model.matrix(mymodel, data=df)
		pred <- 1 / (1 + exp(-X * betas)) 
                cond_1 = cond <- sapply(train_outcome, function(x) x == 1)
                cond_0 = cond <- sapply(train_outcome, function(x) x == 0)
                r2 = data.table(sum(pred[cond_1])/length(pred[cond_1]) - sum(pred[cond_0])/length(pred[cond_0]))

		ll = lf$loglik
                
		# collect data
	        df_coef <- rbind(lf$coefficients)
                rownames(df_coef)[length(rownames(df_coef))] <- select_slice
        	df_prob<- rbind(lf$prob)
        	rownames(df_prob)[length(rownames(df_prob))] <- select_slice
        	df_lowhigh <- rbind(c(lf$ci.lower,lf$ci.upper))#[nrow(df_lowhigh)+1,] = c(lf$ci.lower,lf$ci.upper)
        	rownames(df_lowhigh)[length(rownames(df_lowhigh))] <- select_slice       
		# save file    	    
    	        write.csv(df_coef, file =  coef_file)
                print(prob_file)
        	write.csv(df_prob, file = prob_file)
                q()
            }
	}
    }
}
#}
}
