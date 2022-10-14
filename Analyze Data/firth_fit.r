# requires logistf, available at https://cran.r-project.org/web/packages/logistf/index.html
library(logistf)
library(data.table)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0){
	select_slice = args[2]
	train80 = args[4]
	fold=-1
	if (length(args) > 4){
	    fold = args[6]
	}
}else{
	select_slice = "None"
	train80 = "True"
	fold = 0
}
for (abs_features in c(T)){
    for (sample in c("full")){
	# directory to collect/store data
	dir <- "/project/burghard_687/genetic_data/"
        if (select_slice != "None"  & train80!="True"){
	    train_outcome <- read.csv(paste(dir,"train_outcome_cog27.csv",sep=""))
   	    train_X <- read.csv(paste(dir,"train_X_non_slice_features_cog27",".csv",sep=""))

            train_slice <-read.csv(paste(dir,"slices/train_",select_slice,"_cog27.csv",sep=""))}
	if (train80=="True"){
                train_outcome <- read.csv(paste(dir,"train_outcome_train80_fold=",fold,"_cog27.csv",sep=""))
                train_X <- read.csv(paste(dir,"train_X_non_slice_features_train80_fold=",fold,"_cog27.csv",sep=""))
                if (select_slice != "None"){
                    train_slice <-read.csv(paste(dir,"slices80/train_",select_slice,"_train80_fold=",fold,"_cog27.csv",sep=""))

		}
        }
        if (select_slice == "None"){
            # make feature 0s (effectively remove parameter)
            train_slice <- as.data.frame(matrix(0,length(train_X),1))
            colnames(train_slice) <- c(c('slice'))
	}
        train_X <- cbind(train_X,train_slice)
	num_marriages = length(subset(train_outcome, marriage == 1)$marriage)
	# combine all training features
	df <- cbind(train_X,train_outcome)
        # abs
        if (abs_features){
		df <- abs(df)
	}
	# outcome variable
	y <- "marriage"
	# create name for slice
	slices <- colnames(df)[grepl("slice", colnames(df), fixed = TRUE)]
	# name of non-slice features
	covars_list <- colnames(df)[!grepl("slice", colnames(df), fixed = TRUE)]
	all_covars <- c("place","ethnicity","relig","edu","height","cog27")#"bmi"

    covar_combos <- c("ethnicity_only","pca_distance","edu_only","income","num_marriage","num_div","cog27","relig","height","none","place",paste(all_covars,collapse="+"), paste(c(paste(all_covars,collapse="+"),"nopca"),collapse="+"))##"bmi","cog""place+ethnicity","place+ethnicity+relig",
	for (covar in covar_combos){
            print(covar)
            # file names
	    folder = ""
	    if (abs_features){
		folder = "_abs"
	    }else{
		folder = "_nonabs"
	    }
            coef_file = paste(c( dir,"slices",folder,"_cog27/","firth_",covar,"_coefs_",select_slice,"_",sample,folder,"_no_cog-bmi_cog27.csv"),collapse="")#,"-EV+TestSep.csv"),collapse="")
            ll_file = paste(c( dir,"slices",folder,"_cog27/","firth_",covar,"_loglik_",select_slice,"_",sample,folder,"_no_cog-bmi_cog27.csv"),collapse="")#,"-EV+TestSep.csv"),collapse="")
            r2_file = paste(c( dir,"slices",folder,"_cog27/","firth_",covar,"_r2_",select_slice,"_",sample,folder,"_no_cog-bmi_cog27.csv"),collapse="")#,"-EV+TestSep.csv"),collapse="")
            prob_file = paste(c( dir,"slices",folder,"_cog27/","firth_",covar,"_prob_",select_slice,"_",sample,folder,"_no_cog-bmi_cog27.csv"),collapse="")#"-EV+TestSep.csv"),collapse="")
            lowhigh_file = paste(c( dir,"slices",folder,"_cog27/","firth_",covar,"_lowhigh_",select_slice,"_",sample,folder,"_no_cog-bmi_cog27.csv"),collapse="")#,"-EV+TestSep.csv"),collapse="")

            coef_file_alt = paste(c( dir,"slices",folder,"_cog27_2/","firth_",covar,"_coefs_",select_slice,"_",sample,folder,"_no_cog-bmi_cog27.csv"),collapse="")#,"-EV+TestSep.csv"),collapse="")
            prob_file_alt = paste(c( dir,"slices",folder,"_cog27_2/","firth_",covar,"_prob_",select_slice,"_",sample,folder,"_no_cog-bmi_cog27.csv"),collapse="")#"-EV+TestSep.csv"),collapse="")
	    if (train80 == "True"){
                coef_file = paste(c( dir,"slices",folder,"_cog27/","firth_",covar,"_coefs_",select_slice,"_",sample,"_train80_fold=",fold,folder,"_no_cog-bmi_cog27.csv"),collapse="")#,"-EV+TestSep.csv"),collapse="")
                prob_file = paste(c( dir,"slices",folder,"_cog27/","firth_",covar,"_prob_",select_slice,"_",sample,"_train80_fold=",fold,folder,"_no_cog-bmi_cog27.csv"),collapse="")#,"-EV+TestSep.csv"),collapse="")
                lowhigh_file = paste(c( dir,"slices",folder,"_cog27/","firth_",covar,"_lowhigh_",select_slice,"_",sample,"_train80_fold=",fold,folder,"_no_cog-bmi_cog27.csv"),collapse="")#,"-EV+TestSep.csv"),collapse="")
                ll_file = paste(c( dir,"slices",folder,"_cog27/","firth_",covar,"_loglik_",select_slice,"_",sample,"_train80_fold=",fold,folder,"_no_cog-bmi_cog27.csv"),collapse="")#,"-EV+TestSep.csv"),collapse="")
                r2_file = paste(c( dir,"slices",folder,"_cog27/","firth_",covar,"_r2_",select_slice,"_",sample,"_train80_fold=",fold,folder,"_no_cog-bmi_cog27.csv"),collapse="")#,"-EV+TestSep.csv"),collapse="")

                coef_file_alt = paste(c( dir,"slices",folder,"_cog27_2/","firth_",covar,"_coefs_",select_slice,"_",sample,"_train80_fold=",fold,folder,"_no_cog-bmi_cog27.csv"),collapse="")#,"-EV+TestSep.csv"),collapse="")
                prob_file_alt = paste(c( dir,"slices",folder,"_cog27_2/","firth_",covar,"_prob_",select_slice,"_",sample,"_train80_fold=",fold,folder,"_no_cog-bmi_cog27.csv"),collapse="")#,"-EV+TestSep.csv"),collapse="")
            }

            if ((file.exists(coef_file) && file.exists(prob_file)) || (file.exists(coef_file_alt) && file.exists(prob_file_alt))){
                print('old file')
            	next
            } 
            # collect covariate features
	    covar_features <- c()
	    if (grepl("place",covar)){
	        covar_features <- c(covar_features,c('RABPLACE_1.0','RABPLACE_2.0','RABPLACE_3.0','RABPLACE_4.0','RABPLACE_5.0','RABPLACE_6.0','RABPLACE_7.0','RABPLACE_8.0','RABPLACE_9.0','RABPLACE_10.0','RABPLACE_11.0'))
	    }
            if (grepl("ethnicity_only",covar)){
                covar_features <- c(covar_features,c('ethnicity_1.0','ethnicity_2.0','ethnicity_3.0','ethnicity_4.0'))
            }
            if (grepl("pca_distance",covar)){
                covar_features <- c(covar_features,c('pca_distance'))
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
	    if (grepl("height",covar)){
	        covar_features <- c(covar_features,c('height'))
	    }
            if (grepl("ethnicity",covar) & covar != "ethnicity_only"){
		if(grepl("nopca",covar)){
   		    covar_features <- c(covar_features,c('ethnicity_1.0','ethnicity_2.0','ethnicity_3.0','ethnicity_4.0'))
		}else{
	            covar_features <- c(covar_features,c('ethnicity_1.0','ethnicity_2.0','ethnicity_3.0','ethnicity_4.0','pca_distance'))
		}

	    }

	    if (grepl("edu",covar) & covar != "edu_only"){
	        covar_features <- c(covar_features,c('RAEDYRS','iearn','mdiv','mrct'))#'cogtot'
	    }

            # create data frame
	    df_coef = data.frame(rbind(c(c("intercept"),covar_features,c("slice"))))
	    colnames(df_coef) <- c(c("intercept"),covar_features,c("slice"))
	    df_prob = data.frame(rbind(c(c("intercept"),covar_features,c("slice"))))
	    colnames(df_prob) <- c(c("intercept"),covar_features,c("slice"))
	    df_lowhigh = data.frame(rbind(c(c("intercept"),covar_features,c("slice"),c("intercept"),covar_features,c("slice"))))
	    colnames(df_lowhigh) <- c(c("intercept"),covar_features,c("slice"),c("intercept"),covar_features,c("slice"))
	    for (slice in slices) {
                cleaned_covar_features <- c()
                for(col in c(covar_features)){
			if(sd(df[,col])> 0){
				cleaned_covar_features <- c(cleaned_covar_features,col)
			}
		}
                
                # formulate model
                formula_str <- paste(y, paste(c(cleaned_covar_features,slice), collapse=" + "), sep=" ~ ")
        	mymodel <- as.formula(formula_str)
                print(formula_str)
                # run logistic regression
                lf <- logistf(formula = mymodel,data=df)

		betas <- coef(lf)
                
		X <- model.matrix(mymodel, data=df)
		pred <- 1 / (1 + exp(-X * betas)) 
                cond_1 = cond <- sapply(train_outcome, function(x) x == 1)
                cond_0 = cond <- sapply(train_outcome, function(x) x == 0)
                r2 = data.table(sum(pred[cond_1])/length(pred[cond_1]) - sum(pred[cond_0])/length(pred[cond_0]))

		ll = lf$loglik
                
		# collect data
	        df_coef <- rbind(lf$coefficients)#[nrow(df_coef)+1,] = lf$coefficients
        	rownames(df_coef)[length(rownames(df_coef))] <- select_slice
        	df_prob<- rbind(lf$prob)
        	rownames(df_prob)[length(rownames(df_prob))] <- select_slice
        	df_lowhigh <- rbind(c(lf$ci.lower,lf$ci.upper))#[nrow(df_lowhigh)+1,] = c(lf$ci.lower,lf$ci.upper)
        	rownames(df_lowhigh)[length(rownames(df_lowhigh))] <- select_slice       
		# save file    	    
    	        write.csv(df_coef, file =  coef_file)
        	write.csv(df_prob, file = prob_file)
            }
        }
    }
}
