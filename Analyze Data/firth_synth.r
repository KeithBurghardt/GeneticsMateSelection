library(logistf)

# number of married in HRS
num_married = 1179
# number of random in HRS
num_rand = num_married*num_married

 
# number of links
n <- num_rand + num_married

# marriage probability
y_prob = num_married/n


# p2 is the number of "other" clusters, not counting the first 2
for (p2 in c(10/n,30/n,100/n,500/n,1000/n,5000/n,10000/n)){
        print(p2)
	for (p1 in c(0.01)){
	    coefs = c()
    	    for (i in 1:10){
		set.seed(123+i)
	        x1 = rbinom(n, 1, p1)
        	x3 = rbinom(n, 1, p2)
	        # x1 and x2 are ALMOST but not quite colinear
        	x2 = 1-x1-x3
        	# this ensures x2 is 0 or 1
	        x2[x2 < 0] <- 1
		data <- data.frame(
			outcome = rbinom(n, 1, y_prob),
			x1 =x1,
			x2 = x2
			)
        # logistic regression of nearly colinear data
		model <- logistf(
			formula = outcome ~ x1 + x2,
			data = data
		)
		coefs <- rbind(coefs,coef(model))
	    }
	    write.csv(coefs, file = paste(c("synth_data_",p2,"_",p1,".csv"),collapse = "_"), row.names = FALSE)
	    print(coefs)
	}
}
