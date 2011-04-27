#####
## 2010-4-6 R code to deal with Justin borevitz lab's season flowering data. Postdoc Yan Li give me some RData. 
## this file has function to output its RData and run EMMAX.
##
load("~/script/variation/data/JBLabSeasonFlowering/d4DTF.RData")
m <- read.csv("http://borevitzlab.uchicago.edu/Members/yanli1/250kdata/pred14snpsEig10res_2.csv")
load("~/script/variation/data/JBLabSeasonFlowering/core469_emma_input.RData")

objects();

# d4 is the dataset with phenotype, ecotypeid, planing, loc, shelfHeight, SNPs, PCs. with replicates. around 469*16 rows.
# xs is a SNP Data matrix. rows are SNPs. columns are ecotypes.  213637 X 469
# K is the kinship matrix 469X469.
# m is a data frame storing names & effects of 14 SNPs + 3 interactions

createIndividualToLineIncidenceMatrix = function(lineEcotypeID_vector, individualEcotypeIDvector){
	# ecotype id in the latter is of type integer. the former is of type character.
	no_of_lines = length(lineEcotypeID_vector)
	no_of_replicates = length(individualEcotypeIDvector)
	
	Z = matrix(nrow=no_of_replicates, ncol=no_of_lines)	# incidence matrix
	Z[,] = 0
	
	# turn lineEcotypeID_vector into a matrix of integers
	lineEcotypeID_matrix = matrix(nrow=no_of_lines)
	for( i in 1:(no_of_lines) )
	{
		lineEcotypeID_matrix[i] = as.integer(lineEcotypeID_vector[i])
	}
	#cat("lineEcotypeID_matrix: ", lineEcotypeID_matrix, "\n")
	#cat("length of Z[1][1] is ", length(Z[1][1]), "\n")
	
	for ( i in 1:no_of_replicates )
	{
		individualEcotypeID = as.integer(individualEcotypeIDvector[i])
		for( j in 1:no_of_lines )
		{
			if (individualEcotypeID==lineEcotypeID_matrix[j])
			{
				Z[i,j] <- 1		#Z[i][j] gives warning "number of items to replace is not a multiple of replacement length"
				break
			}
		}
	}
	return (Z)
}

createDesignMatrix = function(d4, m){
	no_of_individuals = dim(d4)[1]
	no_of_snps = dim(m)[1]
	design_matrix = matrix(nrow=no_of_individuals, ncol=no_of_snps+4)
	design_matrix[,] = 0
	design_matrix_in_list = list()
	no_of_singleton_snps = no_of_snps-3 
	for ( i in 1:no_of_singleton_snps)
	{
		snp_id = m$SNP[i]
		design_matrix_in_list[[snp_id]] = d4[as.character(snp_id)]
		design_matrix[,i] = unlist(d4[as.character(snp_id)], use.names=FALSE)
		# d4$snp_id doesn't work. even though d4$"X5_3188328" or d4$X5_3188328 does.
		# design_matrix[,i] is not a 2D matrix, just a vector.
		# d4[as.character(snp_id)] is still a data frame. unlist would cast it into a numeric vector.
		# d4$X5_3188328 is a numeric vector by itself.
	}
	# 3 SNP interactions
	design_matrix[, no_of_singleton_snps+1] = d4$X5_3188328 * d4$X4_264496
	design_matrix[, no_of_singleton_snps+2] = d4$X5_3188328 * d4$X1_24345319
	design_matrix[, no_of_singleton_snps+3] = d4$X1_24345319 * d4$X4_203180
	
	design_matrix_in_list$planting = d4$planting
	
	for (i in 1:no_of_individuals)
	{
		if (d4$planting[i]=="spring")
		{
			design_matrix[i,no_of_singleton_snps+4] = 0
		}
		else
		{
			design_matrix[i,no_of_singleton_snps+4] = 1
		}
	}
	
	design_matrix_in_list$loc = d4$loc
	for (i in 1:no_of_individuals)
	{
		if (d4$loc[i]=="spain")
		{
			design_matrix[i,no_of_singleton_snps+5] = 0
		}
		else
		{
			design_matrix[i,no_of_singleton_snps+5] = 1
		}
	}
	
	design_matrix_in_list$shelfHeight = d4$shelfHeight
	for (i in 1:no_of_individuals)
	{
		if (d4$shelfHeight[i]=="bottom")
		{
			design_matrix[i,no_of_singleton_snps+6] = 0	# top or not
			design_matrix[i,no_of_singleton_snps+7] = 0	# middle or not
		}
		else if (d4$shelfHeight[i]=="bottom-middle")
		{
			design_matrix[i,no_of_singleton_snps+6] = 0
			design_matrix[i,no_of_singleton_snps+7] = 1
		}
		else if (d4$shelfHeight[i]=="top")
		{
			design_matrix[i,no_of_singleton_snps+6] = 1
			design_matrix[i,no_of_singleton_snps+7] = 0
		}
		else if (d4$shelfHeight[i]=="top-middle")
		{
			design_matrix[i,no_of_singleton_snps+6] = 1
			design_matrix[i,no_of_singleton_snps+7] = 1
		}
	}
	return (list(design_matrix=design_matrix, design_matrix_in_list=design_matrix_in_list))
}

# 2010-3-22 output the d4 (selected SNP + phenotype + environment data) + K

write.table(d4, file='d4DTF.tsv', append = FALSE, quote = FALSE, sep = "\t",
	row.names = FALSE, col.names = TRUE, qmethod = c("escape", "double"))
write.table(K, file='K.tsv', append = FALSE, quote = FALSE, sep = "\t",
	row.names =attributes(xs)$dimnames[2][[1]] , col.names = attributes(xs)$dimnames[2][[1]],
	qmethod = c("escape", "double"))

Z = createIndividualToLineIncidenceMatrix(attributes(xs)$dimnames[2][[1]], d4$ecotypeid)

new_K = Z%*%K%*%t(Z)

no_of_replicates = length(d4$ecotypeid)

X0 <- matrix(1, no_of_replicates, 1)	# intercept

source('~/script/variation/src/gwa/emma/R/emma.R')
one_marker_rs <- emma.REMLE(d4$DTF, X0, new_K, cal.pvalue=TRUE)

cat("pvalue=", one_marker_rs$pvalue, " stat=", one_marker_rs$stat,  " vg=", one_marker_rs$vg, 
		" ve=", one_marker_rs$ve,
		" pvalues=", one_marker_rs$pvalues,
		" beta=", one_marker_rs$beta,
		" x_beta_var=", one_marker_rs$x_beta_var,
		" REML=maxLL=", one_marker_rs$REML,
		" delta=", one_marker_rs$delta,
		" genotype_var_perc=", one_marker_rs$genotype_var_perc,
		"\n")

design_matrix_return_data = createDesignMatrix(d4, m)
design_matrix = design_matrix_return_data$design_matrix

design_matrix_in_list = design_matrix_return_data$design_matrix_in_list

design_matrix_frame = as.data.frame(design_matrix_in_list)


design_matrix_frame$ecotypeid = d4$ecotypeid
design_matrix_frame$DTF = d4$DTF

EMMAX <- function()
{	
	source('~/script/variation/src/gwa/emma/R/emma.R')
	one_marker_rs <- emma.REMLE(d4$DTF, X0, new_K, cal.pvalue=TRUE)
	
	cat("pvalue=", one_marker_rs$pvalue, " stat=", one_marker_rs$stat,  " vg=", one_marker_rs$vg, 
			" ve=", one_marker_rs$ve,
			" pvalues=", one_marker_rs$pvalues,
			" beta=", one_marker_rs$beta,
			" x_beta_var=", one_marker_rs$x_beta_var,
			" REML=maxLL=", one_marker_rs$REML,
			" delta=", one_marker_rs$delta,
			" genotype_var_perc=", one_marker_rs$genotype_var_perc,
			"\n")
	# generalized least square
	
	# generate a variance matrix
	I <- matrix(0,nrow=no_of_replicates,ncol=no_of_replicates)
	I[row(I)==col(I)] <- 1
	variance_matrix = one_marker_rs$vg*new_K + one_marker_rs$ve*I
	# cholesky decomposition of variance matrix
	L = t(chol(variance_matrix))
	L_inverse = solve(L)
	
	
	transformed_phenotype_vector = L_inverse %*% as.matrix(d4$DTF, no_of_replicates,1)
	transformed_design_matrix = L_inverse %*% design_matrix
	transformed_data_frame = as.data.frame(cbind(transformed_phenotype_vector, transformed_design_matrix))
	
	lm_result = lm(V1~V2+V3+V4+V5+V6+V7+V8+V9+V10+V11+V12+V13+V14+V15+V16+V17+V18+V19+V20+V21+V22, data=transformed_data_frame)
	cat(summary(lm_result))
	cat(anova(lm_result))

}


# pure EMMA
X <- cbind(X0, design_matrix)
source('~/script/variation/src/gwa/emma/R/emma.R')
one_marker_rs <- emma.REMLE(d4$DTF, X, new_K, cal.pvalue=TRUE)

cat("pvalue=", one_marker_rs$pvalue, " stat=", one_marker_rs$stat,  " vg=", one_marker_rs$vg, 
		" ve=", one_marker_rs$ve,
		" pvalues=", one_marker_rs$pvalues,
		" beta=", one_marker_rs$beta,
		" x_beta_var=", one_marker_rs$x_beta_var,
		" REML=maxLL=", one_marker_rs$REML,
		" delta=", one_marker_rs$delta,
		" genotype_var_perc=", one_marker_rs$genotype_var_perc,
		"\n")