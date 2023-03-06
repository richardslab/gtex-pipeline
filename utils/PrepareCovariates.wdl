version 1.0

task prepare_covariates {
	input {
		File covariate_file
		Array[String] covariate_list
		Array[String] identifier_list
		String identifier_column="BQCID" #column header in the covariate file that will be matched against the identifier list
		Array[String]? individual_list
	}

	File covariate_list_file = write_lines(covariate_list)
	File identifier_list_file = write_lines(identifier_list)
	File individual_list_file = write_lines(select_first([individual_list,[]]))
	
	command <<<
	# need to quote the heredoc word so that the '$' in the R code do not get misunderstood.

	R --vanilla --no-save --file=-<<"EOF" --args "~{covariate_file}" "~{covariate_list_file}" "~{identifier_list_file}" "~{individual_list_file}" 

	args <- commandArgs(trailingOnly = TRUE)

	covariate_file <- args[1]
	cat(sprintf("got the following arguments:\n"))
	cat(sprintf("args[1] (covariate file): %s\n", covariate_file))

	covariate_list <- read.delim(args[2], header=FALSE)[[1]]
	cat(sprintf("args[2] (covariate list): %s\n", paste(collapse=", ", covariate_list)))
	
	identifier_list <- read.delim(args[3], header=FALSE)[[1]]
	cat(sprintf("args[3] (identifier list): %s\n", paste(collapse=", ", identifier_list)))

	individual_list <- read.delim(args[4], header=FALSE)[[1]]
	cat(sprintf("args[3] (individual list): %s\n", paste(collapse=", ", individual_list)))
	
	if (length(individual_list)==0) {
		cat("no individual list provided, using identifier as individual.")
		individual_list <- identifier_list
	}

	if (length(individual_list)!=length(identifier_list)) {
		cat(sprintf("lengths of individual and identifier lists mismatch: %d and %d", length(individual_list),length(identifier_list)))
	}

	ids=data.frame(identifier=identifier_list,individual=individual_list)

	covariates <- read.delim(covariate_file, header = TRUE, sep = '\t')

	cat(names(covariates))
	cat("\n")

	unavailable_columns <- setdiff(covariate_list, names(covariates))

	if (length(unavailable_columns)!=0) {
		sprintf("some requested covariates are unavailable: [%s]", paste(collapse=", ", unavailable_columns))
	}	

	covariates <- subset(merge(y=covariates,x=ids, by.y="~{identifier_column}",by.x="identifier",all.x=TRUE),select=c("individual","identifier", covariate_list))

	show(covariates)

	# make sure that all requested samples are present:
	unavailable_identifiers <- setdiff(identifier_list, covariates$identifier)
	if (length(unavailable_identifiers)>0) {

			stop(sprintf("Got different number of individuals than requested: %d vs. %d.\n Samples requested that were not returned are: %s\n", 
					nrow(covariates), length(individual_list), paste(collapse=", ",unavailable_identifiers)))
		}

	if (any(duplicated(covariates$identifier))){
		cat( covariates[which(duplicated(covariates$identifier),"identifier")])
	
		stop(sprintf("got some row ids more than once: %d", 
			paste(collapse=", ", covariates[which(duplicated(covariates$identifier)),"identifier"])))
	
	}

	if (any(duplicated(covariates$individual))){
		cat( covariates[which(duplicated(covariates$individual),"individual")])
	
		stop(sprintf("got some individuals more than once: %d", 
			paste(collapse=", ", covariates[which(duplicated(covariates$individual)),"individual"])))
	
	}


	# prepare for transpose
	rownames(covariates) <- covariates$individual 
	covariates <- subset(covariates, select=c(-individual,-identifier))

	rotated <- t(covariates)
	rotated <- as.data.frame(rotated)
	rotated$ID <- rownames(rotated)
	# move ID to the first column
	rotated <- rotated[,c(ncol(rotated), seq(ncol(rotated) - 1))]
	write.table(x = rotated, file="extracted_covariates.tsv", sep="\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

	EOF

	>>>
	output {
		File prepared_covariates = "extracted_covariates.tsv"
	}

	runtime {
		docker: "r-base:latest"
		memory: "2GB"
		disks: "local-disk 20 HDD"
	}

}


workflow PrepareCovariates{
	call prepare_covariates{}

	output {
		File covariates = prepare_covariates.prepared_covariates
	}
}
