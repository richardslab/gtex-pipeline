

task prepare_covariates {
	input {
		File covariate_file
		Array[String] covariate_list
		Array[String] individuals_list
		String output_file
	}

	File covariate_list_file = write_lines(covariate_list)
	File individual_list_file = write_lines(individuals_list)

	command <<<
	Rscript <<EOF "~{covariate_file}" "~{covariate_list_file}" "~{individual_list_file}" "~{output_file}"


	args <- commandArgs(trailingOnly = TRUE)

	print("got the following arguments:")
	print("args[1] (covariate file):"+args[1])

	covariate_list=read.delim(args[2])[[1]]
	print("args[2] (covariate list):"+covariate_list)
	
	individual_list=read.delim(args[3])[[1]]
	
	print("args[3] (participants list):"+individual_list)
	print("args[4] (output filename):"+args[4])

	
	covariates <- read.delim(args[1],header = TRUE,sep = '\t',row.names = FALSE)
	covariates <- subset(covariates,select=c(BQCID,covariate_list))
	covariates <- subset(covariates,BQCID %in% individual_list)

	# make sure that all requested samples are present:

	if(nrow(covariates) != length(individual_list)) {
		stop("got different number of individuals than requested.")
	}

	#prepare for transpose
	rownames(covariates) <- covariates$BQCID 
	covariates <- select(covariates,-BQCID)
	rotated <- t(covariates)
	rotated <- as.data.frame(rotated)
	rotated$ID <- rownames(rotated)

	write.table(x = rotated,file=args[4], sep="\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

	EOF

	>>>
	output {
		File prepared_covariates = output_file
	}

	runtime {
		docker: "r-base:latest"
		memory: "2GB"
		disks: "local-disk ~{disk_space} HDD"
	}

}


workflow PrepareCovariates{
	call prepare_covariates{}

	output {
		File covariates = prepare_covariates.prepared_covariates
	}
}
