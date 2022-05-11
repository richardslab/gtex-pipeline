version 1.0 

# RNASeqQC aggregate changes the names of samples by splitting on '.' and taking the first part. 
# This workflow does the same thing to the input string

task SplitOnPeriodAndTakeFirstPart{
	input {
		String string
	}
	command<<<
	python <<EOF 
	
	string="~{string}"
	print(string.split('.')[0])
	
	EOF
	>>>
	output {
		String part=read_string(stdout())
	}
	runtime {
		docker: "python:latest"
		memory: "2GB"
		disks: "local-disk 15 HDD"
	}
}


workflow CreateRNASeqQCSampleName{
	input {
		String string
	}

	call SplitOnPeriodAndTakeFirstPart{
		input:
			string=string
	}

	output {
		String part=SplitOnPeriodAndTakeFirstPart.part
	}	
}


