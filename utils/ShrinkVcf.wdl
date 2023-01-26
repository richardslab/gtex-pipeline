version 1.0
## import "write_array.wdl" as wa
## import "Error.wdl" as e

task ShrinkVcfTask {
	input {
		File vcf
		File vcf_index

		Int threads=1
		String basename
	}
	command <<<
		wget https://raw.githubusercontent.com/broadinstitute/palantir-workflows/main/Scripts/monitoring/cromwell_monitoring_script.sh
		bash ./cromwell_monitoring_script.sh | tee monitoring.log &

		set -euo pipefail

		#remove format field (except genotype), and remove any site that is filtered
		bcftools annotate  -x FORMAT -Oz -o ~{basename}.vcf.gz ~{vcf} 

		bcftools index -t ~{basename}.vcf.gz 

	>>>
	output {
		File vcf_out="~{basename}.vcf.gz"
		File vcf_out_index="~{basename}.vcf.gz.tbi"
		File monitoring_log="monitoring.log"
	}

	runtime {
		docker: "bschiffthaler/bcftools:latest"
		preemptible: 0
		disks: "local-disk " + (2*ceil(size(vcf,"GiB"))+20) + " HDD"
		bootDiskSizeGb: "16"
		memory: 20 + " GB"
		cpu: "~{threads}"
	}
}


workflow ShrinkVcfs{
	input {
		File vcf_in
		File vcf_in_index

		Int threads=1
		String basename
	}


	call ShrinkVcfTask {
		input:
			vcf=vcf_in,
			vcf_index=vcf_in_index,
			threads=threads,
			basename=basename,
	}	

	output { 
		File shrunk_vcf = ShrinkVcfTask.vcf_out
		File shrunk_vcf_index = ShrinkVcfTask.vcf_out_index
	}
}

