version 1.0
##import "write_array.wdl" as wa
##import "Error.wdl" as e

task ConcatVcfsTask {
	input {
		Array[File] vcfs
		Int threads=1
		String basename
	}
	File files=write_lines(vcfs)
	command <<<
		wget https://raw.githubusercontent.com/broadinstitute/palantir-workflows/main/Scripts/monitoring/cromwell_monitoring_script.sh
		bash ./cromwell_monitoring_script.sh | tee monitoring.log &

		set -euo pipefail

		bcftools concat  --threads ~{threads} -l ~{files} -Oz -o ~{basename}.vcf.gz 
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
		disks: "local-disk " + (2*ceil(size(vcfs,"GiB"))+20) + " HDD"
		bootDiskSizeGb: "16"
		memory: 20 + " GB"
		cpu: "~{threads}"
	}
}


workflow ConcatVcfs{
	input {
		Array[File] vcfs_in
		Int threads=1
		String basename
	}

	call ConcatVcfsTask {
		input:
			vcfs=vcfs_in,
			threads=threads,
			basename=basename,
	}
	
	output { 
		File vcf_cat = ConcatVcfsTask.vcf_out
		File vcf_cat_index = ConcatVcfsTask.vcf_out_index
	}
}

