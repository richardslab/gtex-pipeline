version 1.0
##import "write_array.wdl" as wa
##import "Error.wdl" as e

task ConcatVcfsTask {
	input {
		Array[File] vcfs
		String basename
	}
	File files=write_lines(vcfs)
	command <<<
		wget https://raw.githubusercontent.com/broadinstitute/palantir-workflows/main/Scripts/monitoring/cromwell_monitoring_script.sh
		bash ./cromwell_monitoring_script.sh | tee monitoring.log &

		set -euo pipefail

		bcftools concat -f ~{files} -v 1 -n -Oz -o ~{basename}.vcf.gz 
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
		cpu: "1"
	}
}


workflow ConcatVcfs{
	input {
		Array[File] vcfs_in
		String basename
	}

	call ConcatVcfsTask {
		input:
			vcfs=vcfs_in,
			basename=basename,
	}
	
	output { 
		File vcf_cat = ConcatVcfsTask.vcf_out
		File vcf_cat_index = ConcatVcfsTask.vcf_out_index
	}
}

