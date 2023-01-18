version 1.0
##import "write_array.wdl" as wa
##import "Error.wdl" as e

task ConcatVcfsTask {
	input {
		Array[File] vcfs
		String basename
	}
	
	command <<<
		wget https://raw.githubusercontent.com/broadinstitute/palantir-workflows/main/Scripts/monitoring/cromwell_monitoring_script.sh
		bash ./cromwell_monitoring_script.sh | tee monitoring.log &

		set -euo pipefail

		bcftools concat -n -Oz -o ~{basename}.vcf.gz ~{sep=" " vcfs}
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

task UpdateSequencesFromFai {
	input {
		File vcf_in
		File fasta_fai
	}
	String vcf_out_name=basename(vcf_in,".vcf.gz")+".reheadered.vcf.gz"
	command <<<
		wget https://raw.githubusercontent.com/broadinstitute/palantir-workflows/main/Scripts/monitoring/cromwell_monitoring_script.sh
		bash ./cromwell_monitoring_script.sh | tee monitoring.log &

		set -euo pipefail
		
		bcftools reheader -f ~{fasta_fai} -Oz -o ~{vcf_out_name} vcf_in

	>>>
	output {
		File vcf_out=vcf_out_name
		File monitoring_log="monitoring.log"
	}

	runtime {
		docker: "bschiffthaler/bcftools:latest"
		preemptible: 0
		disks: "local-disk " + (2*ceil(size(vcf_in,"GiB"))+20) + " HDD"
		bootDiskSizeGb: "16"
		memory: 20 + " GB"
		cpu: "1"
	}
}

workflow ConcatVcfs{
	input {
		Array[File] vcfs_in
		File? fasta_index #optional_index_for_header
		String basename
	}

	if (defined(fasta_index)){
		scatter(vcf in vcfs_in){
			call UpdateSequencesFromFai{input:
				vcf_in=vcf,
				fasta_fai=select_first([fasta_index])
			}
			
		}
		Array[File] vcf_out_files=UpdateSequencesFromFai.vcf_out
	}

	Array[File] vcfs_to_use=select_first([vcf_out_files,vcfs_in])
	call ConcatVcfsTask {
		input:
			vcfs=vcfs_to_use,
			basename=basename
	}
	
	output { 
		File vcf_cat = ConcatVcfsTask.vcf_out
		File vcf_cat_index = ConcatVcfsTask.vcf_out_index
	}
}

