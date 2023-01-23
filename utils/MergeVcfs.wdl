version 1.0
##import "write_array.wdl" as wa
import "Error.wdl" as e

task MergeVcfsTask {
	input {
		Array[File] vcfs
		Int threads=1
		String basename
		String? region
	}
	command <<<
		wget https://raw.githubusercontent.com/broadinstitute/palantir-workflows/main/Scripts/monitoring/cromwell_monitoring_script.sh
		bash ./cromwell_monitoring_script.sh | tee monitoring.log &

		set -euo pipefail

		#remove format field (except genotype), and remove any site that is filtered
		bcftools merge  --threads ~{threads} -l ~{write_lines(vcfs)} -0 -f -F+ ~{"-r " + region} |\
			bcftools annotate  -x FORMAT -Oz -o ~{basename}.vcf.gz 

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


workflow MergeVcfs{
	input {
		Array[String] names
		Array[File] vcfs_in
		Int threads=1
		String basename
		String? region
	}

	if (length(names) !=length(vcfs_in)){
		call e.Error{input:
			message="length of names and files must be the same!",
			error=1
		}
		
	} 

	if (length(names) ==length(vcfs_in)){

		call MergeVcfsTask {
			input:
				vcfs=vcfs_in,
				threads=threads,
				basename=basename,
				region=region
		}
	}	

	output { 
		File? vcf_merged = MergeVcfsTask.vcf_out
		File? vcf_merged_index = MergeVcfsTask.vcf_out_index
	}
}

