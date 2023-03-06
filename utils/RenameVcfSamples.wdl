version 1.0
import "write_array.wdl" as wa
import "Error.wdl" as e
task RenameVcfSamples{
	input {
		File vcf
		File rename_by
	}
	command <<<
		set -euo pipefail 
		bcftools reheader -s ~{rename_by} -o output.vcf.gz ~{vcf}
		bcftools index -t output.vcf.gz 

	>>>
	output {
		File vcf_out="output.vcf.gz"
		File vcf_out_index="output.vcf.gz.tbi"
	}

	runtime {
		docker: "bschiffthaler/bcftools:latest"
		preemptible: 0
		disks: "local-disk " + (2*ceil(size(vcf,"GiB"))+20) + " HDD"
		bootDiskSizeGb: "16"
		memory: 20 + " GB"
	}
}


workflow RenameVcfSamplesWF{
	input {
		File vcf_in
		Array[String] samples_in_vcf
		Array[String] new_sample_names
	}

	if(length(samples_in_vcf)!=length(new_sample_names)){
		call e.Error {
			input:
				message="Length of arrays do not match: ~{length(samples_in_vcf)} != ~{length(new_sample_names)}",
				error=1
		}
	}
	Array[Array[String]] transposed=transpose([samples_in_vcf,new_sample_names])
	call wa.write_array_to_tsv as array{
		input:
			array=transposed
	}

	call RenameVcfSamples{
		input:
			vcf=vcf_in,
			rename_by=array.tsv
	}

	output{ 
		File vcf_renamed = RenameVcfSamples.vcf_out
		File vcf_renamed_index = RenameVcfSamples.vcf_out_index}
}