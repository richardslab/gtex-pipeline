version 1.0

task ConvertPlinkToVcf {
	input {
		File bim
		File bed
		File fam
	}

	String outbase=basename(bim, '.bim')
	command <<<
		set -euo pipefail

		# make sure that the files all have the same basename
		# will not work if names have spaces in them

		# the difference between using basename and using substitution is that substitution will 
		# maintain the leading path, which is desirable sometimes.

		bim_bash=~{bim}
		bed_bash=~{bed}
		fam_bash=~{fam}

		bim_base=${bim_bash%.bim}
		bed_base=${bed_bash%.bed}
		fam_base=${fam_bash%.fam}

		lines=$(echo $bim_base $bed_base $fam_base | tr ' ' '\n'| uniq | wc -l)

		if [ "${lines}" -eq "1" ]; then
			bim_bash=~{bim}
			plink --bfile "${bim_base}" --recode vcf-iid --out "~{outbase}"
		else
			echo "Error, found too many basenames in the input: $lines" 
			exit 1
		fi

	>>>
	output {
		File vcf="~{outbase}.vcf"
	}


	runtime {
			docker: "dnastack/plink:1.9"
			preemptible: 0
			disks: "local-disk " + ceil(size([bim,bed,fam],"GiB")+20) + " HDD"
			bootDiskSizeGb: "16"
			memory: 20 + " GB"
	}
}


task RenameChrXAndSubsetToSNPs {

	input {
		File vcf_in
		File ref_index
	}
	String outbase=basename(vcf_in, '.vcf')
	
	command <<<
		set -exuo pipefail 

		mkdir temp
		# remove lines that start with '##contig=<ID=' and for the remaining lines that
		# to not start with '#', replace 23 with X and add 'chr' to the begining of the line.
		grep -v '^##contig=<ID=' "~{vcf_in}" | \
			sed  '/#/!{s/^23\t/X\t/; s/^/chr/}'  > raw.vcf
		rm ~{vcf_in}

		bcftools reheader -f ~{ref_index} -o reheadered.unsorted.snps.vcf raw.vcf
		rm raw.vcf

		bcftools view --no-update -v snps -e 'REF=="-"||ALT=="-" || REF=="."||ALT=="."' -Oz -o unsorted.snps.vcf.gz reheadered.unsorted.snps.vcf
		rm reheadered.unsorted.snps.vcf

		bcftools sort -Oz -o "~{outbase}".snps.vcf.gz unsorted.snps.vcf.gz 
		rm unsorted.snps.vcf.gz


		bcftools index -t "~{outbase}".snps.vcf.gz 
	>>>

	output {
		File vcf="~{outbase}.snps.vcf.gz"
		File vcf_index="~{outbase}.snps.vcf.gz.tbi"
	}

	runtime {
			docker: "bschiffthaler/bcftools:latest"
			preemptible: 0
			disks: "local-disk " + (4*ceil(size(vcf_in,"GiB"))+20) + " HDD"
			bootDiskSizeGb: "16"
			memory: 20 + " GB"
	}

}

task ReheaderVcf{
	input {
		File vcf_in
		File vcf_index_in
		File ref_fasta
		File ref_dict
		File ref_index
		String basename
	}

	command <<<
		set -euo pipefail

		java -jar picard.jar UpdateVcfSequenceDictionary \
			-R "~{ref_fasta}" \
			-I "~{vcf_in}" \
			-O "~{basename}.vcf.gz" \
			-SD "~{ref_dict}"

		tabix "~{basename}.vcf.gz"
	>>>

	output {
		File vcf="~{basename}.vcf.gz"
		File vcf_index="~{basename}.vcf.gz.tbi"
	}

	runtime {
			docker: "broadinstitute/picard:2.26.8"
			preemptible: 0
			disks: "local-disk " + (ceil(size([vcf_in,vcf_in,ref_fasta],"GiB"))+20) + " HDD"
			bootDiskSizeGb: "16"
			memory: 20 + " GB"
	}
}

workflow ConvertPlinkToVcfWF {
	call ConvertPlinkToVcf{}

	call RenameChrXAndSubsetToSNPs{
		input:
			vcf_in = ConvertPlinkToVcf.vcf
	}

#	call ReheaderVcf{
#		input:
#		vcf_in=RenameChrXAndSubsetToSNPs.vcf,
#		vcf_index_in=RenameChrXAndSubsetToSNPs.vcf_index
#	}

	output {
		File vcf=RenameChrXAndSubsetToSNPs.vcf
		File vcf_index=RenameChrXAndSubsetToSNPs.vcf_index
	}
}