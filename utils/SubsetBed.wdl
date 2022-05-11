version 1.0


task SubsetBed {
	input {
		File bed
		File contig_list
	}
	command<<<
#		set -euxo pipefail
		set -x
		ls
		file ~{bed}

		gunzip -c ~{bed} | head -n1 > out.bed
		gunzip -c ~{bed} | tail -n+2 | grep -w -F -f ~{contig_list} >> out.bed 

		bgzip out.bed
		tabix -p bed out.bed.gz
		ls

	>>>
	runtime {
		docker: "gcr.io/broad-cga-francois-gtex/gtex_eqtl:V8"	
		memory: "1GB"
		disks: "local-disk 20 HDD"
		cpu: 1
		preemptible: 1
	}

	output {
		File out = "out.bed.gz"
		File out_index = "out.bed.gz.tbi"

	}

	meta {
		author: "Yossi Farjoun"
	}
}

workflow SubsetBed_WF{
	input {
		File bed
		File contig_list
	}

	call SubsetBed {
		input:
			bed=bed,
			contig_list=contig_list
	}

	output {
		File bed_out=SubsetBed.out
		File bed_out_index=SubsetBed.out_index
	}
}
