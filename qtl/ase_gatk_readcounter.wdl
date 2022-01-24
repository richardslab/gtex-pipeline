version 1.0

task ase_gatk_readcounter {
    input {
        File gatk_jar
        File genome_fasta
        File genome_fasta_index
        File genome_fasta_dict
        File het_vcf
        File het_vcf_index
        File bam_file
        File bam_index
        String prefix
        Boolean filter_wasp = false

        Int memory 
        Int disk_space
        Int num_threads = 1 
        Int num_preempt = 1
    }

    command <<<
        set -euo pipefail
        if [[ ~{filter_wasp} = "true" ]]
        then
            date +"[%b %d %H:%M:%S] Filtering out reads with allelic mapping bias"
            samtools view -h ~{bam_file} | grep -v "vW:i:[2-7]" | samtools view -1 > filtered.bam
            samtools index filtered.bam
        else
            ln -s ~{bam_file} filtered.bam
            ln -s ~{bam_index} filtered.bam.bai
        fi
        
        python3 /src/run_GATK_ASEReadCounter.py ~{gatk_jar} ~{genome_fasta} ~{het_vcf} filtered.bam ~{prefix}
        

        # filter out chrX and chrY?
        mv ~{prefix}.readcounts.txt.gz ~{prefix}.readcounts.all.txt.gz
        zcat ~{prefix}.readcounts.all.txt.gz | awk '$1!="chrX" && $1!="X" $1!="chrY" && $1!="Y" {print $0}' | gzip -c > ~{prefix}.readcounts.txt.gz
        zcat ~{prefix}.readcounts.all.txt.gz | awk '$1=="contig" || $1=="chrX" || $1=="X"  {print $0}' | gzip -c > ~{prefix}.readcounts.chrX.txt.gz
        zcat ~{prefix}.readcounts.all.txt.gz | awk '$1=="contig" || $1=="chrY" || $1=="Y"  {print $0}' | gzip -c > ~{prefix}.readcounts.chrY.txt.gz
    >>>

    output {
        File ase_read_counts = "${prefix}.readcounts.txt.gz"
        File ase_read_counts_chrX = "${prefix}.readcounts.chrX.txt.gz"
        File ase_read_counts_chrY = "${prefix}.readcounts.chrY.txt.gz"
    }

    runtime {
        docker: "gcr.io/broad-cga-francois-gtex/gtex_eqtl:V8"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_threads}"
        preemptible: "${num_preempt}"
    }

    meta {
        author: "Francois Aguet"
    }
}


workflow ase_gatk_readcounter_workflow {
    call ase_gatk_readcounter
}
