version 1.0

task markduplicates {
    input { 
        File input_bam
        String prefix
        Int? max_records_in_ram
        Float? sorting_collection_size_ratio

        Float memory
        Int disk_space
        Int num_threads
        Int num_preempt
    }

    String output_bam = sub(basename(input_bam), "\\.bam$", ".md.bam")

    command <<<
        set -euo pipefail
        # taking memory from the variable so that memory increase can happen.
        java_memory=$(( ${MEMORY_SIZE} - 0.5 ))
        
        python3 -u /src/run_MarkDuplicates.py ~{input_bam} ~{prefix} \
            --memory ~{java_memory} \
            ~{"--max_records_in_ram " + max_records_in_ram} \
            ~{"--sorting_collection_size_ratio " + sorting_collection_size_ratio}
        samtools index ~{output_bam}
    >>>

    output {
        File bam_file = "~{output_bam}"
        File bam_index = "~{output_bam}.bai"
        File metrics = "~{prefix}.marked_dup_metrics.txt"
    }

    runtime {
        docker: "gcr.io/broad-cga-francois-gtex/gtex_rnaseq:V10"
        memory: "~{memory}GB"
        disks: "local-disk ~{disk_space} HDD"
        cpu: "~{num_threads}"
        preemptible: "~{num_preempt}"
        maxRetries: 3
    }

    meta {
        author: "Francois Aguet"
    }
}


workflow markduplicates_workflow {
    call markduplicates
}
