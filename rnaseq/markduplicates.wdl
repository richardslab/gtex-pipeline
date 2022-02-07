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
        Float java_memory_discount=0.5
    }

    String output_bam = sub(basename(input_bam), "\\.bam$", ".md.bam")

    command <<<
        set -euo pipefail

        curl -L -O https://github.com/broadinstitute/palantir-workflows/raw/main/Scripts/monitoring/cromwell_monitoring_script.sh 
        chmod a+x cromwell_monitoring_script.sh 
        ./cromwell_monitoring_script.sh &

        # taking memory from the variable so that memory increase can happen.
        # awk is used as a workaround for 'bc' not being available 
        # (thanks https://stackoverflow.com/a/48534957/360496)
        a="~{java_memory_discount}"
        b="${MEM_SIZE}"
        java_memory=$(awk -v a="$a" -v b="$b" 'BEGIN { printf "%s\n", int(b-a) }' </dev/null )
        
        python3 -u /src/run_MarkDuplicates.py ~{input_bam} ~{prefix} \
            --memory "${java_memory}" \
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
        maxRetries: 1
    }

    meta {
        author: "Francois Aguet"
    }
}


workflow markduplicates_workflow {
    call markduplicates
}
