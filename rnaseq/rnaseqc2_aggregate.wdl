version 1.0

task rnaseqc2_aggregate {
    input {
        Array[File]+ tpm_gcts
        Array[File]+ count_gcts
        Array[File]+ exon_count_gcts
        Array[File]+ metrics_tsvs
        String prefix
        Array[File]? insertsize_hists
        String? flags

        Int memory
        Int disk_space
        Int num_threads
        Int num_preempt
    }
    command <<<
        set -euo pipefail
        date +"[%b %d %H:%M:%S] Aggregating RNA-SeQC outputs"

        if [ ~{length(tpm_gcts)} != ~{length(count_gcts)} ] || \
           [ ~{length(tpm_gcts)} != ~{length(count_gcts)} ] || \
           [ ~{length(tpm_gcts)} != ~{length(exon_count_gcts)} ] || \
           [ ~{length(tpm_gcts)} != ~{length(metrics_tsvs)} ] ; then 
            cat << EOF 
            some of the lengths of the inputs are not equal:
            
            count_gcts: ~{length(count_gcts)}
            count_gcts: ~{length(count_gcts)}
            exon_count_gcts: ~{length(exon_count_gcts)}
            metrics_tsvs: ~{length(metrics_tsvs)}
        EOF
        
        fi

        mkdir individual_outputs
        mv ~{sep=' ' tpm_gcts} individual_outputs/
        mv ~{sep=' ' count_gcts} individual_outputs/
        mv ~{sep=' ' exon_count_gcts} individual_outputs/
        mv ~{sep=' ' metrics_tsvs} individual_outputs/
        if [ -n '~{sep=' ' insertsize_hists}' ]; then
            mv ~{sep=' ' insertsize_hists} individual_outputs/
        fi
        touch ~{prefix}.insert_size_hists.txt.gz
        python3 -m rnaseqc aggregate \
            -o . \
            individual_outputs \
            ~{prefix} \
            ~{flags}
        date +"[%b %d %H:%M:%S] done"
    >>>

    output {
        File metrics="~{prefix}.metrics.txt.gz"
        File insert_size_hists="~{prefix}.insert_size_hists.txt.gz"
        File tpm_gct=glob("~{prefix}.gene_tpm.*")[0]
        File count_gct=glob("~{prefix}.gene_reads.*")[0]
        File exon_count_gct=glob("~{prefix}.exon_reads.*")[0]
    }

    runtime {
        docker: "gcr.io/broad-cga-francois-gtex/gtex_rnaseq:V10"
        memory: "~{memory}GB"
        disks: "local-disk ~{disk_space} HDD"
        cpu: "~{num_threads}"
        preemptible: "~{num_preempt}"
    }

    meta {
        author: "Francois Aguet"
    }
}


workflow rnaseqc2_aggregate_workflow {
    input {
        Array[File]+ tpm_gcts
        Array[File]+ count_gcts
        Array[File]+ exon_count_gcts
        Array[File]+ metrics_tsvs
        String prefix
        Array[File]? insertsize_hists
        String? flags

        Int memory=20
        Int? disk_space
        Int num_threads=1
        Int num_preempt=1
    }
    call rnaseqc2_aggregate{
        input:
            tpm_gcts=tpm_gcts,
            count_gcts=count_gcts,
            exon_count_gcts=exon_count_gcts,
            metrics_tsvs=metrics_tsvs,
            prefix=prefix,
            insertsize_hists=insertsize_hists,
            flags=flags,
            memory=memory,
            disk_space=select_first([disk_space, ceil(size(flatten([tpm_gcts,count_gcts,exon_count_gcts,metrics_tsvs]),"GB")*3+20)]),
            num_threads=num_threads,
            num_preempt=num_preempt
    }
    output {
        File metrics=rnaseqc2_aggregate.metrics
        File insert_size_hists=rnaseqc2_aggregate.insert_size_hists
        File tpm_gct=rnaseqc2_aggregate.tpm_gct
        File count_gct=rnaseqc2_aggregate.count_gct
        File exon_count_gct=rnaseqc2_aggregate.exon_count_gct    
    }
}
