version 1.0 

task fastqtl_nominal {
    input {
        File expression_bed
        File expression_bed_index
        File vcf
        File vcf_index
        String prefix
        File covariates

        String? cis_window
        Int? ma_sample_threshold
        Float? maf_threshold
        Int chunks

        Int memory=10
        Int disk_space
        Int num_threads
        Int num_preempt
    }

    command <<<
        set -euo pipefail

        wget https://raw.githubusercontent.com/broadinstitute/palantir-workflows/main/Scripts/monitoring/cromwell_monitoring_script.sh
        bash ./cromwell_monitoring_script.sh | tee monitoring.log &

        touch ~{vcf_index}  # avoid tabix "index older than vcf" error
        touch ~{expression_bed_index}
        # nominal pass
        /opt/fastqtl/python/run_FastQTL_threaded.py ~{vcf} ~{expression_bed} ~{prefix} \
            --covariates ~{covariates} \
            ~{"--window " + cis_window} \
            ~{"--ma_sample_threshold " + ma_sample_threshold} \
            ~{"--maf_threshold " + maf_threshold} \
            --chunks ~{chunks} \
            --threads ~{num_threads}
    >>>

    runtime {
        docker: "gcr.io/broad-cga-francois-gtex/gtex_eqtl:V8"
        memory: "~{memory}GB"
        disks: "local-disk ~{disk_space} HDD"
        cpu: "~{num_threads}"
        preemptible: "~{num_preempt}"
    }

    output {
        File allpairs="~{prefix}.allpairs.txt.gz"
        File allpairs_log="~{prefix}.allpairs.log"
        File monitoring_log="monitoring.log"
    }

    meta {
        author: "Francois Aguet"
    }
}


task fastqtl_permutations_scatter {
    input {
        File expression_bed
        File expression_bed_index
        File vcf
        File vcf_index
        String prefix
        File covariates

        Int current_chunk
        Int chunks
        Int permutations
        String? cis_window
        File? phenotype_groups
        Int? ma_sample_threshold
        Float? maf_threshold

        Int memory
        Int disk_space
        Int num_threads
        Int num_preempt
    }

    command <<<
        set -euo pipefail
        touch ~{vcf_index}  # avoid tabix "index older than vcf" error
        touch ~{expression_bed_index}
        # permutation pass
        /opt/fastqtl/python/run_chunk.py ~{vcf} ~{expression_bed} ~{prefix} ~{current_chunk} ~{chunks}\
            --permute ~{permutations} \
            --covariates ~{covariates} \
            ~{"--window " + cis_window} \
            ~{"--phenotype_groups " + phenotype_groups} \
            ~{"--ma_sample_threshold " + ma_sample_threshold} \
            ~{"--maf_threshold " + maf_threshold}
        mv ~{prefix}_chunk*.txt.gz ~{prefix}_chunk_~{current_chunk}.txt.gz
        mv ~{prefix}_chunk*.log ~{prefix}_chunk_~{current_chunk}.log
    >>>

    runtime {
        docker: "gcr.io/broad-cga-francois-gtex/gtex_eqtl:V8"
        memory: "~{memory}GB"
        disks: "local-disk ~{disk_space} HDD"
        cpu: "~{num_threads}"
        preemptible: "~{num_preempt}"
    }

    output {
        File chunk="~{prefix}_chunk_~{current_chunk}.txt.gz"
        File chunk_log="~{prefix}_chunk_~{current_chunk}.log"
    }

    meta {
        author: "Francois Aguet"
    }
}


task fastqtl_permutations_merge {
    input {
        Array[File] chunks
        Array[File] logs
        String prefix

        Int memory
        Int disk_space
        Int num_threads
        Int num_preempt

        Float? qvalue_lambda
    }
    command <<<
        set -euo pipefail
        /opt/fastqtl/python/merge_chunks.py ~{write_lines(chunks)} ~{write_lines(logs)} ~{prefix}\
            --permute ~{"--qvalue_lambda " + qvalue_lambda}
    >>>

    runtime {
        docker: "gcr.io/broad-cga-francois-gtex/gtex_eqtl:V8"
        memory: "~{memory}GB"
        disks: "local-disk ~{disk_space} HDD"
        cpu: "~{num_threads}"
        preemptible: "~{num_preempt}"
    }

    output {
        File genes="~{prefix}.genes.txt.gz"
        File genes_log="~{prefix}.genes.log"
    }

    meta {
        author: "Francois Aguet"
    }
}


task fastqtl_postprocess {
    input {
        File permutations_output
        File nominal_output
        Float fdr
        File annotation_gtf
        String prefix
        File? variant_lookup

        Int memory
        Int disk_space
        Int num_threads
        Int num_preempt
    }
    command <<<
        set -euo pipefail
        # it seems that fatqtl_postprocess needs the gtf to be uncompressed...
        gunzip -c < ~{annotation_gtf} > annotation.gtf 

        # post-processing
        /opt/fastqtl/python/annotate_outputs.py ~{permutations_output} ~{fdr} annotation.gtf \
        --nominal_results ~{nominal_output} \
        ~{"--snp_lookup " + variant_lookup}

    >>>

    runtime {
        docker: "gcr.io/broad-cga-francois-gtex/gtex_eqtl:V8"
        memory: "~{memory}GB"
        disks: "local-disk ~{disk_space} HDD"
        cpu: "~{num_threads}"
        preemptible: "~{num_preempt}"
    }

    output {
        File genes_annotated="~{prefix}.genes.annotated.txt.gz"
        File signifpairs="~{prefix}.signifpairs.txt.gz"
    }

    meta {
        author: "Francois Aguet"
    }
}


workflow fastqtl_workflow {
    input {
        File expression_bed
        File expression_bed_index
        File vcf
        File vcf_index
        String prefix
        File covariates

        Int permutations
        Int chunks
        String? cis_window
        Int? ma_sample_threshold
        Float? maf_threshold

        # post-processing
        Float fdr
        File annotation_gtf
        File? variant_lookup
    }
    call fastqtl_nominal {
        input:
            chunks=chunks, 
            prefix=prefix,
            expression_bed=expression_bed, 
            expression_bed_index=expression_bed_index, 
            vcf=vcf, 
            vcf_index=vcf_index, 
            covariates=covariates, 
            cis_window=cis_window,
            ma_sample_threshold=ma_sample_threshold, 
            maf_threshold=maf_threshold,
            disk_space=ceil(3*size(vcf,"GB")+200),
            num_threads=4,
            num_preempt=1
    }

    scatter(i in range(chunks)) {
        call fastqtl_permutations_scatter {
            input:
                current_chunk=i+1, 
                chunks=chunks, 
                prefix=prefix, 
                permutations=permutations,
                expression_bed=expression_bed, 
                expression_bed_index=expression_bed_index,
                vcf=vcf, 
                vcf_index=vcf_index,
                covariates=covariates, 
                cis_window=cis_window,
                ma_sample_threshold=ma_sample_threshold, 
                maf_threshold=maf_threshold,
                memory=10,
                disk_space=ceil(size(vcf,"GB")+20),
                num_threads=4,
                num_preempt=1
        }
    }

    call fastqtl_permutations_merge {
        input: 
            chunks=fastqtl_permutations_scatter.chunk, 
            logs=fastqtl_permutations_scatter.chunk_log, 
            prefix=prefix,
            memory=10,
            disk_space=50,
            num_threads=4,
            num_preempt=1
        }

    output {
        File genes_permutation=fastqtl_permutations_merge.genes
        File genes_nominal=fastqtl_nominal.allpairs

   }
}
