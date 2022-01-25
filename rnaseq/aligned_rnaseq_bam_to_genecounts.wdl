version 1.0

import "markduplicates.wdl" as MD
import "rnaseqc2.wdl" as RNAQC

workflow rnaseq_pipeline_bam_workflow {
    input {
        File bam_file
        File bam_index
        String prefix
    }

    call MD.markduplicates as markduplicates {
    input: 
        input_bam=bam_file, 
        prefix=prefix
    }

    call RNAQC.rnaseqc2 as rnaseqc2 {
    input: 
        bam_file=markduplicates.bam_file, 
        sample_id=prefix
    }

    output {
        File deduped_bam = markduplicates.bam_file
        File deduped_bam_index = markduplicates.bam_index
        File duplication_metrics = markduplicates.metrics
        File gene_tpm = rnaseqc2.gene_tpm
        File gene_counts = rnaseqc2.gene_counts
        File exon_counts = rnaseqc2.exon_counts
        File metrics = rnaseqc2.metrics
        File gc_content = rnaseqc2.gc_content
        File insertsize_distr = rnaseqc2.insertsize_distr
    }
}
