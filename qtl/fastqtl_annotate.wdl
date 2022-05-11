version 1.0 
import "fastqtl.wdl" as fq

workflow fastqtl_annotate {
    input {
        File permutations
        File nominal
        String prefix
        Float fdr
        File annotation_gtf
        File? variant_lookup

    }
    
    call fq.fastqtl_postprocess as fastqtl_postprocess{
        input:
            permutations_output=permutations,
            nominal_output=nominal,
            prefix=prefix, 
            fdr=fdr, 
            annotation_gtf=annotation_gtf, 
            variant_lookup=variant_lookup,
            memory=20,
            disk_space=ceil(20+2*(size(permutations,"GB")+size(nominal,"GB"))),
            num_threads=4,
            num_preempt=1
    }

    output {
        
        File genes_annotated=fastqtl_postprocess.genes_annotated
        File signifpairs=fastqtl_postprocess.signifpairs
   }
}
