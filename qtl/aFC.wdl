version 1.0

task aFC {
    input {
        File vcf_file
        File vcf_index
        File expression_bed
        File expression_bed_index
        File covariates_file
        File afc_qtl_file
        String prefix
        String docker="gcr.io/broad-cga-francois-gtex/gtex_eqtl:V8"
        
        Int memory
        Int disk_space
        Int num_threads
        Int num_preempt

        Boolean raw=false
    }
    
    command <<<
        set -euo pipefail
        python3 /opt/aFC/aFC.py \
            --vcf ~{vcf_file} \
            --pheno ~{expression_bed} \
            --qtl ~{afc_qtl_file} \
            --cov ~{covariates_file} \
            ~{true="--log_xform 0" false="--log_xform 1 --log_base 2" raw} \
            --output_se \
            --o ~{prefix}.aFC.txt
        gzip ~{prefix}.aFC.txt
    >>>

    runtime {
        docker: "~{docker}"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_threads}"
        preemptible: "${num_preempt}"
    }

    output {
        File afc_file="${prefix}.aFC.txt.gz"
    }

    meta {
        author: "Francois Aguet"
    }
}

task convert_qtls{
    input {
        File fastQTL_output
        String prefix

        Int memory
        Int disk_space
        Int num_preempt

        Int scattercount
    }

    command <<<
    set -xeuo pipefail 

    cat -<< "EOF" > extract_qtl_ids.py
    #!/usr/bin/env python3
    import pandas as pd
    import argparse
    import numpy as np
    import os

    parser = argparse.ArgumentParser(description='Extract variant-gene pairs from list of associations.')
    parser.add_argument('--input_pairs', help="output from FastQTL.")
    parser.add_argument('--prefix', help='Prefix for output file: <prefix>.extracted_qtls.txt')
    parser.add_argument('-o', '--output_dir', default='.', help='Output directory')
    parser.add_argument('-s', '--scatter_count', default=1, help='Scatter count.')
    args = parser.parse_args()

    print('Loading input')
    ref_pairs_df = pd.read_csv(args.input_pairs, sep='\t', usecols=['variant_id', 'gene_id'], dtype=str)
    ref_pairs_df[["sid_chr", "sid_pos", "sid_ref", "sid_alt"]] = ref_pairs_df['variant_id'].str.split(":", expand=True)
    ref_pairs_df = ref_pairs_df.rename(columns={"variant_id": "sid", "gene_id": "pid"}).reindex(
        columns=["sid", "pid", "sid_chr", "sid_pos"])[["sid", "pid", "sid_chr", "sid_pos"]]
    split_df=np.split(ref_pairs_df,args.scattercount)
    for i,df in enumerate(split_df):

        with open(os.path.join(args.output_dir, f"{args.prefix}_{i}_extracted_qtls.txt"), 'wt') as f:
            df.to_csv(f, sep='\t', na_rep='NA', float_format='%.6g', index=False)

    EOF
    python3 extract_qtl_ids.py --input_pairs "~{fastQTL_output}" --prefix "~{prefix}" --scattercount ~{scattercount}
 
    >>>
    runtime {
        docker: "gcr.io/broad-cga-francois-gtex/gtex_eqtl:V8"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: 1
        preemptible: "${num_preempt}"
    }

    output {
        Array[File] qtl_files=glob("${prefix}_*_extracted_qtls.txt")
    }

    meta {
        author: "Yossi Farjoun"
    }
}

workflow aFC_workflow {
    input {
        Int scattercount=1
    }

    call convert_qtls{
        input:
            scattercount=scattercount
    }
    scatter (qtl_file in convert_qtls.qtl_files){
        call aFC {
            input:
            afc_qtl_file=qtl_file
        }
    }
    output {
        Array[File] afc_file=aFC.afc_file
    }
}
