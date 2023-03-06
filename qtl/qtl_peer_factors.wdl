version 1.0


task get_num_peers_needed {
    input {
        Int num_samples
    }
    command <<<
    python <<EOF
    
    # From the gtex readme: https://github.com/richardslab/gtex-pipeline/blob/master/qtl/README.md

    # 15 factors for N < 150
    # 30 factors for 150 ≤ N < 250
    # 45 factors for 250 ≤ N < 350
    # 60 factors for N ≥ 350

    num_samples = ~{num_samples}
    if num_samples < 150:
        num_peers = 15
    elif num_samples < 250:
        num_peers = 30
    elif num_samples < 350:
        num_peers = 45
    else:
        num_peers = 60

    print(num_peers)

    EOF

    >>>

    output {
        Int num_peers=read_int(stdout())
    }

    runtime {
        docker: "python:latest"
        memory: "2GB"
        disks: "local-disk 15 HDD"
    }
}

task qtl_peer_factors {
    input {
        File expression_file
        String prefix
        Int num_peer

        File? genotype_pcs
        File? add_covariates

        Int memory
        Int disk_space 
        Int num_threads = 1
        Int num_preempt = 1
    }

    command <<<
        set -xeuo pipefail

        cat -<< "EOF" > combine_covariates.py
        #!/usr/bin/env python3
        # Author: Francois Aguet
        import pandas as pd
        import numpy as np
        import argparse
        import os

        parser = argparse.ArgumentParser(description='Combine covariates into a single matrix')
        parser.add_argument('expression_covariates', help='')
        parser.add_argument('prefix', help='')
        parser.add_argument('--genotype_pcs', default=None, help='Genotype PCs')
        parser.add_argument('--add_covariates', default=[], nargs='+', help='Additional covariates')
        parser.add_argument('-o', '--output_dir', default='.', help='Output directory')
        args = parser.parse_args()

        print('Combining covariates ... ', end='', flush=True)
        expression_df = pd.read_csv(args.expression_covariates, sep='\t', index_col=0, dtype=str)
        print(f"there are {len(expression_df)} expression covariates.")
        if args.genotype_pcs is not None:
            print(f"Reading genotype PCs {args.genotype_pcs}")
            genotype_df = pd.read_csv(args.genotype_pcs, sep='\t', index_col=0, dtype=str)
            print(f"There are {len(genotype_df)} genotype PCs.")
            combined_df = pd.concat([genotype_df[expression_df.columns], expression_df], axis=0)
            print(f"There are {len(combined_df)} combined covariates.")
        else:
            combined_df = expression_df
        for c in args.add_covariates:
            print(f"Reading additional covariates {args.genotype_pcs}.")
            additional_df = pd.read_csv(c, sep='\t', index_col=0, dtype=str)
            print(f"There are {len(additional_df)} additional covariates in this file.")
            combined_df = pd.concat([combined_df, additional_df[expression_df.columns]], axis=0)
            print(f"There are now {len(combined_df)} combined covariates.")

        # identify and drop colinear covariates
        C = combined_df.astype(np.float64).T
        Q,R = np.linalg.qr(C-np.mean(C, axis=0))
        colinear_ix = np.abs(np.diag(R)) < np.finfo(np.float64).eps * C.shape[1]
        if np.any(colinear_ix):
            print('Colinear covariates detected:')
            for i in C.columns[colinear_ix]:
                print("  * dropped '{}'".format(i))
            combined_df = combined_df.loc[~colinear_ix]

        combined_df.to_csv(os.path.join(args.output_dir, args.prefix+'.combined_covariates.txt'), sep='\t')#, float_format='%.6g')
        print('done.')

        EOF

        Rscript /src/run_PEER.R ~{expression_file} ~{prefix} ~{num_peer}
        env python3 ./combine_covariates.py ~{prefix}.PEER_covariates.txt ~{prefix} ~{"--genotype_pcs " + genotype_pcs} ~{"--add_covariates " + add_covariates}
    >>>

    runtime {
        docker: "gcr.io/broad-cga-francois-gtex/gtex_eqtl:V8"
        memory: "~{memory}GB"
        disks: "local-disk ~{disk_space} HDD"
        cpu: "~{num_threads}"
        preemptible: "~{num_preempt}"
    }

    output {
        File combined_covariates="~{prefix}.combined_covariates.txt"
        File PEER_covariates="~{prefix}.PEER_covariates.txt"
        File PEER_residuals="~{prefix}.PEER_residuals.txt"
        File PEER_alpha="~{prefix}.PEER_alpha.txt"
    }

    meta {
        author: "Francois Aguet"
    }
}

workflow qtl_peer_factors_workflow {
    input {
        Array[String]? samples_for_count
        Int? num_peers
    }

    if (!defined(num_peers)){
        ## should fail is samples_for_count is missing
        Array[String] samples_for_count_def = select_first([samples_for_count])
        
        call get_num_peers_needed{
            input:
                num_samples=length(samples_for_count_def)
        }
    }

    call qtl_peer_factors{
        input:
            num_peer=select_first([num_peers, get_num_peers_needed.num_peers]),
            disk_space = 20,
            memory = 4
    }


    output {
        File combined_covariates=qtl_peer_factors.combined_covariates
        File alpha=qtl_peer_factors.PEER_alpha
    }
}
