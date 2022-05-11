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
        set -euo pipefail
        Rscript /src/run_PEER.R ~{expression_file} ~{prefix} ~{num_peer}
        /src/combine_covariates.py ~{prefix}.PEER_covariates.txt ~{prefix} ~{"--genotype_pcs " + genotype_pcs} ~{"--add_covariates " + add_covariates}
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
