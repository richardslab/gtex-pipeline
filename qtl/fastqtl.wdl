version 1.0 

task fastqtl_nominal {
    input {
        File check
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
        touch ~{prefix}.allpairs.log 
        touch monitoring.log
        
        # nominal pass

        cat -<< "EOF" > run_FastQTL_threaded.py
        #!/usr/bin/env python3
        # Author: Francois Aguet
        import argparse
        import os
        import numpy as np
        import subprocess
        import gzip
        import multiprocessing as mp
        import contextlib
        from datetime import datetime
        import tempfile
        import glob

        @contextlib.contextmanager
        def cd(cd_path):
            saved_path = os.getcwd()
            os.chdir(cd_path)
            yield
            os.chdir(saved_path)

        def get_cmd(args, chunk):
            cmd = os.path.join(fastqtl_dir, 'bin', 'fastQTL')+' --vcf '+args.vcf+' --bed '+args.bed+' --window '+args.window \
                +' --maf-threshold '+args.maf_threshold \
                +' --ma-sample-threshold '+args.ma_sample_threshold \
                +' --interaction-maf-threshold '+args.interaction_maf_threshold
            if args.covariates:
                cmd += ' --cov '+args.covariates
            if args.phenotype_groups:
                cmd += ' --grp '+args.phenotype_groups
            if args.threshold:
                cmd += ' --threshold '+args.threshold
            if args.permute:
                cmd += ' --permute '+' '.join([str(p) for p in args.permute])
            if args.interaction:
                cmd += ' --interaction '+args.interaction
            if args.best_variant_only:
                cmd += ' --report-best-only'
            if args.seed:
                cmd += ' --seed '+args.seed
            if args.exclude_samples:
                cmd += ' --exclude-samples '+args.exclude_samples
            if args.exclude_sites:
                cmd += ' --exclude-sites '+args.exclude_sites
            cmd += ' --chunk '+str(chunk)+' '+args.chunks\
                + ' --out '+args.prefix+'_chunk{0:03d}.txt.gz'.format(chunk)\
                + ' --log '+args.prefix+'_chunk{0:03d}.log'.format(chunk)
            return cmd

        def perm_worker(inputs):
            args = inputs[0]
            chunk = inputs[1]
            cmd = get_cmd(args, chunk)
            print('Processing chunk '+str(chunk), flush=True)
            s = subprocess.check_call(cmd, shell=True, executable='/bin/bash', stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            print('Finished chunk '+str(chunk), flush=True)
            return s


        parser = argparse.ArgumentParser(description='Run FastQTL')
        parser.add_argument('vcf', help='Genotypes in VCF 4.1 format')
        parser.add_argument('bed', help='Phenotypes in UCSC BED extended format')
        parser.add_argument('prefix', help='Prefix for output file name')
        parser.add_argument('--covariates', default='', help='Covariates')
        parser.add_argument('--phenotype_groups', default='', help='File with mapping of phenotype_id to group_id (gene_id)')
        parser.add_argument('--chunks', default='100', help='Number of chunks, minimum: #chromosomes')
        parser.add_argument('--permute', default=None, type=str, nargs='+', help='Number of permutations, e.g. [1000, 10000] (adaptive). Default: None (run nominal pass)')
        parser.add_argument('--interaction', default=None, type=str, help='Interaction term')
        parser.add_argument('--best_variant_only', action='store_true')
        parser.add_argument('--window', default='1e6', help='Cis-window size. Default values is 1Mb (1e6).')
        parser.add_argument('--threshold', default='', help='Output only significant phenotype-variant pairs with a p-value below threshold (default 1)')
        parser.add_argument('--maf_threshold', default='0.0', help='Include only genotypes with minor allele frequency >=maf_threshold (default 0)')
        parser.add_argument('--ma_sample_threshold', default='0', help='Include only genotypes with >=ma_sample_threshold samples carrying the minor allele (default 0)')
        parser.add_argument('--interaction_maf_threshold', default='0', help='MAF threshold for interactions, applied to lower and upper half of samples')
        parser.add_argument('--fdr', default=0.05, type=np.double)
        parser.add_argument('--seed', default=None, help='Random number generator seed')
        parser.add_argument('--exclude_samples', default=None, help='')
        parser.add_argument('--exclude_sites', default=None, help='')
        parser.add_argument('--qvalue_lambda', default=None, help='lambda parameter for pi0est in qvalue.')
        parser.add_argument('-t', '--threads', default=8, type=int, help='Number of threads')
        parser.add_argument('-o', '--output_dir', default='.', help='Output directory')
        args = parser.parse_args()
        fastqtl_dir =  "/opt/fastqtl"

        if not os.path.exists(args.output_dir):
            os.makedirs(args.output_dir)

        print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Running FastQTL on {0:d} threads.'.format(args.threads), flush=True)

        with cd(args.output_dir):
            with mp.Pool(processes=args.threads) as pool:
                pdata_res = [pool.map_async(perm_worker, ((args,k),)) for k in np.arange(1,int(args.chunks)+1)]
                pool.close()
                pool.join()
            print(pdata_res)

            for res in pdata_res:  # check exit status
                if res is None:
                    print("one of the worker threads returned a None!")
                    assert False
                try:
                    resget=res.get()
                except Exception as e:
                    print("An exception got raised in one of the jobs:")
                    print(res)
                    print(e)
                    raise(e)

                if resget is None:
                    print("one of the worker threads res.get() returned None!")
                    print(res)
                if len(res.get())==0:
                    print("one of the worker threads returned an empty array:")
                    print(res)
                    assert False  
                if res.get()[0]!=0:
                    print("one of the worker threads returned a non-zero value:")
                    print(res)
                    assert False



            with tempfile.NamedTemporaryFile(mode='w+') as chunk_list_file, \
                 tempfile.NamedTemporaryFile(mode='w+') as log_list_file:

                # write chunk and log paths to file
                chunk_files = sorted(glob.glob(args.prefix+'_chunk*.txt.gz'))
                chunk_list_file.write('\n'.join(chunk_files)+'\n')
                chunk_list_file.flush()
                log_files = sorted(glob.glob(args.prefix+'_chunk*.log'))
                log_list_file.write('\n'.join(log_files)+'\n')
                log_list_file.flush()

                # merge chunks
                cmd = 'python3 '+os.path.join(fastqtl_dir, 'python', 'merge_chunks.py') \
                    +' {} {} {} --fdr {} -o .'.format(chunk_list_file.name, log_list_file.name, args.prefix, args.fdr)
                if args.qvalue_lambda:
                    cmd += ' --qvalue_lambda {}'.format(args.qvalue_lambda)
                if args.permute:
                    cmd += ' --permute'
                subprocess.check_call(cmd, shell=True)

                # remove chunk files
                for f in chunk_files + log_files:
                    os.remove(f)
        EOF

        env python3 run_FastQTL_threaded.py ~{vcf} ~{expression_bed} ~{prefix} \
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
        Array[File] logs = glob("*.log")
        File? monitoring_log="monitoring.log"  #optional so that delocalization doesn't stop if missing
        File? allpairs_log="~{prefix}.allpairs.log"  #optional so that delocalization doesn't stop if missing
        File? allpairs="~{prefix}.allpairs.txt.gz" #optional so that delocalization doesn't stop if missing
    }

    meta {
        author: "Francois Aguet"
    }
}


task CheckInputs{
    input {
        File expression_bed
        File vcf
        File covariates
    }

  parameter_meta {
    expression_bed: {
      localization_optional: true
    }
    vcf: {
      localization_optional: true
    }
    covariates: {
      localization_optional: true
    }
  }

  command <<<

    set -ueo pipefail 


    gsutil cat ~{covariates} | head -n 1 | cut 2- > covariates.samples

    gsutil cat ~{expression_bed} | zcat | head -n 1 | cut 5- > expression.samples

    gsutil cat ~{vcf} | zcat | grep -m 1 CHROM | cut 10- > vcf.samples


    cat -<< "EOF" > check_inputs.R
        
    vcf_samples=names(read.csv("vcf.samples",sep="\t"))
    
    covariates_samples=names(read.csv("covariates.samples",sep="\t"))
    
    expression_samples=names(read.csv("expression.samples",sep='\t'))

    diff=setdiff(combined_samples, expression_samples)
    
    if(length(diff)){
      stop(paste("problem found. There are samples in combined_covariates that are not present in the expression data: ",paste(diff,collapse = " ")))
    }

    diff=setdiff(expression_samples, combined_samples)
    if(length(diff)){
      stop(paste("problem found. There are samples in the expression data that are not present in the combined_covariates: ",paste(diff,collapse = " ")))
    }

    diff=setdiff(expression_samples, vcf_samples)
    if(length(diff)){
      stop(paste("problem found. There are samples in the expression data that are not present in the vcf: ",paste(diff,collapse = ", ")))
    }

    EOF
    
    Rscript check_inputs.R

    touch _success
  >>>

  output {
    File success="_success"
  }
}


task fastqtl_permutations_scatter {
    input {
        File check
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

    call CheckInputs {
        input: 
            expression_bed=expression_bed,
            covariates=covariates,
            vcf=vcf    
    }

    call fastqtl_nominal {
        input:
            check=CheckInputs.success,
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
                check=CheckInputs.success,
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
        File? genes_nominal=fastqtl_nominal.allpairs

   }
}
