version 1.0

task GetGeneIdTask {
    input {
        File splice_bed
        String prefix
    }
    Int disk_size = 10 + ceil(size([splice_bed], "GB"))

    command <<<
        zcat "~{splice_bed}" | cut -f 4 | awk -F ':' 'NR>1{print $0"\t"$5}' > "~{prefix}.splice_to_gene.txt"
    >>>
    output {
        File map_out="~{prefix}.splice_to_gene.txt"
    }

    
    runtime {
        docker: "ubuntu:latest"
        memory: "1GB"
        disks: "local-disk " + disk_size + " HDD"
        cpu: 1
    }
}


workflow GetGeneId {
    call GetGeneIdTask{}

    output {
        File map=GetGeneIdTask.map_out
    }
}