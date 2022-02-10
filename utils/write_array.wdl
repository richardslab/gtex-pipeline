version 1.0

task write_array_to_tsv{
	input {
		Array[Array[String]] array
	}

	command <<<
	>>>

	output {
		File tsv=write_tsv(array)
	}
 	runtime {
        docker: "python:latest"
        memory: "2GB"
        disks: "local-disk 20 HDD"
    }
    meta {
        author: "Yossi Farjoun"
    }
}