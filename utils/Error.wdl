version 1.0

task Error{
	input {
		String message
		Int error
	}

	command <<<
		echo '~{message}'

		exit ~{error}
	>>>

 	runtime {
        docker: "ubuntu:latest"
        memory: "1GB"
        disks: "local-disk 20 HDD"
        cpu: 1
    }
    meta {
        author: "Yossi Farjoun"
    }
}