version 1.0

task echo_task {
	input {
		String message
	}
	command <<<
		echo ${message}
	>>>

	runtime { 
		docker: "broadgdac/firecloud-ubuntu:15.10"
	}

	output {
		String echoed=read_string(stdout())
	}
}

workflow echo {
	call echo_task
}