version 1.0 


workflow CreateSampleParticipantMap{
	input {
		Array[String] samples
		Array[String] participants
	}

	Array[Array[String]] sample_participant_array =  transpose([samples,participants])
	call write_array_to_tsv {
		input:
			array=sample_participant_array
	}

	output {
		File map=write_array_to_tsv.tsv
	}	
}


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