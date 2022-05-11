version 1.0 

import "write_array.wdl" as wa

workflow CreateSampleParticipantMap{
	input {
		Array[String] samples
		Array[String] participants

		Array[String] header=["sample","participant"]
	}

	Array[Array[String]] sample_participant_array = flatten([[header], transpose([samples,participants])])

	call wa.write_array_to_tsv {
		input:
			array=sample_participant_array
	}

	output {
		File map=write_array_to_tsv.tsv
	}	
}


