version 1.0


import "../qtl/aFC.wdl" as aFC
import "Error.wdl" as E

workflow CrossSampleAFC {
	input { 

		Array[File] expression_beds
		Array[File] expression_bed_indexs
		Array[File] covariates_files
		
		File vcf_file
		File vcf_index

		Array[File] afc_qtl_files
		Array[String] prefixes

		String docker="gcr.io/broad-cga-francois-gtex/gtex_eqtl:V8"
		
		Int memory
		Int disk_space
		Int num_threads
		Int num_preempt

	}

	if ( length(expression_beds) != 2 ){
		call E.Error as error1 {input: message="Current version only accepts input of length 2", error=1}
	}


	if ( length(expression_bed_indexs) != length(expression_beds) ){
		call E.Error as error2 {input: message="length of bed index array was different from length of bed array", error=1}
	}

	if ( length(afc_qtl_files) != length(expression_beds) ){
		call E.Error as error3 {input: message="length of qtl array was different from length of bed array", error=1}
	}

	if ( length(prefixes) != length(expression_beds) ){
		call E.Error as error4 {input: message="length of prefixes array was different from length of bed array", error=1}
	}

	scatter (i in [1,2]){

		call aFC.convert_qtls as convert{
			input: 
				fastQTL_output=afc_qtl_files[i]
			}

		call aFC.aFC as acf_call {
			input:
				expression_bed=expression_beds[3-i],
				expression_bed_index=expression_bed_indexs[3-i],
				covariates_file=covariates_files[3-i],
				prefix=prefixes[i]+ "_in_" + prefixes[3-i], 
				afc_qtl_file=convert.qtl_file		
		}

		File afc_file=acf_call.afc_file

	}

	output {
		File one_in_two_afc=afc_file[1]
		File two_in_one_afc=afc_file[2]
	}
}