version 1.0


import "../qtl/aFC.wdl" as aFC
import "Error.wdl" as E

workflow CrossSampleAFC {
	input { 

		Array[File] expression_beds
		Array[File] expression_bed_indexs
		Array[File] covariates_files
		
		File vcf_file
		File vcf_file_index

		Array[File] afc_qtl_files
		Array[String] prefixes

		String docker="gcr.io/broad-cga-francois-gtex/gtex_eqtl:V8"
		

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

	scatter (this in [0,1]) {
		Int other=2-this

		call aFC.convert_qtls as convert{
			input: 
				prefix=prefixes[this],
				fastQTL_output=afc_qtl_files[this]
			}


		call aFC.aFC as acf_call {
			input:
				expression_bed=expression_beds[other],
				expression_bed_index=expression_bed_indexs[other],
				covariates_file=covariates_files[other],
				prefix=prefixes[this]+ "_in_" + prefixes[other], 
				afc_qtl_file=convert.qtl_file,
				docker=docker,
				vcf_file=vcf_file,
				vcf_index=vcf_file_index
		}

		File afc_file=acf_call.afc_file

	}

	output {
		File one_in_two_afc=afc_file[0]
		File two_in_one_afc=afc_file[1]
	}
}