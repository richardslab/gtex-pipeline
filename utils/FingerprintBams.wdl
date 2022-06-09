version 1.0

import "../genotype/participant_vcfs.wdl" as participant_vcfs
import "../utils/Fingerprinting.wdl" as fp
import "../utils/IdentifySample.wdl" as id

workflow FingerprintBam{
	input {
		File hapMap
		File input_bam
		File input_bam_index

		File input_vcf
		String input_participant_id
		
	}	

	call participant_vcfs.participant_vcfs as get_het_vcfs{
		input:
			vcf_file=input_vcf,
     		participant_id=input_participant_id,

	    	memory=4,
    		disk_space=ceil(size(input_vcf,"GB")+10),
     		num_threads=1,
    		num_preempt=1
	} 


	call fp.CrossCheckSample as fingerprint {
		input:
		first=input_bam,
		first_index=input_bam_index,
		second=get_het_vcfs.snps_vcf,
		second_index=get_het_vcfs.snps_vcf_index,
		hapMap = hapMap	
	}

	call id.IdentifySampleWF as identifySample {
		input:
		sample=input_bam,
		sample_index=input_bam_index,
		hapMap=hapMap
	}

	output {
		File fingerprint_metrics=fingerprint.metrics
		File clustered_metrics=identifySample.fp_clustered
		File fingerprint_matrix=identifySample.fp_metrics
		String fingerprint_match=identifySample.match_group	
	}
}