

## eQTL and sQTL pipeline overview.


- (per sample) fastq files are aligned with sqtl/sQTL.wdl which aligns (2-pass star) with WASP correction and fingerprints against the gwas and WGS data.

### Expression
- (per sample) sequencing data is DuplicateMarked and expression evaluated (per sample) by rnaseqc2 in the aligned_rnaseq_bam_to_genecounts.wdl pipeline 
- (jointly) expression is combined via rnaseqc2_aggregate.wdl (CombineRNASeqc2Files in Terra) 
- (jointly) expression is normalized via eqtl_prepare_expression.wdl (PrepareExpressionData in Terra)


### Splicing
- (per sample) junctions are discovered and counted via leafcutter_bam_to_junc.wdl which runs on the non-deduped BAM.
- (jointly) junction files are combined with leafcutter_cluster.wdl (ClusterJunctions in Terra)
- (jointly) get the junction->gene mapping file with GetGeneId


### Per cohort/arm
- (jointly*) covariates (Sex, PCs) are prepared via PrepareCovariates.wdl
- (jointly*) extra covariates (age, extraction site, and sequencing run) are also prepared via PrepareCovariates.wdl, but with a different configuration.
-

### per QTL : cohort/arm
- (jointly) PEER factors are calculated via qtl_peer_factors.wdl (expressionQTL_PEERFactors and spliceQTL_PEERFactors in terra)
- (jointly) QTLs found via fastqtl.wdl (FastQTL and Splice_FastQTL in Terra)
- (jointly) annotate and filter QTLs vis fastqtl_annotate.wdl

	