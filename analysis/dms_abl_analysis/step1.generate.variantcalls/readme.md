The variant caller provided in this directory was used for variant calling on consensus-called single-strand consensus sequencing files.

Inputs:
The input files for this variant caller are present in the following directories:
dms_abl_analysis/step2.generate.netgrowthrates/data/Consensus_Data/Novogene_lane18/sample2/l298l/sscs
	Sample 2 is our baseline day 6 library prep of genomic DNA extracted from BaF3s transduced with the BCRABL library after 6 days of il3 independent growth.
dms_abl_analysis/step2.generate.netgrowthrates/data/Consensus_Data/Novogene_lane18/sample3/l298l/sscs
	Sample 3 is our baseline day 6 library prep of genomic DNA extracted from BaF3s transduced with the BCRABL library after 6 days of il3 independent growth.
dms_abl_analysis/step2.generate.netgrowthrates/data/Consensus_Data/Novogene_lane18/sample10/l298l/sscs
	Sample 10 is our baseline day 0 library prep of genomic DNA extracted from BaF3s transduced with the BCRABL library.


Outputs:
The output files from our variant calling are present in the following direcoties:
dms_able_analysis/step2.generate.netgrowthrates/data/Consensus_Data/Novogene_lane18/sample2/sscs/variant_caller_outputs
dms_able_analysis/step2.generate.netgrowthrates/data/Consensus_Data/Novogene_lane18/sample3/sscs/variant_caller_outputs
dms_able_analysis/step2.generate.netgrowthrates/data/Consensus_Data/Novogene_lane18/sample10/sscs/variant_caller_outputs


Informatics to generate consensus-called reads:
This variant caller requires single-strand consensus read files. These are similar to the .BAM outputs of bwa-mem2.
A custom pipeline, described below, was used to generate these BAM files. This description is also present in the methods section of the manuscript.

One-hundred-fifty-nucleotide paired-end sequencing of the UMIs and mutagenized region was
done on an Illumina NovaSeq 6000. For each sample, dunovo was used to generate error-corrected
single strand consensus from the UMI barcodes. Then bwa-mem2 was used to align the census to
human ABL cDNA. After filtering out for mouse ABL reads, aligned reads with less than 5 mismatches
would undergo variant calling and annotation using a custom R script. Briefly, for each alignment,
variants were converted from the MDZ read tag and normalized to read depth at that position. Mutant
growth rates were calculated using exponential growth equation and the mutant allele frequency (MAF).
