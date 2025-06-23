Files:

    BAM and BAI files for each biological replicate are named <Sample>.{bam,bai} 

    manifest.txt: The manifest for each sample in this study. The format is as follows:
        TwinstrandId : The internal TwinStrand ID.
        CustomerName : The ID provided by the customer.
        Annotation   : Human readable sample name.

    all.mut: The MUT files for each sample in the study. The format is as follows:
        Chromosome                      : hg38 Chromosome.
        Start                           : Base Start.
        End                             : Base Stop (unless indel this will just be the next base).
        Sample                          : TwinstrandId.
        VariationType                   : Type of variation, either, snv/snp, indel, or mnv.
        REF                             : Reference allele.
        ALT                             : Alternate allele.
        AltDepth                        : Alternate allele depth.
        Depth                           : Number of reads at base (exclusive of Ns).
        tki_resistant_mutation          : Variant is in database of known TKI-R mutations
        tki_resistant_mutation_evidence : Source of TKI-R mutation evidence 

    md5sum.txt: The md5sums of each BAM file.