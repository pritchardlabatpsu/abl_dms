merge_samples_twosamples=function(inputx,inputy){
  # Note: there is a newer version of this function called merge_samples.R
  # The only reason I kept this functino is because it can handle dataframes at inputs.
  #This function takes variant from two samples, merges them, and adds counts and allele frequencies
  # inputx="Novogene_lane11/sample1/duplex"
  # inputy="Novogene_lane11/sample2/sscs"
  # inputx="Novogene_lane14/Sample10_combined/sscs" #Input the directory names for the sequencing Lane and the sample name
  # inputy="Novogene_lane15/sample_3/sscs" #Input the directory names for the sequencing Lane and the sample name
  samplex=inputx
  if(class(inputx)=="character"){ #aka saying that if inputx is not already present as a dataframe
    samplex=read.csv(paste("data/Consensus_Data/",inputx,"/variant_caller_outputs/variants_unique_ann.csv",sep=""),header=T,stringsAsFactors = F)
  }

  sampley=read.csv(paste("data/Consensus_Data/",inputy,"/variant_caller_outputs/variants_unique_ann.csv",sep=""),header=T,stringsAsFactors = F)

  samplex_depths=samplex%>%group_by(protein_start)%>%summarize(depth.x=mean(depth))
  samplex=samplex%>%dplyr::select(-depth)
  sampley_depths=sampley%>%group_by(protein_start)%>%summarize(depth.y=mean(depth))
  sampley=sampley%>%dplyr::select(-depth)
  samples_xy=merge(samplex,sampley,by=c("type",
                                        "alt_start_pos",
                                        "alt_end_pos",
                                        "ref",
                                        "ref_codon",
                                        "alt",
                                        "alt_codon",
                                        "frame_pos",
                                        "protein_start",
                                        "protein_end",
                                        "ref_aa",
                                        "alt_aa",
                                        "amino_acids",
                                        "consequence_terms"),all=T)

  samples_xy=merge(samples_xy,samplex_depths,by="protein_start")
  samples_xy=merge(samples_xy,sampley_depths,by="protein_start")
  #If a mutant isn't seen in a sample, it's count is 0
  samples_xy[samples_xy$ct.x%in%NA,"ct.x"]=0
  samples_xy[samples_xy$ct.y%in%NA,"ct.y"]=0
  #If a mutant isn't seen in a sample, the depth at that site is the average depth seen at that residue
  samples_xy[samples_xy$depth.x%in%NA,"depth.x"]=0
  samples_xy[samples_xy$depth.y%in%NA,"depth.y"]=0
  samples_xy=samples_xy%>%
    mutate(ct=ct.x+ct.y,depth=depth.x+depth.y)%>%
    dplyr::select("type",
                  "alt_start_pos",
                  "alt_end_pos",
                  "ref",
                  "ref_codon",
                  "alt",
                  "alt_codon",
                  "frame_pos",
                  "ct",
                  "depth",
                  "protein_start",
                  "protein_end",
                  "ref_aa",
                  "alt_aa",
                  "amino_acids",
                  "consequence_terms")

  samples_xy
}
