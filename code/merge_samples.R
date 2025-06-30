merge_samples <- function(...){
  #This function takes variant from two samples, merges them, and adds counts and allele frequencies
  ###########Difference between this function and merge_samples_twosamples.R#########
  # This function is the exact same as merge_samples_twosamples.R
  # but it is better because it can handle more than two samples.
  # i.e. when, for sample, imatinib d0, was prepped using 3 separate lirbary preps on 3 different sequencing lanes
  # The only reason I kept the original merge samples in the code repository is because it can
  # handle a use case where the user enters a dataframe in sample 1, and not just a path to the csv.
  # As of right now, this function can not handle dataframe inputs, it can only handle locations of directories.
  # if(class(inputx)=="character"){ #aka saying that if inputx is not already present as a dataframe
  #   samplex=read.csv(paste("data/Consensus_Data/",inputx,"/variant_caller_outputs/variants_unique_ann.csv",sep=""),header=T,stringsAsFactors = F)
  # }
  #############################################
  ##################Inputs###########################
  #Input the directory names for the sequencing Lane and the sample name
  # samples_xy=merge_samples(paste(before_screen1_counts$dirname_input[1],"/sscs",sep=""),
  #                          paste(before_screen1_counts$dirname_input[2],"/sscs",sep=""),
  #                          paste(before_screen1_counts$dirname_input[2],"/sscs",sep=""))
  # Each of these looks like something this: "novogene_lane19_20/Ln19a_Sample11/sscs"

  # Get the list of file paths
  # file_paths=filepaths #feeds as an input from compare_screens.R
  # file_paths=list(x)
  # file_paths=list(paste(before_screen1_counts$dirname_input[1],"/sscs",sep=""),
  #                 paste(before_screen1_counts$dirname_input[2],"/sscs",sep=""))
  file_paths <- list(...)
  # Reformat file paths so that they represent paths to directories
  file_paths=paste("data/Consensus_Data/",file_paths,"/variant_caller_outputs/variants_unique_ann.csv",sep="")

  # Read each CSV file into a list of dataframes
  dataframes <- lapply(file_paths, read.csv)
  # samplex=dataframes[[1]]
  # sampley=dataframes[[2]]
  # Function to merge two dataframes by 'type' and 'protein_start', etc
  merge_function <- function(x, y){
    samplex=x
    sampley=y
    # merge(x, y, by = c("type", "protein_start"), all = TRUE)
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

    samples_xy=merge(samples_xy,samplex_depths,by="protein_start",all = T)
    samples_xy=merge(samples_xy,sampley_depths,by="protein_start",all = T)
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

  # Iteratively merge all dataframes
  merged_df <- Reduce(merge_function, dataframes)
  # Return the merged dataframe
  return(merged_df)
}

