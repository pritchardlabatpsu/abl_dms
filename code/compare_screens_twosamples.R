compare_screens_twosamples=function(before_screen1_identifier,
                         after_screen1_identifier,
                         before_screen2_identifier,
                         after_screen2_identifier,
                         cellcounts_matrix_location){

  # before_screen1_identifier=before_screen1_identifier
  # after_screen1_identifier=after_screen1_identifier
  # before_screen2_identifier=before_screen2_identifier
  # after_screen2_identifier=after_screen2_identifier
  # cellcounts_matrix_location=cellcounts_matrix_dir


  # Inputs: 1) Condition 1 baseline, 2) Condition 1 treated, 3) Condition 2 Baseline, 4) Condition 2 treated
  # 5)The location of the cell counts table
  # These inputs can have a single sample in the name or specify two samples to merge
  # This function takes in two screens and makes comparisons across them.
  # This function uses the compare samples function, that compares the before treatment and the treated timepoints of a screen
  # This function compares two screens
  # Inputs:
  # Which conditions would you like to compare:
  # before_screen1_identifier=comparisons$before_screen1_identifier[i]
  # after_screen1_identifier=comparisons$after_screen1_identifier[i]
  # before_screen2_identifier=comparisons$before_screen2_identifier[i]
  # after_screen2_identifier=comparisons$after_screen2_identifier[i]
  # cellcounts_matrix_location=cellcounts_matrix_dir

  # before_screen1_identifier="KTRMD0"
  # before_screen1_identifier=c("KTRLD0","KTRLD0_FT")
  # after_screen1_identifier=c("KTRL1_D6","KTRL2_D6")
  # after_screen1_identifier="KTRM1_D6"
  # before_screen2_identifier="KTRMD0"
  # after_screen2_identifier=c("KTRM1_D6","KTRM2_D6")
  # after_screen2_identifier="KTRM2_D6"


  # before_screen1_identifier="BTR3ID0"
  # after_screen1_identifier="BTR3M_D6_A"
  # before_screen2_identifier="BTR3ID0"
  # after_screen2_identifier="BTR3M_D6_B"
  # cellcounts_matrix_location="data/Consensus_Data/novogene_lane23/TwistAllRegions_CellCounts_Matrix.csv"
  source("code/res_residues_adder.R")
  source("code/resmuts_adder.R")
  source("code/merge_samples_twosamples.R")
  source("code/variants_parser.R")
  source("code/compare_samples.R")
  source("code/depth_finder.R")
  # after_screen1_counts=cell_counts_table%>%filter(Identifier%in%after_screen1_identifier)

  # The cell counts and the exact deltat is derived from a spreadsheet in the novogene lane 18 folder
  # cellcounts_matrix_location=cellcounts_matrix_dir
  cell_counts_table=read.csv(cellcounts_matrix_location)
  cell_counts_table=cell_counts_table%>%filter(!skip_flag%in%"TRUE")

  before_screen1_counts=cell_counts_table%>%filter(Identifier%in%before_screen1_identifier)
  after_screen1_counts=cell_counts_table%>%filter(Identifier%in%after_screen1_identifier)
  before_screen2_counts=cell_counts_table%>%filter(Identifier%in%before_screen2_identifier)
  after_screen2_counts=cell_counts_table%>%filter(Identifier%in%after_screen2_identifier)


  delta_t=mean(after_screen1_counts$Time_Hours)-mean(before_screen1_counts$Time_Hours) #hours
  cells_before=mean(before_screen1_counts$TotalCells) #total cells at before time point
  cells_after=mean(after_screen1_counts$TotalCells) #total cells at after time point
  netgr=mean(after_screen1_counts$Netgr_wodrug)


  # Running merge samples on all the directories specified for this sample
  # This if statements sees if there are two samples specified at a timepoint, and if so, merges those samples
  # Note that the code below can only do at max two separate sequencing samples. It can be modified in the future to allow multiple samples
  if(length(before_screen1_counts[,1])>=2){
    before_timepoint=merge_samples_twosamples(paste(before_screen1_counts$dirname_input[1],"/sscs",sep=""),paste(before_screen1_counts$dirname_input[2],"/sscs",sep=""))
  } else {
    before_timepoint=read.csv(paste("data/Consensus_Data/",before_screen1_counts$dirname_input[1],"/sscs/variant_caller_outputs/variants_unique_ann.csv",sep=""))
  }
  # filepaths=paste(before_screen1_counts$dirname_input,"/sscs",sep="")
  # before_timepoint=do.call(merge_samples,as.list(filepaths))

  before_timepoint=variants_parser(before_timepoint,cells_before)

  # Running merge samples on all the directories specified for this sample
  # This if statements sees if there are two samples specified at a timepoint, and if so, merges those samples
  if(length(after_screen1_counts[,1])>=2){
    after_timepoint=merge_samples_twosamples(paste(after_screen1_counts$dirname_input[1],"/sscs",sep=""),paste(after_screen1_counts$dirname_input[2],"/sscs",sep=""))
  } else {
    after_timepoint=read.csv(paste("data/Consensus_Data/",after_screen1_counts$dirname_input[1],"/sscs/variant_caller_outputs/variants_unique_ann.csv",sep=""))
  }
  # filepaths=paste(after_screen1_counts$dirname_input,"/sscs",sep="")
  # after_timepoint=do.call(merge_samples,as.list(filepaths))


  after_timepoint=variants_parser(after_timepoint,cells_after)

  before_after=compare_samples(before_timepoint,after_timepoint,netgr,delta_t)

  # before_after=before_after%>%filter(ct.x>=3)
  before_after_screen1=before_after


  ###### TSII Day 0 vs Day 4 #########





  delta_t=mean(after_screen2_counts$Time_Hours)-mean(before_screen2_counts$Time_Hours) #hours
  cells_before=mean(before_screen2_counts$TotalCells) #total cells at before time point
  cells_after=mean(after_screen2_counts$TotalCells) #total cells at after time point
  netgr=mean(after_screen2_counts$Netgr_wodrug)


  # Running merge samples on all the directories specified for this sample
  # This if statements sees if there are two samples specified at a timepoint, and if so, merges those samples
  # if(length(before_screen2_counts[,1])>=2){
  #   before_timepoint=merge_samples(paste(before_screen2_counts$dirname_input[1],"/sscs",sep=""),paste(before_screen2_counts$dirname_input[2],"/sscs",sep=""))
  # } else {
  #   before_timepoint=read.csv(paste("data/Consensus_Data/",before_screen2_counts$dirname_input[1],"/sscs/variant_caller_outputs/variants_unique_ann.csv",sep=""))
  # }
  filepaths=paste(before_screen2_counts$dirname_input,"/sscs",sep="")
  before_timepoint=do.call(merge_samples,as.list(filepaths))

  before_timepoint=variants_parser(before_timepoint,cells_before)

  # Running merge samples on all the directories specified for this sample
  # This if statements sees if there are two samples specified at a timepoint, and if so, merges those samples
  # if(length(after_screen2_counts[,1])>=2){
  #   after_timepoint=merge_samples(paste(after_screen2_counts$dirname_input[1],"/sscs",sep=""),paste(after_screen2_counts$dirname_input[2],"/sscs",sep=""))
  # } else {
  #   after_timepoint=read.csv(paste("data/Consensus_Data/",after_screen2_counts$dirname_input[1],"/sscs/variant_caller_outputs/variants_unique_ann.csv",sep=""))
  # }
  filepaths=paste(after_screen2_counts$dirname_input,"/sscs",sep="")
  after_timepoint=do.call(merge_samples,as.list(filepaths))

  after_timepoint=variants_parser(after_timepoint,cells_after)

  before_after=compare_samples(before_timepoint,after_timepoint,netgr,delta_t)

  # before_after=before_after%>%filter(ct.x>=3)
  before_after_screen2=before_after



  ############### Adding Count, Depth, and MAF columns for Ivan 2.20.23###########
  # a=before_after_screen1
  # b=before_after_screen2
  # before_after_screen1=a
  # before_after_screen2=b

  before_after_screen1=before_after_screen1%>%
    mutate(ct_screen1_before=ct.x,
           depth_screen1_before=depth.x,
           maf_screen1_before=maf.x,
           ct_screen1_after=ct.y,
           depth_screen1_after=depth.y,
           maf_screen1_after=maf.y)%>%
    dplyr::select(-c("ct.x","depth.x","maf.x","totalcells.x","totalmutant.x","ct.y","depth.y","maf.y","totalcells.y","totalmutant.y"))

  before_after_screen2=before_after_screen2%>%
    mutate(ct_screen2_before=ct.x,
           depth_screen2_before=depth.x,
           maf_screen2_before=maf.x,
           ct_screen2_after=ct.y,
           depth_screen2_after=depth.y,
           maf_screen2_after=maf.y)%>%
    dplyr::select(-c("ct.x","depth.x","maf.x","totalcells.x","totalmutant.x","ct.y","depth.y","maf.y","totalcells.y","totalmutant.y"))

  screen_compare=merge(before_after_screen1,
                       before_after_screen2,
                       by=c("ref_aa","protein_start","alt_aa","ref","alt","alt_start_pos","consequence_terms","ref_codon","alt_codon","alt_codon_shortest","n_nuc_min"),
                       all.x=T,
                       all.y=T,suffixes = c("_screen1","_screen2"))

  screen_compare=resmuts_adder(screen_compare)
  screen_compare=res_residues_adder(screen_compare)
  screen_compare=screen_compare%>%mutate(species=paste(ref_aa,protein_start,alt_aa,sep = ""))

  screen_compare_means=screen_compare%>%
    filter(!score_screen1%in%NA,!score_screen2%in%NA,!score_screen1%in%NaN,!score_screen2%in%NaN)%>%
    rowwise()%>%
    mutate(score_mean=mean(c(score_screen1,score_screen2)),
           netgr_obs_mean=mean(c(netgr_obs_screen1,netgr_obs_screen2)))

  # ###8/30/23 Change for Ivan (this update avoids remove NA counts on D0):
  screen_compare_means=screen_compare%>%
    # filter(!score_screen1%in%NA,!score_screen2%in%NA,!score_screen1%in%NaN,!score_screen2%in%NaN)%>%
    rowwise()%>%
    mutate(score_mean=mean(c(score_screen1,score_screen2)),
           netgr_obs_mean=mean(c(netgr_obs_screen1,netgr_obs_screen2)))

  # write.csv(screen_compare_means,"Imat_Enrichment_bgmerged_2.22.23.csv")
  ##############

  # before_after_screen1=before_after_screen1%>%
  #   mutate(ct_screen1=ct.x)%>%
  #   dplyr::select(-c("ct.x","depth.x","maf.x","totalcells.x","totalmutant.x","ct.y","depth.y","maf.y","totalcells.y","totalmutant.y"))
  #
  #
  # before_after_screen2=before_after_screen2%>%
  #   mutate(ct_screen2=ct.x)%>%
  #   dplyr::select(-c("ct.x","depth.x","maf.x","totalcells.x","totalmutant.x","ct.y","depth.y","maf.y","totalcells.y","totalmutant.y"))
  #
  # screen_compare=merge(before_after_screen1,
  #                      before_after_screen2,
  #                      by=c("ref_aa","protein_start","alt_aa","ref","alt","alt_start_pos","consequence_terms","ref_codon","alt_codon","alt_codon_shortest","n_nuc_min"),
  #                      all.x=T,
  #                      all.y=T)
  #
  # screen_compare=resmuts_adder(screen_compare)
  # screen_compare=screen_compare%>%mutate(species=paste(ref_aa,protein_start,alt_aa,sep = ""))
  screen_compare_means
}
