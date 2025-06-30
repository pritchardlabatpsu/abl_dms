# The following function takes in output from the variant caller, and filters it
# Filters for twist codon
# Calculates MAFs
# Calculates total mutant cells/mL based on total cells/mL and MAF. The total cells value defaults to 1.
# Parses out the ref amino acid and the alt amino acid
# Figures out which codons would realistically be the shortest codons (based on hamming distance)
# Also filters for twist's 19 expected codons according to twist's (as opposed to all possible codons)
# Update 07.28.24:
# Twist used one set of codons for mutagenizing ABL residues 242:494 and another set for 495:512
# This function was updated so that it counts a unique set of alt codons residues 495 onwards
# Note that this update assumes that the updated twist alternate codons only occur during residues 495 through 512.
# This change also affects the heatmap_plotting_function code
#
# input_df=before_timepoint
# totalcells=cells_before
# a=variants_parser(before_timepoint,cells_before)

variants_parser=function(input_df,totalcells=1,codon_table=read.csv("data/codon_table.csv")){
  # input_df=after_timepoint
  # totalcells=cells_after
  # tic()
  # codon_table=read.csv("data/codon_table.csv")
  codon_table2=codon_table%>%filter(Twist2%in%T)
  codon_table=codon_table%>%filter(Twist%in%T)

  # Filtering for twist codons:
  # input_df_default=input_df%>%filter(protein_start<=494,alt_codon%in%codon_table$Codon)
  #
  # if(max(input_df$protein_start)>=495){
  #   input_df_extended=input_df%>%filter(protein_start>=495,alt_codon%in%codon_table2$Codon)
  #   input_df_nonextended=input_df%>%filter(protein_start<=494,alt_codon%in%codon_table$Codon)
  #   input_df_default=rbind(input_df_nonextended,input_df_extended)
  # }
  # # a=input_df_extended%>%filter(protein_start%in%497,alt_aa%in%c("A","D","G","P"))
  # input_df=input_df_default



  ##############################Not parrallel option for variants parser###########
  # input_df=input_df%>%
  #   rowwise()%>%
  #   mutate(alt_codon_shortest=shortest_codon_finder(ref_codon,alt_aa)[[1]][1],
  #          n_nuc_min=shortest_codon_finder(ref_codon,alt_aa)[[2]][1])
  ############################################################
  ##############################Parrallel option for variants parser###########
  library(foreach)
  library(doParallel)
  library(tictoc)
  cores=detectCores()
  cl= makeCluster(cores[1]-1)
  registerDoParallel(cl)

  # i=1
  alt_codon_shortest=foreach(i=seq(1,nrow(input_df)),combine=rbind,.packages = c("dplyr",
                                                                       "httr",
                                                                       "xml2")) %dopar% {
                                                                         source("code/shortest_codon_finder.R")
                                                                         # input_df$alt_codon_shortest[i]=shortest_codon_finder(input_df[i,"ref_codon"],input_df[i,"alt_aa"])[[1]][1]
                                                                         # input_df$n_nuc_min[i]=shortest_codon_finder(input_df[i,"ref_codon"],input_df[i,"alt_aa"])[[2]][1]
                                                                         # input_df
                                                                         alt_codon_shortest=shortest_codon_finder(input_df[i,"ref_codon"],input_df[i,"alt_aa"])
                                                                         alt_codon_shortest
                                                                       }

  stopCluster(cl)
  alt_codon_shortest=as.data.frame(do.call(rbind,alt_codon_shortest))
  input_df$alt_codon_shortest=alt_codon_shortest$V1
  input_df$n_nuc_min=alt_codon_shortest$V2
  ############################################################

  input_df=input_df%>%mutate(maf=ct/depth)
  input_df$totalcells=totalcells
  input_df=input_df%>%mutate(totalmutant=maf*totalcells)
  input_df_simple=input_df%>%dplyr::select(alt_start_pos,protein_start,ref,alt,ref_aa,alt_aa,ref_codon,alt_codon,alt_codon_shortest,n_nuc_min,consequence_terms,ct,depth,maf,totalcells,totalmutant)
  input_df_simple
  # toc()
}
