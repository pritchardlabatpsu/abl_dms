depth_finder=function(input_df,colname){
  library(tidyr)
  #This function returns the mean depth at a desired residue
  #Mean depth is the average of the bwa-mem depth of all mutants seen at that residue
  # input_df=input_df_xy
  # colname="depth.y"
  # input_df=il3D0.D2
  # resi=c(242:250)
  # colname="depth.y"
  input_df$coverage2=input_df[[colname]]
  depths=input_df%>%
    group_by(protein_start)%>%
    summarize(depth_mean=mean(coverage2,na.rm=T))%>%
    filter(!protein_start%in%NA)

    # The following piece of code finds any residue at which the depth is NaN and replcaes it with the depth that was seen at the previous non NaN residue.
    # The depth at a residue would be NaN if no mutants were seen at that residue. This can be a problem for low coverage samples.
    depths <- depths %>%
    mutate(depth_mean = ifelse(is.nan(depth_mean), NA, depth_mean)) %>%
    fill(depth_mean, .direction = "down")

  # x=input_df%>%filter(protein_start%in%"321")
  # protein_start=242
  # depths[depths$protein_start==protein_start,"depth_mean"][[1]]
  input_df$coverage2=as.numeric(input_df$coverage2)
  input_df=input_df%>%
    rowwise()%>%
    mutate(coverage2=case_when(coverage2%in%NA~
                                 depths[depths$protein_start==protein_start,"depth_mean"][[1]],
                                                 T~coverage2))

  # class(input_df$coverage2)
  input_df[,eval(colname)]=input_df$coverage2
  input_df%>%dplyr::select(-coverage2)
}

