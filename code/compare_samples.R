compare_samples=function(input_df_x,
                         input_df_y,
                         netgr_wodrug=0.055, #In inverse hours,
                         duration=48){  #In hours

  # The following function compares before and after samples and calculates enrichment scores
  # input_df_x=before_timepoint
  # input_df_y=after_timepoint
  # netgr_wodrug=netgr
  # duration=delta_t

  input_df_xy=merge(input_df_x%>%filter(consequence_terms%in%"missense_variant"),input_df_y%>%filter(consequence_terms%in%"missense_variant"),by=c("ref_aa","protein_start","alt_aa","ref","alt","alt_start_pos","consequence_terms","ref_codon","alt_codon","alt_codon_shortest","n_nuc_min"),all.x = T)
  input_df_xy$alt_aa=factor(input_df_xy$alt_aa,levels=c("P","G","Y","W","F","V","L","I","A","T","S","Q","N","M","C","E","D","R","K","H"))
  input_df_xy[input_df_xy$ct.y%in%NA,"ct.y"]=.5

  input_df_xy=depth_finder(input_df_xy,"depth.y")
  input_df_xy=depth_finder(input_df_xy,"totalcells.y")
  input_df_xy=input_df_xy%>%mutate(maf.y=ct.y/depth.y,
                                   totalmutant.y=maf.y*totalcells.y)

  input_df_xy[input_df_xy$ct.x%in%NA,"ct.x"]=.5

  input_df_xy=input_df_xy%>%mutate(score=log2(maf.y/maf.x))
  # il3D0.D2[il3D0.D2$score%in%NA,"score"]=-6

  input_df_xy[input_df_xy$totalmutant.y%in%NA,"totalmutant.y"]=0
  input_df_xy=input_df_xy%>%mutate(netgr_obs=log(totalmutant.y/totalmutant.x)/duration)
  input_df_xy[input_df_xy$netgr_obs%in%NA,"netgr_obs"]=-netgr_wodrug
  input_df_xy
}
