cosmic_data_adder=function(inputdf){
  #The following function adds cosmic mutational status and counts to mutation data.
  #Inputs: An input dataframe with a column called "species"
  #Output: Same as input dataframe with cosmic count and cosmic mutational status
  #Reading and parsing cosmic data
  cosmic_data=read.table("data/Cosmic_ABL/ABL_Cosmic_Gene_mutations.tsv",sep="\t",header = T,stringsAsFactors = )
  cosmic_data=cosmic_data%>%mutate(AA.Mutation=gsub("p.","",AA.Mutation))
  cosmic_data=cosmic_data[!grepl("ins|del",cosmic_data$CDS.Mutation),]
  cosmic_data=cosmic_data[grepl("Missense",cosmic_data$Type),]
  cosmic_data=cosmic_data%>%filter(!Type%in%"Substitution - coding silent")
  # cosmic_data=cosmic_data%>%filter(Position<=500,Position>=64,Count>=2)
  cosmic_data=cosmic_data%>%filter(Position<=500,Position>=64)
  # write.csv(cosmic_data,"cosmic_abl.csv")
  cosmic_simple=cosmic_data%>%dplyr::select(subs_name=AA.Mutation,cosmic_count=Count)
  cosmic_simple=cosmic_simple%>%group_by(subs_name)%>%summarize(cosmic_count=sum(cosmic_count))
  cosmic_simple$cosmic_present=T
  #Merging input df and cosmic data
  data_cosmic=merge(inputdf,cosmic_simple,by.x="species",by.y = "subs_name",all.x = T)
  data_cosmic[data_cosmic$cosmic_present%in%NA,"cosmic_present"]=F
  data_cosmic[data_cosmic$cosmic_count%in%NA,"cosmic_count"]=0
  data_cosmic
}
