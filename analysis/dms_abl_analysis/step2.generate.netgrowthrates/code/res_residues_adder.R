res_residues_adder=function(input_df){
  # library(dplyr)
  input_df=input_df%>%mutate(resresids=case_when(protein_start%in%253~T,
                                                 protein_start%in%255~T,
                                                 protein_start%in%486~T,
                                                 protein_start%in%396~T,
                                                 protein_start%in%255~T,
                                                 protein_start%in%315~T,
                                                 protein_start%in%252~T,
                                                 protein_start%in%253~T,
                                                 protein_start%in%250~T,
                                                 protein_start%in%359~T,
                                                 protein_start%in%351~T,
                                                 protein_start%in%355~T,
                                                 protein_start%in%317~T,
                                                 protein_start%in%359~T,
                                                 protein_start%in%355~T,
                                                 protein_start%in%459~T,
                                                 protein_start%in%276~T,
                                                 protein_start%in%299~T,

                                                 T~F))
  input_df
}
