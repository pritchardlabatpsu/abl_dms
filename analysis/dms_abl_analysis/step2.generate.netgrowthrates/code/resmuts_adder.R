resmuts_adder=function(input_df){
  # This funciton adds a true or false flag of whether mutations in a dataframe fall amongst the clinically observed imatinib resistance mutations
  # Input: dataframe with the following columns: proteien_start, alt_aa
  # Output: dataframe with resmuts added
  input_df=input_df%>%
    mutate(resmuts=case_when(
      protein_start%in%248&alt_aa%in%"V"~T,
      protein_start%in%253&alt_aa%in%"H"~T,
      protein_start%in%255&alt_aa%in%"V"~T,
      protein_start%in%486&alt_aa%in%"S"~T,
      protein_start%in%396&alt_aa%in%"P"~T,
      protein_start%in%255&alt_aa%in%"K"~T,
      protein_start%in%315&alt_aa%in%"I"~T,
      protein_start%in%252&alt_aa%in%"H"~T,
      protein_start%in%253&alt_aa%in%"F"~T,
      protein_start%in%250&alt_aa%in%"E"~T,
      protein_start%in%359&alt_aa%in%"C"~T,
      protein_start%in%351&alt_aa%in%"T"~T,
      protein_start%in%355&alt_aa%in%"G"~T,
      protein_start%in%317&alt_aa%in%"L"~T,
      protein_start%in%359&alt_aa%in%"I"~T,
      protein_start%in%355&alt_aa%in%"A"~T,
      protein_start%in%459&alt_aa%in%"K"~T,
      protein_start%in%276&alt_aa%in%"G"~T,
      # protein_start%in%311&alt_aa%in%"L"~T,
      # protein_start%in%311&alt_aa%in%"I"~T,
      protein_start%in%244&alt_aa%in%"V"~T,
      T~F))
  input_df
}
