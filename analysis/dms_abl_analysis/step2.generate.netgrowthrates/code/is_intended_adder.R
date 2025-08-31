is_intended_adder=function(input_df,codon_table=read.csv("data/codon_table.csv")){
  # 12.19.2024
  # This function adds a column is_intended to a table.
  # An intended codon is what twist made for our ABL library.
  # Inputs: input_df should have the following two columns: 1) protein_start, 2) alt_codon
  # input_df=screen_compare
  codon_table2=codon_table%>%filter(Twist2%in%T)
  codon_table=codon_table%>%filter(Twist%in%T)
  # Marking twist codons:
  input_df_default=input_df%>%
    mutate(is_intended=case_when(alt_codon%in%codon_table$Codon~1,
                                 T~0))

  # This next part says that if you're at residue <494, use codons from table 2 instead.
  if(max(input_df$protein_start)>=495){
  input_df_default=input_df%>%
    mutate(is_intended=case_when(protein_start<=494&alt_codon%in%codon_table$Codon~1,
                                 protein_start>=495&alt_codon%in%codon_table2$Codon~1,
                                 T~0))
  }
  input_df=input_df_default
  input_df
}
