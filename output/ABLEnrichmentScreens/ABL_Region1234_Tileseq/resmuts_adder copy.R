resmuts_adder=function(input_df){
  # This funciton adds a true or false flag of whether mutations in a dataframe fall amongst the clinically observed imatinib resistance mutations
  # Input: dataframe with the following columns: proteien_start, alt_aa
  # Output: dataframe with resmuts added
  input_df=input_df%>%
    mutate(resmuts=case_when(species%in%c("E255V",
                                          "Y253H",
                                          "T315I",
                                          "H396P",
                                          "H396R",
                                          "F486S",
                                          "G250E",
                                          "E255K",
                                          "Y253F",
                                          "Q252H",
                                          "L248V",
                                          "F359C",
                                          "F359I",
                                          "M351T",
                                          "F317L",
                                          "D276G",
                                          "M244V",
                                          "E459K")~T,
                             T~F))
  input_df
}
