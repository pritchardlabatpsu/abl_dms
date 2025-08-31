# inputs:
# codon_table
# ref_codon="AAT"
# alt_aa="A"

#outputs:
#shortest codon nucleotides
#shortest alt codn
# codon_table=read.csv("data/codon_table.csv",header = T,stringsAsFactors = F)
# For a given reference codon and a given ALT residue, this function will find the shortest possible ALT codon. If multiple codons are available
shortest_codon_finder=function(ref_codon,alt_aa,codon_table=read.csv("data/codon_table.csv",header = T,stringsAsFactors = F)){


  #What is the twist codon for the given alt amino acid?
  # twist_codons=codon_table%>%filter(Twist%in%T)
  # twist_alt_codon=twist_codons[twist_codons$Letter%in%alt_aa,"Codon"]

  #Calculating hamming distance
  #Below, I calculate hamming distance by looking for a mismatch between the REF codon and the alt codon at each one of the 3 positions of the codon and adding them up.
  alt_codons=codon_table%>%
    filter(Letter%in%alt_aa)%>%
    rowwise()%>%
    mutate(distance=as.numeric(!substr(ref_codon,1,1)%in%substr(Codon,1,1))+
             as.numeric(!substr(ref_codon,2,2)%in%substr(Codon,2,2))+
             as.numeric(!substr(ref_codon,3,3)%in%substr(Codon,3,3)))%>%
    ungroup()


  #Sorting through the shortest hamming distance
  mindist=min(alt_codons$distance)
  alt_codons=alt_codons%>%filter(distance==mindist)

  #If there are multiple rows at the shortest distance, and one of the codons is the twist codon, choose the twist codon. Aka if no codon is shorter than twist, choose twist
  #Otherwise if there are multiple shortest codons and none of them are the twist codon, then return both of them
  if(sum(as.numeric(alt_codons$Twist%in%T))!=0){
    #The above boolean is saying "if you find twist codon amongst the shortest codons"
    shortest_alt_codon=as.character(alt_codons[alt_codons$Twist%in%T,"Codon"])
    # shortest_alt_codon=alt_codons$Codon
    # class(shortest_alt_codon)
  } else{
    shortest_alt_codon=alt_codons$Codon[1]
    # Note that when there are multiple options for non-twist shortest ALT codons that code for the same residue, I just pick one of them. This happens for the following:
    # "AGG,CGG"
    # "CTC,TTA,TTG"
    # "CTA,TTA"
    # "AGT,TCT"
    # To view all shortest alt codons, use:
    # shortest_alt_codon=paste(alt_codons$Codon,collapse = ",")
  }
  # return(list(shortest_alt_codon))
  return(list(shortest_alt_codon,mindist))
}
