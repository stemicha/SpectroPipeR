#' ratio judgement function
#'
#' @param log2_peptide_ratio log2 peptide ratios
#' @param log2_median_protein_data_ratio log2 median protein ratios
#' @param log2_protein_ratios__peptide_ratios log2 protein to peptide ratios
#' @param log2_cut ratio cutoff (e.g. 1)
#' @return ratio judgement character output
#' @keywords internal
#' @noRd
#'
ratio_judgement_function <- function(log2_peptide_ratio = NULL,
                                     log2_median_protein_data_ratio = NULL,
                                     log2_protein_ratios__peptide_ratios = NULL,
                                     log2_cut = 1){

  if(sum(!is.na(log2_peptide_ratio), !is.na(log2_median_protein_data_ratio))==2){
    if( abs(log2_protein_ratios__peptide_ratios)>log2_cut){
      if(log2_peptide_ratio<0 & log2_median_protein_data_ratio<0){
        if(abs(log2_peptide_ratio)>abs(log2_median_protein_data_ratio)){
          return("protein ratio to low")
        }
        if(abs(log2_peptide_ratio)<abs(log2_median_protein_data_ratio)){
          return("protein ratio to high")
        }
      }
      if(log2_peptide_ratio>0 & log2_median_protein_data_ratio>0){
        if(abs(log2_peptide_ratio)>abs(log2_median_protein_data_ratio)){
          return("protein ratio to low")
        }
        if(abs(log2_peptide_ratio)<abs(log2_median_protein_data_ratio)){
          return("protein ratio to high")
        }
      }
      if(log2_peptide_ratio<=0 & log2_median_protein_data_ratio>0){
        return("protein ratio to high")
      }
      if(log2_peptide_ratio>=0 & log2_median_protein_data_ratio<0){
        return("protein ratio to low")
      }

    }else{
      return("OK")
    }


  }else{
    return(NA)
  }


}
