#' effect size estimation function
#'
#' @param data effect size input
#' @param condition_comparisons condition comparison
#' @return effect size estimation
#' @keywords internal
#' @noRd
#'
#' @import tidyverse
#' @import effsize
effect_size_estimation <- function(data = NULL,
                                   condition_comparisons = NULL){
  # get unique protein groups
  unique_PG <- unique(data$PG.ProteinGroups)

  effect_size_out <- c()
  # condition comp. loop
  for(i in 1:ncol(condition_comparisons)){
    statusprogressbar(run = i, max.run = ncol(condition_comparisons), width = 60L, info = "calc. effect sizes...")
    # filter data for per condition comparison
    tmp_data <-data %>% dplyr::filter(.data$R.Condition%in%condition_comparisons[,i])
    # generate peptide mean scaled data
    tmp_data <- tmp_data %>%
      group_by(.data$PG.ProteinGroups,
               .data$PEP.StrippedSequence) %>%
      dplyr::mutate(mean_scaled_peptide_intensity= .data$peptide_intensity/mean(.data$peptide_intensity,na.rm = T)) %>%
      ungroup() %>%
      dplyr::mutate(R.Condition = factor(.data$R.Condition,
                                         levels = c(condition_comparisons[1,i],condition_comparisons[2,i])))

    # protein loop
    tmp_single_out_bind <- c()
    for(k in 1:length(unique_PG)){
      tmp_data_single <- tmp_data %>%
        dplyr::filter(.data$PG.ProteinGroups==unique_PG[k])
      # calc. single effect size
      tmp_cohend<- effsize::cohen.d(d = tmp_data_single$mean_scaled_peptide_intensity,
                                     f = tmp_data_single$R.Condition,
                                     pooled=TRUE,
                                     paired=FALSE,
                                     na.rm=FALSE,
                                     mu=0,
                                     hedges.correction=FALSE,
                                     conf.level=0.95,
                                     noncentral=FALSE,
                                     within=T)
      #bind single
      tmp_single_out_bind <- bind_rows(tmp_single_out_bind,
                                       tibble(PG.ProteinGroups = unique_PG[k],
                                              slr_ratio_meta = paste(condition_comparisons[1,i],"/",condition_comparisons[2,i],sep=""),
                                              effect_size_method = tmp_cohend$method,
                                              d = tmp_cohend$estimate,
                                              d_pooled_SD = tmp_cohend$sd,
                                              d_95CI_lower = tmp_cohend$conf.int[1],
                                              d_95CI_upper = tmp_cohend$conf.int[2],
                                              d_magnitute = tmp_cohend$magnitude)
      )
    } # end loop over proteins

    #bind comparisons
    effect_size_out <- bind_rows(effect_size_out,
                                 tmp_single_out_bind)

  }# end loop over condition comparisons
  return(effect_size_out)
} # end function
