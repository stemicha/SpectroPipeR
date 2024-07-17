#' SpectroPipeR: MVA module
#' @description
#' Function for conducting multivariate analysis, which represents the third step in the pipeline and builds upon the normalization and quantification module from step 2.
#'
#' @param SpectroPipeR_data_quant it is the SpectroPipeR_data_quant list object from norm_quant_module() object e.g. `data_input_norm_quant` see example below
#' @param HCPC_analysis boolean; should a HCPC be performed or not
#' @param costum_colors if you would like to use your own colors please provide a named color vector (e.g. c(condition1 = "black", condition2 = "grey")); names should have the same naming and length like the conditions set in Spectronaut
#' @return saves MVA analysis results in output folder
#' @import ggplot2
#' @import tidyverse
#' @import readr
#' @import tidyr
#' @importFrom methods is
#' @importFrom GGally ggpairs ggally_points ggally_barDiag
#' @import FactoMineR
#' @import factoextra
#' @import umap
#' @import ggrepel
#' @import ggcorrplot
#' @import tibble
#' @importFrom stats median sd lm as.formula quantile cor na.omit
#' @md
#'
#' @returns SpectroPipeR_data_MVA list object with the MVA analysis data in addition to the automatically saved tables and plots
#'  For the description of the generated figures and tables please read the manual & vignettes
#'
#' | <u> __list element__ </u> | <u> __description__ </u>              |
#' |:--------------------------|:--------------------------------------|
#' | PCA_peptide_intensity          | *PCA list element:* PCA list element of peptide int. |
#' | PCA_protein_intensity          | *PCA list element:* PCA list element of protein int. |
#' | UMAP_protein_intensity         | *umap element:* UMAP element of protein int. |
#' | peptide_intensity_correlation  | *matrix:* Spearman correlation scores of peptide int. |
#' | protein_intensity_correlation  | *matrix:* Spearman correlation scores of protein int. |
#'
#' @export
#'
#' @examples
#'#load library
#'library(SpectroPipeR)
#'
#'# use default parameters list
#'params <- list(output_folder = "../SpectroPipeR_test_folder")
#'
#'# example input file
#'example_file_path <- system.file("extdata",
#'                                 "SN_test_HYE_mix_file.tsv",
#'                                 package="SpectroPipeR")
#'
#'# step 1: load Spectronaut data module
#'SpectroPipeR_data <- read_spectronaut_module(file = example_file_path,
#'                                             parameter = params,
#'                                             print.plot = FALSE)
#'
#'# step 2: normalize & quantification module
#'SpectroPipeR_data_quant <- norm_quant_module(SpectroPipeR_data = SpectroPipeR_data)
#'
#'# step 3: MVA module
#'SpectroPipeR_MVA <- MVA_module(SpectroPipeR_data_quant = SpectroPipeR_data_quant,
#'           HCPC_analysis = FALSE)

MVA_module <- function(SpectroPipeR_data_quant = NULL,
                       HCPC_analysis = FALSE,
                       costum_colors = NULL){


  # check if SpectroPipeR_norm_quant object was provided
  if(!is(SpectroPipeR_data_quant,"SpectroPipeR_norm_quant")){
    stop_text <- paste("no valid SpectroPipeR_norm_quant object was provided: please process data with norm_quant_module() to generate one!")
    message_function(text = stop_text,color = "red",log_file_name = log_file_name)
    stop()
  }

  #output folder
  out_folder <- SpectroPipeR_data_quant$parameter$output_folder

  #log file name
  log_file_name <- SpectroPipeR_data_quant$parameter$log_file_name

  message_function(text = "",
                   color = "white",
                   log_file_name = log_file_name)
  message_function(text = "#*****************************************",
                   color = "white",
                   log_file_name = log_file_name)
  message_function(text = "# MVA MODULE",
                   color = "white",
                   log_file_name = log_file_name)
  message_function(text = "#*****************************************",
                   color = "white",
                   log_file_name = log_file_name)
  message_function(text = "",
                   color = "white",
                   log_file_name = log_file_name)

  message_function(text = "perfoming multivariate analysis ...",
                   color = "blue",
                   log_file_name = log_file_name)

  #number of samples detected  ====
  sample_length <- length(unique(SpectroPipeR_data_quant$data_input_normalized$R.FileName))
  condition_length <- length(unique(SpectroPipeR_data_quant$data_input_normalized$R.Condition))

  #colors setup ====
  if(is.null(costum_colors)){
    general.colors <- c("#1B5E20","#C62828","#0D47A1","#4E342E","#E65100","#000000","#6A1B9A","#4F8F00","#89694E","#00695C","#827717","#5E5E5E")
    general.colors.gradient <- colorRampPalette(c("black","#0057e7","#d62d20","#008744","#ffa700","grey"))
    condition_elements_length <- length(unique(SpectroPipeR_data_quant$data_input_normalized$R.Condition))
    if(condition_elements_length > length(general.colors)){
      condition_colors <- general.colors.gradient(condition_elements_length)
    }else{
      condition_colors <- general.colors[1:condition_elements_length]
    }
    names(condition_colors) <- unique(SpectroPipeR_data_quant$data_input_normalized$R.Condition)

  }else{
    if(is.vector(costum_colors)==T){
      number_of_conditions <- unique(SpectroPipeR_data_quant$data_input_normalized$R.Condition)
      if(length(costum_colors)==length(number_of_conditions)){
        message_function(text = "use custom colors ...",
                         color = "green",
                         log_file_name = log_file_name)
        condition_colors <- costum_colors
      }
    }else{
      stop_text <- paste("custom colors is not a valid vector")
      message_function(text = stop_text,color = "red",log_file_name = log_file_name)
      stop()
    }
  }

  #create measurement order color vector ====
  # Spectral colors
  measurement_colors <- c("#9E0142","#D53E4F","#F46D43","#FDAE61","#FEE08B","#E6F598","#ABDDA4","#66C2A5","#3288BD","#5E4FA2")

  #create file number specific folder ====
  if(dir.exists(paste0(out_folder,"/","04_multivariate_analysis/",sample_length,"_sample_analysis"))){
    message_function(text = paste0(out_folder,"/","04_multivariate_analysis/",sample_length,"_sample_analysis"," ATTENTION !!! ---> folder already exists - files will be replaced !!!"),
                     color = "magenta",
                     log_file_name = log_file_name)
  }else{
    dir.create(paste0(out_folder,"/","04_multivariate_analysis/"),recursive = T,showWarnings = F)
    dir.create(paste0(out_folder,"/","04_multivariate_analysis/",sample_length,"_sample_analysis"),recursive = T,showWarnings = F)

  }



  #habilliage ====
  # add measurement order
  message_function(text = "generating raw file / condition file ...",
                   color = "blue",
                   log_file_name = log_file_name)

  habil <- SpectroPipeR_data_quant$data_input_normalized %>%
            dplyr::distinct(.data$R.Condition,.data$R.FileName,.data$R.Replicate,.data$`R.Run Date`) %>%
            dplyr::arrange(.data$R.Condition) %>%
            dplyr::mutate(`R.Run Date` = lubridate::as_datetime(.data$`R.Run Date`)) %>%
            dplyr::mutate(measurement_order = dense_rank(.data$`R.Run Date`))


  # save habil data
  readr::write_csv(x = habil,file = paste0(out_folder,"/","04_multivariate_analysis/",sample_length,"_sample_analysis/sample_to_condition_file.csv"))


  message_function(text = "PCA analysis on peptide level ...",
                   color = "blue",
                   log_file_name = log_file_name)

  #PCA analysis on peptide level ====
  peptide_data_PCA   <-   SpectroPipeR_data_quant$peptide_intensity %>%
                          dplyr::select(.data$R.FileName,
                                        .data$EG.ModifiedPeptide,
                                        .data$peptide_intensity) %>%
    tidyr::pivot_wider(names_from = .data$EG.ModifiedPeptide,
                       values_from  = .data$peptide_intensity)
  peptide_data_PCA <- base::as.data.frame(peptide_data_PCA)
  base::rownames(peptide_data_PCA)<-peptide_data_PCA$R.FileName #set rownames
  peptide_data_PCA <- peptide_data_PCA[,-1] #remove file name column
  peptide_data_PCA_res <-  FactoMineR::PCA(log2(peptide_data_PCA), graph=FALSE,scale.unit=T,ncp = 10)
  readr::write_rds(x = peptide_data_PCA_res,
                   file = paste(out_folder,"/","04_multivariate_analysis/",sample_length,"_sample_analysis/PCA_analysis__peptide_level.RDS",sep=""))

  message_function(text = "PCA analysis on protein level ...",
                   color = "blue",
                   log_file_name = log_file_name)

  #PCA analysis on protein level ====
  protein_data_PCA   <-   SpectroPipeR_data_quant$protein_data %>%
    dplyr::select(.data$R.FileName,
                  .data$PG.ProteinGroups,
                  .data$protein_intensity) %>%
    tidyr::pivot_wider(names_from = .data$R.FileName,
                values_from = .data$protein_intensity)
  protein_data_PCA <- base::as.data.frame(protein_data_PCA)
  base::rownames(protein_data_PCA)<-protein_data_PCA$PG.ProteinGroups #set rownames
  protein_data_PCA <- protein_data_PCA[,-1] #remove file name column
  proteins_counts_pre <- dim(protein_data_PCA)[1]
  protein_data_PCA <- stats::na.omit(protein_data_PCA) # remove NA values / maxLFQ can produce NAs
  proteins_counts_post <- dim(protein_data_PCA)[1]
  protein_data_PCA <- t(protein_data_PCA) #transpose

  # proteins removed proteins removed due to NAs logging
  message_function(text = paste("...",proteins_counts_pre-proteins_counts_post,"proteins removed due to NAs"),
                   color = "blue",
                   log_file_name = log_file_name)

  protein_data_PCA_res <-  FactoMineR::PCA(log2(protein_data_PCA), graph=FALSE,scale.unit=T,ncp = 10)
  readr::write_rds(x = protein_data_PCA_res,file = paste(out_folder,"/","04_multivariate_analysis/",sample_length,"_sample_analysis/PCA_analysis__protein_level.RDS",sep=""))

  # resolving only 1 replicate per condition issue
  # >> PCA has only 1 dimension on peptide level
  if(ncol(protein_data_PCA_res$ind$coord)==1){
    protein_data_PCA_res_only_1dim <- TRUE
  }else{
    protein_data_PCA_res_only_1dim <- FALSE
  }

  if(ncol(peptide_data_PCA_res$ind$coord)==1){
    peptide_data_PCA_res_only_1dim <- TRUE
  }else{
    peptide_data_PCA_res_only_1dim <- FALSE
  }

  #Screeplot ====
  message_function(text = "generating Scree plots ...",
                   color = "blue",
                   log_file_name = log_file_name)

  protein_data_PCA_res_screePlot <- factoextra::fviz_eig(protein_data_PCA_res,
                                                         choice = c("variance"),
                                                         addlabels=TRUE,
                                                         hjust = -0.3,
                                                         ncp = 10) +
    labs(title="Scree plot (protein level)", subtitle="explained variances")+
    ylim(c(0,100))+
    theme_bw(base_size = 18)

  peptide_data_PCA_res_screePlot <- factoextra::fviz_eig(peptide_data_PCA_res,
                                                       choice = c("variance"),
                                                       addlabels=TRUE,
                                                       hjust = -0.3,
                                                       ncp = 10) +
    labs(title="Scree plot (peptide level)", subtitle="explained variances")+
    ylim(c(0,100))+
    theme_bw(base_size = 18)


  #Biplot top 10 elements ====
  #contribution plot
  message_function(text = "generating Top10 variable contribution plot ...",
                   color = "blue",
                   log_file_name = log_file_name)

  if(peptide_data_PCA_res_only_1dim == T){
    peptide_data_PCA_res_biplot <- ggplot()+
      labs(title=paste("TOP", 10," >>> NOT POSSIBLE with provided data!"),
           subtitle="contributions to PCA (peptide level)")+
      theme_bw(base_size = 18)
  }else{
    peptide_data_PCA_res_biplot <- factoextra::fviz_pca_var(peptide_data_PCA_res,
                                                            select.var = list(contrib = 10),
                                                            col.var = "contrib",
                                                            col.circle = "black")+
      labs(title=paste("TOP", 10), subtitle="contributions to PCA (peptide level)")+
      theme_bw(base_size = 18)
  }



  # resolving only 1 replicate per condition issue
  # >> PCA has only 1 dimension on peptide level
  if(protein_data_PCA_res_only_1dim == T){
    protein_data_PCA_res_biplot <- ggplot()+
      labs(title=paste("TOP", 10," >>> NOT POSSIBLE with provided data!"),
           subtitle="contributions to PCA (protein level)")+
      theme_bw(base_size = 18)

  }else{
  protein_data_PCA_res_biplot <- factoextra::fviz_pca_var(protein_data_PCA_res,
                                                          select.var = list(contrib = 10),
                                                          col.var = "contrib",
                                                          col.circle = "black")+
    labs(title=paste("TOP", 10), subtitle="contributions to PCA (protein level)")+
    theme_bw(base_size = 18)

  }

  #3D PCA plot  =====
  message_function(text = "generating 3D PCA plot ...",
                   color = "blue",
                   log_file_name = log_file_name)
    #protein level
  plotly3d_dat_protein <- as.data.frame(protein_data_PCA_res$ind$coord)

  #PCA condition lookup ====

  habil<- habil %>%
    dplyr::arrange(
      match(.data$R.FileName,
            rownames(plotly3d_dat_protein)
            )
      )

  ## 3D PCA plot individuals

  # resolving only 1 replicate per condition issue
  # >> PCA has only 1 dimension on peptide level
  if(dim(protein_data_PCA_res$ind$coord)[2]==1){

    #add dummy data
    plotly3d_dat_protein$Dim.2 <- 0
    plotly3d_dat_protein$Dim.3 <- 0
    plotly3d_ind_plot_protein <- plotly::plot_ly(plotly3d_dat_protein,
                                                 x = ~Dim.1,
                                                 y = ~Dim.2,
                                                 z = ~Dim.3,
                                                 color = ~as.factor(habil$R.Condition),
                                                 colors = condition_colors) %>%
      plotly::add_markers() %>%
      plotly::add_text(text = rownames(plotly3d_dat_protein), textposition = "top right",visible = "legendonly") %>%
      plotly::layout(title = "PCA (protein level)",
                     scene = list(xaxis = list(title = paste("1st PC (", round(protein_data_PCA_res$eig[1,2]), "%)", sep="")),
                                  yaxis = list(title = paste("2nd PC (", NA, "%)", sep="")),
                                  zaxis = list(title = paste("3rd PC (", NA, "%)", sep=""))))

  }else{
    plotly3d_ind_plot_protein <- plotly::plot_ly(as.data.frame(protein_data_PCA_res$ind$coord),
                                                 x = ~Dim.1,
                                                 y = ~Dim.2,
                                                 z = ~Dim.3,
                                                 color = ~as.factor(habil$R.Condition),
                                                 colors = condition_colors) %>%
      plotly::add_markers() %>%
      plotly::add_text(text = rownames(plotly3d_dat_protein), textposition = "top right",visible = "legendonly") %>%
      plotly::layout(title = "PCA (protein level)",
                     scene = list(xaxis = list(title = paste("1st PC (", round(protein_data_PCA_res$eig[1,2]), "%)", sep="")),
                                  yaxis = list(title = paste("2nd PC (", round(protein_data_PCA_res$eig[2,2]), "%)", sep="")),
                                  zaxis = list(title = paste("3rd PC (", round(protein_data_PCA_res$eig[3,2]), "%)", sep=""))))

  }

  #peptide level
  plotly3d_dat_peptide <- as.data.frame(peptide_data_PCA_res$ind$coord)

  # resolving only 1 replicate per condition issue
  # >> PCA has only 1 dimension on peptide level
  if(dim(protein_data_PCA_res$ind$coord)[2]==1){

    #add dummy data
    plotly3d_dat_peptide$Dim.2 <- 0
    plotly3d_dat_peptide$Dim.3 <- 0
    plotly3d_ind_plot_peptide <- plotly::plot_ly(plotly3d_dat_peptide,
                                                 x = ~Dim.1,
                                                 y = ~Dim.2,
                                                 z = ~Dim.3,
                                                 color = ~as.factor(habil$R.Condition),
                                                 colors = condition_colors) %>%
      plotly::add_markers() %>%
      plotly::add_text(text = rownames(plotly3d_dat_protein), textposition = "top right",visible = "legendonly") %>%
      plotly::layout(title = "PCA (peptide level)",
                     scene = list(xaxis = list(title = paste("1st PC (", round(peptide_data_PCA_res$eig[1,2]), "%)", sep="")),
                                  yaxis = list(title = paste("2nd PC (", NA, "%)", sep="")),
                                  zaxis = list(title = paste("3rd PC (", NA, "%)", sep=""))))

  }else{
    plotly3d_ind_plot_peptide <- plotly::plot_ly(as.data.frame(peptide_data_PCA_res$ind$coord),
                                                 x = ~Dim.1,
                                                 y = ~Dim.2,
                                                 z = ~Dim.3,
                                                 color = ~as.factor(habil$R.Condition),
                                                 colors = condition_colors) %>%
      plotly::add_markers() %>%
      plotly::add_text(text = rownames(plotly3d_dat_protein), textposition = "top right",visible = "legendonly") %>%
      plotly::layout(title = "PCA (peptide level)",
                     scene = list(xaxis = list(title = paste("1st PC (", round(peptide_data_PCA_res$eig[1,2]), "%)", sep="")),
                                  yaxis = list(title = paste("2nd PC (", round(peptide_data_PCA_res$eig[2,2]), "%)", sep="")),
                                  zaxis = list(title = paste("3rd PC (", round(peptide_data_PCA_res$eig[3,2]), "%)", sep=""))))

  }

  #PCA plot 1st vs. 2nd dimension ====
  message_function(text = "generating PCA plot 1st vs. 2nd dimension ...",
                   color = "blue",
                   log_file_name = log_file_name)

  if(peptide_data_PCA_res_only_1dim==TRUE){
    pca_plot_peptide <- ggplot()+
      theme_classic(base_size = 18)+
      labs(title="PCA (peptide level) >>> could not be calc. \ndue to number of replicates",
           subtitle="individuals coordinates",color="condition")+
      guides(fill="none")
  }else{
    pca_plot_peptide<- factoextra::fviz_pca_ind(peptide_data_PCA_res,axes = c(1, 2),
                                                mean.point=F,pointsize=3,
                                                pointshape=19,
                                                label="none",
                                                ellipse.type = "convex",
                                                ellipse.alpha=0.2,
                                                addEllipses = T,
                                                ellipse.level=0.95,
                                                habillage = as.factor(habil$R.Condition))+
      scale_color_manual(values = condition_colors)+
      scale_fill_manual(values = condition_colors)+
      theme_classic(base_size = 18)+
      labs(title="PCA (peptide level)",subtitle="individuals coordinates",color="condition")+
      guides(fill="none")
  }

  if(protein_data_PCA_res_only_1dim==TRUE){
    pca_plot_protein <- ggplot()+
      theme_classic(base_size = 18)+
      labs(title="PCA (protein level) >>> could not be calc. \ndue to number of replicates",
           subtitle="individuals coordinates",color="condition")+
      guides(fill="none")
  }else{
    pca_plot_protein<- factoextra::fviz_pca_ind(protein_data_PCA_res,
                                                axes = c(1, 2),
                                                mean.point=F,
                                                pointshape=19,
                                                pointsize=3,
                                                label="none",
                                                ellipse.type = "convex",
                                                ellipse.alpha=0.2,
                                                addEllipses = T,
                                                ellipse.level=0.95,
                                                habillage = as.factor(habil$R.Condition))+
      scale_color_manual(values = condition_colors)+
      scale_fill_manual(values = condition_colors)+
      theme_classic(base_size = 18)+
      labs(title="PCA (protein level)",subtitle="individuals coordinates",color="condition")+
      guides(fill="none")
  }



  # protein PCA plot 1vs2 dimension with replicates labels

  # get individuals data
  protein_data_PCA_res_ind_data <- as.data.frame(protein_data_PCA_res$ind$coord)
  protein_data_PCA_res_ind_data$R.FileName <- rownames(protein_data_PCA_res_ind_data)
  protein_data_PCA_res_ind_data <- tibble::as_tibble(protein_data_PCA_res_ind_data)

  # add condition file info
  protein_data_PCA_res_ind_data <- dplyr::left_join(protein_data_PCA_res_ind_data,habil,by = "R.FileName")

  # generate the convex hull
  if(protein_data_PCA_res_only_1dim==TRUE){
    protein_data_PCA_res_hull_group <- NA
  }else{
    protein_data_PCA_res_hull_group <- protein_data_PCA_res_ind_data %>%
      group_by(.data$R.Condition) %>%
      slice(chull(.data$Dim.1,
                  .data$Dim.2))
  }


  # slice --> dplyr function
  # chull --> grDevice function
  ## Computes the subset of points which lie on the convex hull of the set of points specified.


  # Define the scatterplot
  if(protein_data_PCA_res_only_1dim==TRUE){
    pca_plot_protein_replicate <- ggplot(data = protein_data_PCA_res_ind_data,
                                         mapping = aes(x = .data$Dim.1,
                                                       y = rep(0,nrow(protein_data_PCA_res_ind_data)),
                                                       color = .data$R.Condition,
                                                       fill = .data$R.Condition)) +
      geom_hline(yintercept = 0,linetype = "dashed")+
      geom_vline(xintercept = 0,linetype = "dashed")+
      theme_classic(base_size = 16)+
      geom_point(size=3)+ # add point shape R.Replicate
      scale_color_manual(values = condition_colors)+
      scale_fill_manual(values = condition_colors)+
      labs(title = "PCA (protein level)",
           subtitle = "protein intensities",
           caption = "replicates highlighted as text labels",
           x = paste("Dim.1 (",round(protein_data_PCA_res$eig[1,2],digits = 2)," %)",sep=""),# labeling of axis
           y = paste("Dim.2 (",NA," %)",sep=""),
           # get expl. variance  Dim.xxx data by pca$eig[xxx,2]
           color = "Conditions",
           fill = "Conditions")+
      ggrepel::geom_text_repel(mapping = aes(label = .data$R.Replicate),
                               size = 4,
                               min.segment.length = 0.01,
                               box.padding = 1)+
      guides(color = "none")

  }else{
    pca_plot_protein_replicate <- ggplot(data = protein_data_PCA_res_ind_data,
                                         mapping = aes(x = .data$Dim.1,
                                                       y = .data$Dim.2,
                                                       color = .data$R.Condition,
                                                       fill = .data$R.Condition)) +
      geom_hline(yintercept = 0,linetype = "dashed")+
      geom_vline(xintercept = 0,linetype = "dashed")+
      theme_classic(base_size = 16)+
      geom_point(size=3)+ # add point shape R.Replicate
      geom_polygon(data = protein_data_PCA_res_hull_group,
                   alpha = 0.5)+ # add convex shape
      scale_color_manual(values = condition_colors)+
      scale_fill_manual(values = condition_colors)+
      labs(title = "PCA (protein level)",
           subtitle = "protein intensities",
           caption = "replicates highlighted as text labels",
           x = paste("Dim.1 (",round(protein_data_PCA_res$eig[1,2],digits = 2)," %)",sep=""),# labeling of axis
           y = paste("Dim.2 (",round(protein_data_PCA_res$eig[2,2],digits = 2)," %)",sep=""),
           # get expl. variance  Dim.xxx data by pca$eig[xxx,2]
           color = "Conditions",
           fill = "Conditions")+
      ggrepel::geom_text_repel(mapping = aes(label = .data$R.Replicate),
                               size = 4,
                               min.segment.length = 0.01,
                               box.padding = 1)+
      guides(color = "none")
  }


  ggsave_pdf_png(filename = paste0(out_folder,"/","04_multivariate_analysis/",sample_length,"_sample_analysis/PCA_plot_protein_level_replicates_marked"),
         plot = pca_plot_protein_replicate,
         width = 12,
         height = 9)

  # protein PCA plot 1vs2 dimension with measurement order
  if(protein_data_PCA_res_only_1dim==TRUE){
    pca_plot_protein_measurement_order <- ggplot(data = protein_data_PCA_res_ind_data,
                                                 mapping = aes(x = .data$Dim.1,
                                                               y = rep(0,nrow(protein_data_PCA_res_ind_data)),
                                                               color = .data$measurement_order)) +
      geom_hline(yintercept = 0,linetype = "dashed")+
      geom_vline(xintercept = 0,linetype = "dashed")+
      theme_classic(base_size = 16)+
      geom_point(size=3)+ # add point shape R.Replicate
      scale_color_gradientn(colours = measurement_colors)+
      labs(title = "PCA (protein level)",
           subtitle = "protein intensities",
           x = paste("Dim.1 (",round(protein_data_PCA_res$eig[1,2],digits = 2)," %)",sep=""),# labeling of axis
           y = paste("Dim.2 (",NA," %)",sep=""),
           # get expl. variance  Dim.xxx data by pca$eig[xxx,2]
           color = "measurement\norder")
  }else{
    pca_plot_protein_measurement_order <- ggplot(data = protein_data_PCA_res_ind_data,
                                                 mapping = aes(x = .data$Dim.1,
                                                               y = .data$Dim.2,
                                                               color = .data$measurement_order)) +
      geom_hline(yintercept = 0,linetype = "dashed")+
      geom_vline(xintercept = 0,linetype = "dashed")+
      theme_classic(base_size = 16)+
      geom_point(size=3)+ # add point shape R.Replicate
      scale_color_gradientn(colours = measurement_colors)+
      labs(title = "PCA (protein level)",
           subtitle = "protein intensities",
           x = paste("Dim.1 (",round(protein_data_PCA_res$eig[1,2],digits = 2)," %)",sep=""),# labeling of axis
           y = paste("Dim.2 (",round(protein_data_PCA_res$eig[2,2],digits = 2)," %)",sep=""),
           # get expl. variance  Dim.xxx data by pca$eig[xxx,2]
           color = "measurement\norder")
  }

  ggsave_pdf_png(filename = paste0(out_folder,"/","04_multivariate_analysis/",sample_length,"_sample_analysis/PCA_plot_protein_level_measurement_order"),
                 plot = pca_plot_protein_measurement_order,
                 width = 12,
                 height = 9)

  # protein PCA plot 1vs2 dimension with condition labels
  if(protein_data_PCA_res_only_1dim==TRUE){
    pca_plot_protein_condition <- ggplot(data = protein_data_PCA_res_ind_data,
                                         mapping = aes(x = .data$Dim.1,
                                                       y = rep(0,nrow(protein_data_PCA_res_ind_data)),
                                                       color = .data$R.Condition,
                                                       fill = .data$R.Condition)) +
      geom_hline(yintercept = 0,linetype = "dashed")+
      geom_vline(xintercept = 0,linetype = "dashed")+
      theme_classic(base_size = 16)+
      geom_point(size=3)+ # add point shape R.Replicate
      scale_color_manual(values = condition_colors)+
      scale_fill_manual(values = condition_colors)+
      labs(title = "PCA (protein level)",
           subtitle = "protein intensities",
           caption = "condition highlighted as text labels",
           x = paste("Dim.1 (",round(protein_data_PCA_res$eig[1,2],digits = 2)," %)",sep=""),# labeling of axis
           y = paste("Dim.2 (",NA," %)",sep=""),
           # get expl. variance  Dim.xxx data by pca$eig[xxx,2]
           color = "Conditions",
           fill = "Conditions")+
      ggrepel::geom_label_repel(data = protein_data_PCA_res_ind_data %>%
                                  group_by(.data$R.Condition) %>%
                                  dplyr::summarise(mean_Dim.1 = mean(.data$Dim.1),
                                                   mean_Dim.2 = 0) %>%
                                  ungroup(),
                                mapping = aes(x = .data$mean_Dim.1,
                                              y = 0,
                                              label = .data$R.Condition),
                                color="white",
                                segment.colour = "black",
                                size = 4,
                                alpha = 0.8,
                                min.segment.length = 0.01,
                                box.padding = 1,show.legend = FALSE)+
      guides(color = "none", text = "none")
  }else{
    pca_plot_protein_condition <- ggplot(data = protein_data_PCA_res_ind_data,
                                         mapping = aes(x = .data$Dim.1,
                                                       y = .data$Dim.2,
                                                       color = .data$R.Condition,
                                                       fill = .data$R.Condition)) +
      geom_hline(yintercept = 0,linetype = "dashed")+
      geom_vline(xintercept = 0,linetype = "dashed")+
      theme_classic(base_size = 16)+
      geom_point(size=3)+ # add point shape R.Replicate
      geom_polygon(data = protein_data_PCA_res_hull_group,
                   alpha = 0.5)+ # add convex shape
      scale_color_manual(values = condition_colors)+
      scale_fill_manual(values = condition_colors)+
      labs(title = "PCA (protein level)",
           subtitle = "protein intensities",
           caption = "condition highlighted as text labels",
           x = paste("Dim.1 (",round(protein_data_PCA_res$eig[1,2],digits = 2)," %)",sep=""),# labeling of axis
           y = paste("Dim.2 (",round(protein_data_PCA_res$eig[2,2],digits = 2)," %)",sep=""),
           # get expl. variance  Dim.xxx data by pca$eig[xxx,2]
           color = "Conditions",
           fill = "Conditions")+
      ggrepel::geom_label_repel(data = protein_data_PCA_res_ind_data %>%
                                  group_by(.data$R.Condition) %>%
                                  dplyr::summarise(mean_Dim.1 = mean(.data$Dim.1),
                                                   mean_Dim.2 = mean(.data$Dim.2)) %>%
                                  ungroup(),
                                mapping = aes(x = .data$mean_Dim.1,
                                              y = .data$mean_Dim.2,
                                              label = .data$R.Condition),
                                color="white",
                                segment.colour = "black",
                                size = 4,
                                alpha = 0.8,
                                min.segment.length = 0.01,
                                box.padding = 1,show.legend = FALSE)+
      guides(color = "none", text = "none")
  }


  ggsave_pdf_png(filename = paste0(out_folder,"/","04_multivariate_analysis/",sample_length,"_sample_analysis/PCA_plot_protein_level_conditions_marked"),
                 plot = pca_plot_protein_condition,
                 width = 12,
                 height = 9)



  #PCA plot 1st to 5th dimension ====
  message_function(text = "generating PCA plot 1st to 5th dimension ...",
                   color = "blue",
                   log_file_name = log_file_name)

  #prepare data for ggpairs >>> peptide
  if(peptide_data_PCA_res_only_1dim==TRUE){
    peptide_PCA_1_to_5dim <- as.data.frame(peptide_data_PCA_res$ind$coord)
    peptide_PCA_1_to_5dim <- dplyr::bind_cols(peptide_PCA_1_to_5dim,
                                          tibble::tibble(Dim.2 = rep(0,nrow(peptide_PCA_1_to_5dim)),
                                              Dim.3 = rep(0,nrow(peptide_PCA_1_to_5dim)),
                                              Dim.4 = rep(0,nrow(peptide_PCA_1_to_5dim)),
                                              Dim.5 = rep(0,nrow(peptide_PCA_1_to_5dim)),
                                              )
                                       )
    peptide_PCA_1_to_5dim$R.FileName <- rownames(peptide_PCA_1_to_5dim)
    peptide_PCA_1_to_5dim <- left_join(tibble::as_tibble(peptide_PCA_1_to_5dim),
                                       habil,
                                       by="R.FileName") %>%
      dplyr::select(-.data$R.FileName)
    ncol_peptide <- 5
  }else{
    ncol_peptide <- ifelse(test = ncol(peptide_data_PCA_res$ind$coord)>5,yes = 5,no = ncol(peptide_data_PCA_res$ind$coord))
    peptide_PCA_1_to_5dim <- as.data.frame(peptide_data_PCA_res$ind$coord)[,1:ifelse(ncol_peptide<=5,ncol_peptide,5)]
    peptide_PCA_1_to_5dim$R.FileName <- rownames(peptide_PCA_1_to_5dim)
    peptide_PCA_1_to_5dim <- left_join(tibble::as_tibble(peptide_PCA_1_to_5dim),
                                       habil,
                                       by="R.FileName") %>%
      dplyr::select(-.data$R.FileName)
  }


  #peptide ggpairs plot 1-5 dim.
  ggplotPCA_peptide_PCA_1_to_5dim <- GGally::ggpairs(peptide_PCA_1_to_5dim,
                                         columns = 1:ncol_peptide,
                                         legend = c(1,1), #legend form second row and first column
                                         title = paste0("PCA 1st to ",ncol_peptide,"th dim. (peptide level)"),
                                         mapping = aes(color = .data$R.Condition,
                                                       fill = .data$R.Condition),
                                         upper  = list(continuous = "blank"),
                                         lower = list(continuous = function(data, mapping, ...) {
                                           GGally::ggally_points(data = data, mapping = mapping, alpha = .8,size=2) +
                                             scale_colour_manual(values = condition_colors)+
                                             geom_vline(xintercept = 0,linetype="dotted")+
                                             geom_hline(yintercept = 0,linetype="dotted")
                                                }
                                              ),
                                         diag = list(continuous = function(data, mapping, ...) {
                                           ggplot(data = data, mapping = mapping)+geom_density(alpha = .4) +
                                             scale_colour_manual(values = condition_colors) +
                                             scale_fill_manual(values = condition_colors)
                                                }
                                              )
                                          )+
                                        theme_bw(base_size = 18)+
                                        theme(legend.position = "right")

  #peptide ggpairs plot 1-5 dim. // measurement order
  ggpairs_colors_measurement_order <- colorRampPalette(colors = measurement_colors)(max(habil$measurement_order))
  names(ggpairs_colors_measurement_order) <- seq(1,max(habil$measurement_order))

  ggplotPCA_peptide_PCA_1_to_5dim_measurement_order <- GGally::ggpairs(peptide_PCA_1_to_5dim,
                                             columns = 1:ncol_peptide,
                                             #legend = c(1,1), #legend form second row and first column
                                             title = paste0("PCA 1st to ",ncol_peptide,"th dim. (peptide level)"),
                                             mapping = aes(color=as.factor(.data$measurement_order),
                                                           fill=as.factor(.data$measurement_order)),
                                             upper  = list(continuous = "blank"),
                                             lower = list(continuous = function(data, mapping, ...) {
                                               GGally::ggally_points(data = data, mapping = mapping, alpha = .8,size=2) +
                                                 scale_colour_manual(values = ggpairs_colors_measurement_order)+
                                                 geom_vline(xintercept = 0,linetype="dotted")+
                                                 geom_hline(yintercept = 0,linetype="dotted")
                                             }
                                             ),
                                             diag = list(continuous = function(data, mapping, ...) {
                                               GGally::ggally_barDiag(data = data, mapping = mapping,bins = 15,rescale = F,color=NA) +
                                                 scale_fill_manual(values = ggpairs_colors_measurement_order)

                                             })
  )+
    theme_bw(base_size = 18)+
    theme(legend.position = "right")


  #prepare data for ggpairs >>> protein
  if(protein_data_PCA_res_only_1dim==TRUE){
    protein_PCA_1_to_5dim <- as.data.frame(protein_data_PCA_res$ind$coord)
    protein_PCA_1_to_5dim <- dplyr::bind_cols(protein_PCA_1_to_5dim,
                                       tibble::tibble(Dim.2 = rep(0,nrow(protein_PCA_1_to_5dim)),
                                              Dim.3 = rep(0,nrow(protein_PCA_1_to_5dim)),
                                              Dim.4 = rep(0,nrow(protein_PCA_1_to_5dim)),
                                              Dim.5 = rep(0,nrow(protein_PCA_1_to_5dim)),
                                       )
    )
    protein_PCA_1_to_5dim$R.FileName <- rownames(protein_PCA_1_to_5dim)
    protein_PCA_1_to_5dim <- left_join(tibble::as_tibble(protein_PCA_1_to_5dim),
                                       habil,
                                       by="R.FileName") %>%
      dplyr::select(-.data$R.FileName)
    ncol_protein <- 5
  }else{
    ncol_protein <- ifelse(test = ncol(protein_data_PCA_res$ind$coord)>5,yes = 5,no = ncol(protein_data_PCA_res$ind$coord))
    protein_PCA_1_to_5dim <- as.data.frame(protein_data_PCA_res$ind$coord)[,1:ifelse(ncol_protein<=5,ncol_protein,5)]
    protein_PCA_1_to_5dim$R.FileName <- rownames(protein_PCA_1_to_5dim)
    protein_PCA_1_to_5dim <- left_join(tibble::as_tibble(protein_PCA_1_to_5dim),
                                       habil,
                                       by="R.FileName") %>%
      dplyr::select(-.data$R.FileName)
  }


  #peptide ggpairs plot 1-5 dim.
  ggplotPCA_protein_PCA_1_to_5dim <- GGally::ggpairs(protein_PCA_1_to_5dim,
                                             columns = 1:ncol_protein,
                                             title = paste0("PCA 1st to ",ncol_protein,"th dim. (protein level)"),
                                             legend = c(1,1), #legend form sencond row and first column
                                             mapping = aes(color = .data$R.Condition,
                                                           fill = .data$R.Condition),
                                             upper  = list(continuous = "blank"),
                                             lower = list(continuous = function(data, mapping, ...) {
                                               GGally::ggally_points(data = data, mapping = mapping, alpha = .8,size=2) +
                                                 scale_colour_manual(values = condition_colors)+
                                                 geom_vline(xintercept = 0,linetype="dotted")+
                                                 geom_hline(yintercept = 0,linetype="dotted")
                                             }
                                             ),
                                             diag = list(continuous = function(data, mapping, ...) {
                                               ggplot(data = data, mapping = mapping)+
                                                 geom_density(alpha = .4) +
                                                 scale_colour_manual(values = condition_colors) +
                                                 scale_fill_manual(values = condition_colors)
                                             }
                                             )
  )+
    theme_bw(base_size = 18)+
    theme(legend.position = "right")

  ggplotPCA_protein_PCA_1_to_5dim_measurement_order <- GGally::ggpairs(protein_PCA_1_to_5dim,
                                                               columns = 1:ncol_protein,
                                                               #legend = c(1,1), #legend form second row and first column
                                                               title = paste0("PCA 1st to ",ncol_protein,"th dim. (protein level)"),
                                                               mapping = aes(color=as.factor(.data$measurement_order),
                                                                             fill=as.factor(.data$measurement_order)),
                                                               upper  = list(continuous = "blank"),
                                                               lower = list(continuous = function(data, mapping, ...) {
                                                                 GGally::ggally_points(data = data, mapping = mapping, alpha = .8,size=2) +
                                                                   scale_colour_manual(values = ggpairs_colors_measurement_order)+
                                                                   geom_vline(xintercept = 0,linetype="dotted")+
                                                                   geom_hline(yintercept = 0,linetype="dotted")
                                                               }
                                                               ),
                                                               diag = list(continuous = function(data, mapping, ...) {
                                                                 GGally::ggally_barDiag(data = data, mapping = mapping,bins = 15,rescale = F,color=NA) +
                                                                   scale_fill_manual(values = ggpairs_colors_measurement_order)

                                                               })
  )+
    theme_bw(base_size = 18)+
    theme(legend.position = "right")


  #saving plot
  message_function(text = "saving PCA plot 1st to 5th dimension ...",
                   color = "blue",
                   log_file_name = log_file_name)

  ggsave_pdf_png(plot = ggplotPCA_peptide_PCA_1_to_5dim,
                 filename = paste0(out_folder,"/","04_multivariate_analysis/",sample_length,"_sample_analysis/PCA_plot_1st_to_5th_dimension__peptide_level"),
                 width = 20,
                 height = 15)
  ggsave_pdf_png(plot = ggplotPCA_protein_PCA_1_to_5dim,
                 filename = paste0(out_folder,"/","04_multivariate_analysis/",sample_length,"_sample_analysis/PCA_plot_1st_to_5th_dimension__protein_level"),
                 width = 20,
                 height = 15)

  ggsave_pdf_png(plot = ggplotPCA_peptide_PCA_1_to_5dim_measurement_order,
                 filename = paste0(out_folder,"/","04_multivariate_analysis/",sample_length,"_sample_analysis/PCA_plot_1st_to_5th_dimension__peptide_level_measurement_order"),
                 width = 20,
                 height = 15)
  ggsave_pdf_png(plot = ggplotPCA_protein_PCA_1_to_5dim_measurement_order,
                 filename = paste0(out_folder,"/","04_multivariate_analysis/",sample_length,"_sample_analysis/PCA_plot_1st_to_5th_dimension__protein_level_measurement_order"),
                 width = 20,
                 height = 15)

  #Biplot of Top10 variables and individuals ====
  message_function(text = "generating Biplot of Top10 variables and individuals ...",
                   color = "blue",
                   log_file_name = log_file_name)

  if(protein_data_PCA_res_only_1dim==TRUE){
    pca_plot_protein_ind_var <- ggplot()+
      scale_color_manual(values = condition_colors)+
      scale_fill_manual(values = condition_colors)+
      theme_classic(base_size = 18)+
      labs(title="Biplot of Top10 variables and individuals\n>>> calc. could not be performed because of N of replicates",
           subtitle="protein level",color="condition")+
      guides(fill = "none")

  }else{
    pca_plot_protein_ind_var <- factoextra::fviz_pca_biplot(protein_data_PCA_res,
                                                            label="var",
                                                            axes = c(1, 2),
                                                            mean.point=F,
                                                            pointsize=3,
                                                            pointshape=19,
                                                            ellipse.type = "convex",
                                                            ellipse.alpha=0.2,
                                                            addEllipses = T,
                                                            ellipse.level=0.95,
                                                            select.var = list(contrib = 10),
                                                            habillage = as.factor(habil$R.Condition))+
      scale_color_manual(values = condition_colors)+
      scale_fill_manual(values = condition_colors)+
      theme_classic(base_size = 18)+
      labs(title="Biplot of Top10 variables and individuals",subtitle="protein level",color="condition")+
      guides(fill = "none")
  }

  if(peptide_data_PCA_res_only_1dim==TRUE){
    pca_plot_peptide_ind_var <- ggplot()+
      scale_color_manual(values = condition_colors)+
      scale_fill_manual(values = condition_colors)+
      theme_classic(base_size = 18)+
      labs(title="Biplot of Top10 variables and individuals\n>>> calc. could not be performed because of N of replicates",
           subtitle="peptide level",color="condition")+
      guides(fill = "none")
  }else{
    pca_plot_peptide_ind_var <- factoextra::fviz_pca_biplot(peptide_data_PCA_res,
                                                            label="var",
                                                            axes = c(1, 2),
                                                            mean.point=F,
                                                            pointsize=3,
                                                            pointshape=19,
                                                            ellipse.type = "convex",
                                                            ellipse.alpha=0.2,
                                                            addEllipses = T,
                                                            ellipse.level=0.95,
                                                            select.var = list(contrib = 10),
                                                            habillage = as.factor(habil$R.Condition))+
      scale_color_manual(values = condition_colors)+
      scale_fill_manual(values = condition_colors)+
      theme_classic(base_size = 18)+
      labs(title="Biplot of Top10 variables and individuals",subtitle="peptide level",color="condition")+
      guides(fill = "none")
  }




  #individuals contribution plot  ====
  message_function(text = "generating individuals contribution plot ...",
                   color = "blue",
                   log_file_name = log_file_name)

  if(peptide_data_PCA_res_only_1dim==TRUE){
    contri_ind_peptide <- factoextra::get_pca_ind(peptide_data_PCA_res)$contrib
    contri_ind_peptide <- cbind(contri_ind_peptide,
                                    data.frame(Dim.2 = rep(0,nrow(contri_ind_peptide)),
                                               Dim.3 = rep(0,nrow(contri_ind_peptide)),
                                               Dim.4 = rep(0,nrow(contri_ind_peptide)),
                                               Dim.5 = rep(0,nrow(contri_ind_peptide))
                                                   )
    )
  }else{
    contri_ind_peptide <- factoextra::get_pca_ind(peptide_data_PCA_res)$contrib[,1:ifelse(ncol_peptide<=5,ncol_peptide,5)]
  }
  contri_ind_peptide <- contri_ind_peptide[unlist(sapply(habil$R.FileName,function(x)which(x==rownames(contri_ind_peptide)))),] #reordering
  contri_ind_peptide_plot<-ggcorrplot_edited(t(contri_ind_peptide),method="circle", colors = rev(c("#a50f15", "#fb6a4a", "grey")),
                                   legend.title = "contribution",leg.lim=c(min(contri_ind_peptide),max(contri_ind_peptide)))+
    theme_classic(base_size = 18)+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5,hjust = 1))+
    labs(title="contribution plot of individuals",subtitle="peptide level",y="",x="")

  if(protein_data_PCA_res_only_1dim==TRUE){
    contri_ind_protein <- factoextra::get_pca_ind(protein_data_PCA_res)$contrib
    contri_ind_protein <- cbind(contri_ind_protein,
                                           data.frame(Dim.2 = rep(0,nrow(contri_ind_protein)),
                                                      Dim.3 = rep(0,nrow(contri_ind_protein)),
                                                      Dim.4 = rep(0,nrow(contri_ind_protein)),
                                                      Dim.5 = rep(0,nrow(contri_ind_protein))
                                           )
    )
  }else{
    contri_ind_protein <- factoextra::get_pca_ind(protein_data_PCA_res)$contrib[,1:ifelse(ncol_protein<=5,ncol_protein,5)]

  }

  contri_ind_protein <- contri_ind_protein[unlist(sapply(habil$R.FileName,function(x)which(x==rownames(contri_ind_protein)))),] #reordering
  contri_ind_protein_plot<-ggcorrplot_edited(t(contri_ind_protein),method="circle", colors = rev(c("#a50f15", "#fb6a4a", "grey")),
                                      legend.title = "contribution",leg.lim=c(min(contri_ind_protein),max(contri_ind_protein)))+
    theme_classic(base_size = 18)+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5,hjust = 1))+
    labs(title="contribution plot of individuals",subtitle="protein level",y="",x="")

  #generate final individuals contribution plot & save
  contri_ind_final <- contri_ind_peptide_plot+contri_ind_protein_plot
  ggsave_pdf_png(filename = paste0(out_folder,"/","04_multivariate_analysis/",sample_length,"_sample_analysis/PCA_individuals_contribution_plot"),
         plot = contri_ind_final,
         width = 18,
         height = if((0.35*sample_length)<10){10}else{0.35*sample_length},
         limitsize = F)


  #aggregating all plots to a final plot matrix with PATCHWORK ====
  message_function(text = "generating PCA plot matrix ...",
                   color = "blue",
                   log_file_name = log_file_name)

  PCA_final_plots <- (peptide_data_PCA_res_screePlot + protein_data_PCA_res_screePlot)/
                    (peptide_data_PCA_res_biplot + protein_data_PCA_res_biplot)/
                    (pca_plot_peptide + pca_plot_protein)/
                    (pca_plot_peptide_ind_var + pca_plot_protein_ind_var)

  message_function(text = "saving PCA plot matrix ...",
                   color = "blue",
                   log_file_name = log_file_name)

  ggsave_pdf_png(filename = paste0(out_folder,"/","04_multivariate_analysis/",sample_length,"_sample_analysis/PCA_plots"),
                 plot = PCA_final_plots,
                 width = 25,
                 height = 25)

  #save 3D PCA plots  ====
  message_function(text = "saving 3D PCA plots ...",
                   color = "blue",
                   log_file_name = log_file_name)

  htmlwidgets::saveWidget(plotly3d_ind_plot_peptide,
                          selfcontained = T,
                          file = paste0(out_folder,"/","04_multivariate_analysis/",sample_length,"_sample_analysis/3D_PCA_plot__peptide_level.html"))
  htmlwidgets::saveWidget(plotly3d_ind_plot_protein,
                          selfcontained = T,
                          file = paste0(out_folder,"/","04_multivariate_analysis/",sample_length,"_sample_analysis/3D_PCA_plot__protein_level.html"))



  #performing HCPC analysis  ====
  if(HCPC_analysis == T){
    message_function(text = "performing HCPC analysis (this might take some time) ...",
                     color = "blue",
                     log_file_name = log_file_name)

    peptide_data_PCA_res_hcpc <- FactoMineR::HCPC(peptide_data_PCA_res,
                                                  graph = FALSE,
                                                  nb.clust=-1,
                                                  metric = "euclidean",
                                                  method = "complete")
    protein_data_PCA_res_hcpc <- FactoMineR::HCPC(protein_data_PCA_res,
                                                  graph = FALSE,
                                                  nb.clust=-1,
                                                  metric = "euclidean",
                                                  method = "complete")

    message_function(text = "saving HCPC plots ...",
                     color = "blue",
                     log_file_name = log_file_name)

    #save HCPC peptide level
    plot_save_png(filename = paste0(out_folder,"/","04_multivariate_analysis/",sample_length,"_sample_analysis/HCPC_plots_peptide_level__3D_map"),
                      plot = plot(peptide_data_PCA_res_hcpc,
                                  choice = "3D.map",
                                  title="HCPC analysis 3D plot (peptide level)"),
                      width = 15,
                      height = 15)


    plot_save_png(filename = paste0(out_folder,"/","04_multivariate_analysis/",sample_length,"_sample_analysis/HCPC_plots_peptide_level__2D_map"),
                      plot =   plot(peptide_data_PCA_res_hcpc,
                                    choice = "map",
                                    title="HCPC analysis tree plot (peptide level)"),
                      width = 15,
                      height = 15)

    plot_save_png(filename = paste0(out_folder,"/","04_multivariate_analysis/",sample_length,"_sample_analysis/HCPC_plots_peptide_level__barplot"),
                      plot =   plot(peptide_data_PCA_res_hcpc,
                                    choice = "bar",
                                    title="HCPC analysis tree plot (peptide level)"),
                      width = 15,
                      height = 15)

    plot_save_png(filename = paste0(out_folder,"/","04_multivariate_analysis/",sample_length,"_sample_analysis/HCPC_plots_peptide_level__tree_plot"),
                      plot =   plot(peptide_data_PCA_res_hcpc,
                                    choice = "tree",
                                    title="HCPC analysis tree plot (peptide level)"),
                      width = 15,
                      height = 15)

    #save HCPC protein level
    plot_save_png(filename = paste0(out_folder,"/","04_multivariate_analysis/",sample_length,"_sample_analysis/HCPC_plots_protein_level__3D_map"),
                      plot = plot(protein_data_PCA_res_hcpc,
                                  choice = "3D.map",
                                  title="HCPC analysis 3D plot (protein level)"),
                      width = 15,
                      height = 15)

    plot_save_png(filename = paste0(out_folder,"/","04_multivariate_analysis/",sample_length,"_sample_analysis/HCPC_plots_protein_level__2D_map"),
                      plot = plot(protein_data_PCA_res_hcpc,
                                  choice = "map",
                                  title="HCPC analysis tree plot (protein level)"),
                      width = 15,
                      height = 15)

    plot_save_png(filename = paste0(out_folder,"/","04_multivariate_analysis/",sample_length,"_sample_analysis/HCPC_plots_protein_level__barplot"),
                      plot = plot(protein_data_PCA_res_hcpc,
                                  choice = "bar",
                                  title="HCPC analysis tree plot (protein level)"),
                      width = 15,
                      height = 15)

    plot_save_png(filename = paste0(out_folder,"/","04_multivariate_analysis/",sample_length,"_sample_analysis/HCPC_plots_protein_level__tree_plot"),
                      plot = plot(protein_data_PCA_res_hcpc,
                                  choice = "tree",
                                  title="HCPC analysis tree plot (protein level)"),
                      width = 15,
                      height = 15)

  }



  #correlation plot ====
  message_function(text = "performing correlation analysis ...",
                   color = "blue",
                   log_file_name = log_file_name)

  corr_protein <- round(stats::cor(t(log2(protein_data_PCA)),
                            method = "spearman",
                            use = "complete.obs"), 2)
  corr_peptide <- round(stats::cor(t(log2(peptide_data_PCA)),
                            method = "spearman",
                            use = "complete.obs"), 2)

  message_function(text = "performing correlation analysis - calculate p-values ...",
                   color = "blue",
                   log_file_name = log_file_name)

  p_mat_protein <- ggcorrplot::cor_pmat(t(log2(protein_data_PCA)))
  p_mat_peptide <- ggcorrplot::cor_pmat(t(log2(peptide_data_PCA)))


  # using hierarchical clustering ====
  message_function(text = "generating correlation plots ...",
                   color = "blue",
                   log_file_name = log_file_name)

  correlation_plot_protein<-ggcorrplot_edited(corr_protein, colors = c("white", "#feb24c", "#de2d26"),
                               leg.lim = c(min(corr_protein),max(corr_protein)),
                               legend.title = "correlation score",
                               sig.level = 0.05,hc.order = T,hc.method = "complete",
                               method = "square",
                               type = "full",
                               lab = F,
                               p.mat = p_mat_protein,
                               insig = "blank")+
    theme_classic(base_size = 12)+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5,hjust = 1))+
    labs(title="correlation plot (protein level)",
         subtitle=paste("spearman corr. & reord. (linkage:","complete",")"),
         y="",
         x="",
         caption="insignificant correlation = blank squares")

  correlation_plot_peptide<-ggcorrplot_edited(corr_peptide, colors = c("white", "#feb24c", "#de2d26"),
                                       leg.lim = c(min(corr_peptide),max(corr_peptide)),
                                       legend.title = "correlation score",
                                       sig.level = 0.05,hc.order = T,hc.method = "complete",
                                       method = "square",
                                       type = "full",
                                       lab = F,
                                       p.mat = p_mat_peptide,
                                       insig = "blank")+
    theme_classic(base_size = 12)+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5,hjust = 1))+
    labs(title="correlation plot (peptide level)",
         subtitle=paste("spearman corr. & reord. (linkage:","complete",")"),
         y="",
         x="",
         caption="insignificant correlation = blank squares")


  correlation_plot_final <- correlation_plot_peptide+correlation_plot_protein

  message_function(text = "saving correlation plots ...",
                   color = "blue",
                   log_file_name = log_file_name)

  ggsave_pdf_png(filename = paste0(out_folder,"/","04_multivariate_analysis/",sample_length,"_sample_analysis/correlation_plots"),
         plot = correlation_plot_final,
         width = if((0.4*sample_length)<13){20}else{0.5*sample_length},
         height = if((0.4*sample_length)<10){10}else{0.25*sample_length},
         limitsize = F,
         dpi = 100)

  # correlation plot condition added ====

  # peptide level
  corr_peptide_tidy <- as.data.frame(corr_peptide)
  corr_peptide_tidy_rows <- rownames(corr_peptide_tidy)
  corr_peptide_tidy$rows <- corr_peptide_tidy_rows
  corr_peptide_tidy <- corr_peptide_tidy %>%
    pivot_longer(colnames(corr_peptide_tidy)[-which(colnames(corr_peptide_tidy)=="rows")],
                 names_to = "cols",
                 values_to = "values") %>%
    mutate(cols = factor(.data$cols,levels = corr_peptide_tidy_rows)) %>%
    mutate(rows = factor(.data$rows,levels = corr_peptide_tidy_rows))

  #add conditions
  corr_peptide_tidy <- left_join(corr_peptide_tidy,
                                 habil %>%
                                   dplyr::select(.data$R.FileName,
                                                 .data$R.Condition) %>%
                                   dplyr::rename(rows_cond = .data$R.Condition),
                                 by = c("rows" = "R.FileName"))
  corr_peptide_tidy <- left_join(corr_peptide_tidy,
                                 habil %>%
                                   dplyr::select(.data$R.FileName,
                                                 .data$R.Condition) %>%
                                   dplyr::rename(cols_cond = .data$R.Condition),
                                 by = c("cols" = "R.FileName"))

  correlation_plot_peptide_condition <-  ggplot(data = corr_peptide_tidy,
                                                mapping = aes(x = .data$rows,
                                                y = .data$cols,
                                                fill = .data$values))+
    geom_tile(color = "grey90")+
    theme_light()+
    scale_y_discrete(expand = c(0,0))+
    scale_x_discrete(expand = c(0,0))+
    scale_fill_gradient2(name = "correlation score",
                         low = "white",
                         mid = "#feb24c",
                         high = "#de2d26",
                         midpoint = mean(c(min(corr_peptide),max(corr_peptide))),
                         oob = scales::squish)+
    facet_grid(cols_cond ~ rows_cond,scales = "free", switch = "both")+
    theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1),
          panel.spacing = unit(0, "npc"),
          legend.position = "bottom",
          panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.ticks.y.left = element_blank(),
          legend.key.width = unit(2.5, "cm"),
          strip.text = element_text(size = 5),
          strip.background = element_rect(color = c("black")))+
    labs(title="correlation plot (peptide level)",
         subtitle = "spearman corr.",
         y="",
         x="",
         caption="")

  #write out peptide correlation data
  write_csv(x = corr_peptide_tidy %>%
              dplyr::distinct(.data$rows,
                              .data$cols,
                              .data$rows_cond,
                              .data$cols_cond,
                              .data$values) %>%
              dplyr::rename(R.FileName_1 = .data$rows,
                            R.FileName_2 = .data$cols,
                            R.Condition_1 = .data$rows_cond,
                            R.Condition_2 = .data$cols_cond,
                            correlation_score = .data$values),
            file = paste0(out_folder,"/","04_multivariate_analysis/",
                          sample_length,"_sample_analysis/correlation_spearman__peptide_level.csv"))

  # protein level
  corr_protein_tidy <- as.data.frame(corr_protein)
  corr_protein_tidy_rows <- rownames(corr_protein_tidy)
  corr_protein_tidy$rows <- corr_protein_tidy_rows
  corr_protein_tidy <- corr_protein_tidy %>%
    pivot_longer(colnames(corr_protein_tidy)[-which(colnames(corr_protein_tidy)=="rows")],
                 names_to = "cols",
                 values_to = "values") %>%
    mutate(cols = factor(.data$cols,levels = corr_protein_tidy_rows)) %>%
    mutate(rows = factor(.data$rows,levels = corr_protein_tidy_rows))

  #add conditions
  corr_protein_tidy <- left_join(corr_protein_tidy,
                                 habil %>%
                                   dplyr::select(.data$R.FileName, .data$R.Condition) %>%
                                   dplyr::rename(rows_cond = .data$R.Condition),
                                 by = c("rows" = "R.FileName"))
  corr_protein_tidy <- left_join(corr_protein_tidy,
                                 habil %>%
                                   dplyr::select(.data$R.FileName, .data$R.Condition) %>%
                                   dplyr::rename(cols_cond = .data$R.Condition),
                                 by = c("cols" = "R.FileName"))

  correlation_plot_protein_condition <-  ggplot(data = corr_protein_tidy,mapping = aes(x = .data$rows,
                                                                                       y = .data$cols,
                                                                                       fill = .data$values))+
    geom_tile()+
    theme_light()+
    scale_y_discrete(expand = c(0,0))+
    scale_x_discrete(expand = c(0,0))+
    scale_fill_gradient2(name = "correlation score",
                         low = "white",
                         mid = "#feb24c",
                         high = "#de2d26",
                         midpoint = mean(c(min(corr_protein),max(corr_protein))),
                         oob = scales::squish)+
    facet_grid(cols_cond ~ rows_cond,scales = "free", switch = "both")+
    theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1),
          panel.spacing = unit(0, "npc"),
          axis.ticks.y.left = element_blank(),
          panel.grid = element_blank(),
          panel.border = element_blank(),
          legend.position = "bottom",
          legend.key.width = unit(2.5, "cm"),
          strip.text = element_text(size = 5),
          strip.background = element_rect(color = c("black")))+
    labs(title="correlation plot (protein level)",
         subtitle = "spearman corr.",
         y="",
         x="",
         caption="")

  #write out peptide correlation data
  write_csv(x = corr_protein_tidy %>%
              dplyr::distinct(.data$rows,
                              .data$cols,
                              .data$rows_cond,
                              .data$cols_cond,
                              .data$values) %>%
              dplyr::rename(R.FileName_1 = .data$rows,
                            R.FileName_2 = .data$cols,
                            R.Condition_1 = .data$rows_cond,
                            R.Condition_2 = .data$cols_cond,
                            correlation_score = .data$values),
            file = paste0(out_folder,"/","04_multivariate_analysis/",
                          sample_length,"_sample_analysis/correlation_spearman__protein_level.csv"))

  correlation_plot_final_condition <- correlation_plot_peptide_condition+correlation_plot_protein_condition

  message_function(text = "saving correlation plots with condition info ...",
                   color = "blue",
                   log_file_name = log_file_name)

  ggsave_pdf_png(filename = paste0(out_folder,"/","04_multivariate_analysis/",sample_length,"_sample_analysis/correlation_plots_condition"),
                 plot = correlation_plot_final_condition,
                 width = if((0.4*sample_length)<13){20}else{0.4*sample_length},
                 height = if((0.4*sample_length)<10){10}else{0.2*sample_length},
                 limitsize = F,
                 dpi = 100)

  #umap analysis on protein level ====
  message_function(text = "perform UMAP analysis on protein level ...",
                   color = "blue",
                   log_file_name = log_file_name)

  umap_data<- SpectroPipeR_data_quant$protein_data %>%
    dplyr::select(.data$R.FileName,
                  .data$PG.ProteinGroups,
                  .data$protein_intensity) %>%
    pivot_wider(names_from = .data$PG.ProteinGroups,
                values_from = .data$protein_intensity)

  # umap config
  custom.settings = umap::umap.defaults
  custom.settings$n_neighbors = 5

  protein_data_PCA_umap <- c()

  message_function(text = "try UMAP with default settings ...",
                   color = "blue",
                   log_file_name = log_file_name)

  protein_data_PCA_umap <- try(umap::umap(log2(protein_data_PCA),
                                          config = umap::umap.defaults,
                                          method = "naive"),
                               silent = T)
  #if umap fails due to low complex data
  if(as.character(class(protein_data_PCA_umap))=="try-error"){

    message_function(text = "due to low complexity try UMAP with n_neighbors = 5 setting ...",
                     color = "blue",
                     log_file_name = log_file_name)

    protein_data_PCA_umap <- try(umap::umap(log2(protein_data_PCA),
                                            config = custom.settings,
                                            method = "naive"),
                                 silent = T)
  }
  #if umap still fails
  if(as.character(class(protein_data_PCA_umap))=="try-error"){
      message_function(text = "UMAP failed ...",
                     color = "blue",
                     log_file_name = log_file_name)
  }else{
    write_rds(x = protein_data_PCA_umap,
              file = paste0(out_folder,"/","04_multivariate_analysis/",sample_length,"_sample_analysis/UMAP_analysis__protein_level.RDS"))

    protein_data_PCA_umap_extract <- as.data.frame(protein_data_PCA_umap$layout)
    colnames(protein_data_PCA_umap_extract) <- c("x","y")
    protein_data_PCA_umap_extract$R.FileName <- rownames(protein_data_PCA_umap$layout)
    protein_data_PCA_umap_extract <- left_join(protein_data_PCA_umap_extract,
                                               habil,
                                               by = "R.FileName")

    # generate the convex hull for umap
    protein_data_PCA_umap_extract_hull_group <- protein_data_PCA_umap_extract %>%
      group_by(.data$R.Condition) %>%
      slice(chull(.data$x, .data$y))

    # slice --> dplyr function
    # chull --> grDevice function
    ## Computes the subset of points which lie on the convex hull of the set of points specified.


    # Define the scatterplot
    umap_protein_plot <- ggplot(protein_data_PCA_umap_extract, aes(x = .data$x,
                                                                   y = .data$y,
                                                                   color = .data$R.Condition,
                                                                   fill = .data$R.Condition)) +
      geom_hline(yintercept = 0,linetype = "dashed")+
      geom_vline(xintercept = 0,linetype = "dashed")+
      theme_classic(base_size = 16)+
      geom_point(size=3)+ # add point shape R.Replicate
      geom_polygon(data = protein_data_PCA_umap_extract_hull_group, alpha = 0.5)+ # add convex shape
      scale_color_manual(values = condition_colors)+
      scale_fill_manual(values = condition_colors)+
      labs(title = "UMAP (protein level)",
           subtitle = "protein intensities",
           caption = "",
           x = "UMAP 1",# labeling of axis
           y = "UMAP 2",
           # get expl. variance  Dim.xxx data by pca$eig[xxx,2]
           color = "Conditions",
           fill = "Conditions")+
      guides(color="none")+
      ggrepel::geom_label_repel(data = protein_data_PCA_umap_extract %>%
                         group_by(.data$R.Condition) %>%
                         summarise(mean_x = mean(.data$x),
                                   mean_y = mean(.data$y)) %>%
                         ungroup(),
                       mapping = aes(x = .data$mean_x,
                                     y = .data$mean_y,
                                     label = .data$R.Condition),
                       color="white",
                       segment.colour = "black",
                       size = 4,
                       alpha = 0.8,
                       min.segment.length = 0.01,
                       box.padding = 1,show.legend = FALSE)

    ggsave_pdf_png(filename = paste0(out_folder,"/","04_multivariate_analysis/",sample_length,"_sample_analysis/UMAP_plot_protein_level_conditions_marked"),
                   plot = umap_protein_plot,
                   width = 12,
                   height = 9)

  }

# generate outputs --------------------------------------------------------


  output_list <- list(PCA_peptide_intensity = peptide_data_PCA_res,
                      PCA_protein_intensity = protein_data_PCA_res,
                      UMAP_protein_intensity = protein_data_PCA_umap,
                      peptide_intensity_correlation = corr_peptide,
                      protein_intensity_correlation = corr_protein
                      )

  message_function(text = paste0("MVA module done --> please check outputs in folder: ", out_folder,"/","04_multivariate_analysis/"),
                   color = "blue",
                   log_file_name = log_file_name)

  class(output_list) <- "SpectroPipeR_data_MVA"
  return(output_list)




}
