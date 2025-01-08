
#' SpectroPipeR UI function
#' @description
#' launches a user interface in the browser, providing access to the majority of SpectroPipeR's functionalities.
#' @return
#' launches a user interface in the browser
#' @export
#' @import shiny
#' @importFrom bslib bs_theme
#' @import tidyverse
#' @importFrom future multisession future
#' @import parallel
#' @examples
#' \donttest{
#' SpectroPipeR_ui()
#'}
#'
SpectroPipeR_ui<- function(){

  future::plan(future::multisession)
  # Define the UI
  ui <- shiny::navbarPage(
    title = shiny::tagList(shiny::tags$img(src = "https://github.com/stemicha/SpectroPipeR/raw/main/vignettes/figures/SpectroPipeR_hexbin_logo.png",
                             height = "80px"),
                    "SpectroPipeR"),
    theme = bslib::bs_theme(
      bootswatch = "flatly"
    ),

    # color of buttons --------------------------------------------------------
    header = shiny::tags$style(shiny::HTML("
    #add_comparison { background-color: #005493; color: white; border: none; }
    #run { background-color: #28a745; color: white; border: none; }
    #restart { background-color: #dc3545; color: white; border: none; }
  ")),
    shiny::tabPanel(
      title = "GUI",
      shiny::sidebarLayout(

        # side bar panel ----------------------------------------------------------
        shiny::sidebarPanel(

          # add parameters ----------------------------------------------------------
          shiny::textInput(inputId = "SN_SpectroPipeR_file",
                    label = "paste absolute path to Spectronaut report in SpectroPipeR report layout (.tsv)",
                    placeholder = "/Users/stephanmichalik/Spectronaut_SpectroPipeR_report_file.tsv"),
          shiny::textInput(inputId = "output_folder",
                    label = "paste absolute path to output folder",
                    placeholder = "/Users/stephanmichalik/SpectroPipeR_output_folder"),
          shiny::h5("general parameter:"),
          shiny::fluidRow(
            shiny::column(width = 6,
                          shiny::numericInput(inputId = "ion_q_value_cutoff",
                                label = "Spectronaut ion Q-value cutoff",
                                value = 0.01)),
            shiny::column(width = 6,
                          shiny::numericInput(inputId = "max_chars_file_name_capping",
                                label = "character number for raw file name capping",
                                value = 25))
          ),

          shiny::fluidRow(
            shiny::column(width = 6,
                          shiny::selectInput(inputId = "ID_condition_filtering",
                               label = "condition wise ion filtering",
                               choices = c(TRUE,FALSE),
                               selected = FALSE)
            ),
            shiny::column(width = 6,
                          shiny::numericInput(inputId = "ID_condition_filtering_percent",
                                label = "condition wise ion filtering proportion",
                                min = 0,
                                max = 1,
                                value = 0.5)
            )
          ),
          shiny::fluidRow(
            shiny::column(width = 4,
                          shiny::selectInput(inputId = "filter_oxidized_peptides",
                               label = "filter oxidized peptides",
                               choices = c(TRUE,FALSE),
                               selected = TRUE)
            ),
            shiny::column(width = 4,
                          shiny::numericInput(inputId = "id_drop_cutoff",
                                label = "ion ID drop threshold",
                                value = 0.3)
            ),
            shiny::column(width = 4,
                   shiny::numericInput(inputId = "normalization_factor_cutoff_outlier",
                                label = "norm. factor threshold",
                                value = 4)
            )
          ),
          shiny::fluidRow(
            shiny::column(width = 6,
                          shiny::selectInput(inputId = "build_HTML_report",
                               label = "build HTML report file",
                               choices = c(TRUE,FALSE),
                               selected = FALSE)
            ),
            shiny::column(width = 6,
                          shiny::selectInput(inputId = "protein_intensity_estimation",
                               label = "protein intensity aggregation method",
                               choices = c("MaxLFQ","Hi3"),
                               selected = "MaxLFQ")
            )
          ),





          # condition comparison selection
          shiny::hr(),
          shiny::h4("selection of comparisons for statistics"),
          shiny::fileInput("condition_file",
                    "first ... upload Spectronaut condition setup file (*_ConditionSetup.tsv)",
                    accept = ".tsv"),
          shiny::uiOutput("condition_selector"),
          shiny::actionButton("add_comparison", "Add Comparison"),
          shiny::br(),
          shiny::uiOutput("comparison_list"),
          shiny::br(),
          shiny::fluidRow(
            shiny::column(width = 4,
                          shiny::selectInput(inputId = "stat_test",
                               label = "statistical test method used for statistical testing",
                               choices = c("rots","modt"),
                               selected = "rots")
            ),
            shiny::column(width = 4,
                          shiny::selectInput(inputId = "type_slr",
                               label = "stat. analysis peptide ratio aggregation to protein level",
                               choices = c("median","tukey"),
                               selected = "median")
            ),
            shiny::column(width = 4,
                          shiny::selectInput(inputId = "paired",
                               label = "should a paired statistical analysis be conducted?",
                               choices = c(TRUE,FALSE),
                               selected = FALSE)
            )
          ),

          shiny::numericInput(inputId = "fold_change",
                       label = "fold-change threshold for Volcano plots",
                       value = 1.5),
          shiny::numericInput(inputId = "p_value_cutoff",
                       label = "p- and q-value threshold for Volcano plots",
                       value = 0.05),
          shiny::uiOutput("number_of_cores_UI"),
          shiny::hr(),
          shiny::actionButton("run", "run SpectroPipeR",width = "100%",icon = shiny::icon("play")),
          shiny::br(),
          shiny::br(),
          shiny::actionButton("restart","restart Sessions",icon = shiny::icon("power-off"),width = "100%")
          #,
          #style = "background-color: #f2f2f2;" # background
        ),
        shiny::mainPanel(
          shiny::h3("After pressing 'run SpectroPipeR' the data data processing starts! Please be patient..."),
          shiny::verbatimTextOutput("log_output")
          #style = "background-color: #f2f2f2;" # background
        )
      )
    )
  )

  # Define the server logic
  server <- function(input, output, session) {



    # cores detect render UI ---------------------------------------------------

    output$number_of_cores_UI<- shiny::renderUI({
      shiny::sliderInput(inputId = "number_of_cores_statistics",
                  label = "number of cores used for stat. analysis",
                  min = 2,
                  max = parallel::detectCores()-2,
                  value = 2)
    })




    # server: restart session -------------------------------------------------
    shiny::observeEvent(input$restart, {
      # reload
      session$reload()
    })


    # server: condition selection  --------------------------------------------
    comparisons <- shiny::reactiveVal(list())
    # Function to get unique conditions from the uploaded file
    get_conditions <- function(df) {
      if (!is.null(df) && "Condition" %in% colnames(df)) {
        unique(df$Condition)
      } else {
        NULL
      }
    }

    # Reactive expression to read the uploaded file
    userFileCondFile <- shiny::reactive({
      shiny::req(input$condition_file)
      readr::read_delim(input$condition_file$datapath,delim = "\t")
    })

    # Generate UI for condition selection
    output$condition_selector <- shiny::renderUI({
      shiny::req(userFileCondFile())
      conditions <- get_conditions(userFileCondFile())
      if (!is.null(conditions)) {
        shiny::selectInput("conditions", "Select Conditions for Comparison",
                    choices = conditions, multiple = TRUE, selectize = TRUE)
      }
    })

    # Add the selected comparison to the list
    shiny::observeEvent(input$add_comparison, {
      shiny::req(length(input$conditions) == 2)
      current_comparisons <- comparisons()
      new_comparison <- list(input$conditions)
      comparisons(append(current_comparisons, new_comparison))
    })

    # Display the list of comparisons
    output$comparison_list <- shiny::renderUI({
      shiny::req(comparisons())
      shiny::tags$ul(
        lapply(comparisons(), function(comp) {
          shiny::tags$li(paste(comp, collapse = " vs. "))
        })
      )
    })

    comparisons_SpectroPipeR <- shiny::reactive({
      if(length(comparisons())==0){
        NULL
      }else{
        do.call(cbind,comparisons())
      }
    })


    # server: run SpectroPipeR ------------------------------------------------

    shiny::observeEvent(input$run, {
      req(input$SN_SpectroPipeR_file, input$output_folder);

      # Clear previous log output
      output$log_output <- renderText({ "Starting SpectroPipeR..." })

      # Run SpectroPipeR in other session
      future::future({
        SpectroPipeR_analysis <- SpectroPipeR::SpectroPipeR(file = shiny::isolate(gsub("\\\\", "/", input$SN_SpectroPipeR_file)),
                                              parameter = list(output_folder = shiny::isolate(gsub("\\\\", "/", input$output_folder)),
                                                               ion_q_value_cutoff = shiny::isolate(as.numeric(input$ion_q_value_cutoff)),
                                                               id_drop_cutoff = shiny::isolate(as.numeric(input$id_drop_cutoff)),
                                                               normalization_factor_cutoff_outlier  = shiny::isolate(as.numeric(input$normalization_factor_cutoff_outlier)),
                                                               filter_oxidized_peptides = shiny::isolate(as.logical(input$filter_oxidized_peptides)),
                                                               protein_intensity_estimation = shiny::isolate(as.character(input$protein_intensity_estimation)),
                                                               stat_test = shiny::isolate(as.character(input$stat_test)),
                                                               type_slr = shiny::isolate(as.character(input$type_slr)),
                                                               fold_change = shiny::isolate(as.numeric(input$fold_change)),
                                                               p_value_cutoff = shiny::isolate(as.numeric(input$p_value_cutoff)),
                                                               paired = shiny::isolate(as.logical(input$paired))
                                              ),
                                              max_chars_file_name_capping = shiny::isolate(as.numeric(input$max_chars_file_name_capping)),
                                              ID_condition_filtering = shiny::isolate(as.logical(input$ID_condition_filtering)),
                                              ID_condition_filtering_percent = shiny::isolate(as.numeric(input$ID_condition_filtering_percent)),
                                              condition_comparisons = shiny::isolate(comparisons_SpectroPipeR()),
                                              number_of_cores_statistics = shiny::isolate(as.numeric(input$number_of_cores_statistics)),
                                              build_HTML_report = shiny::isolate(as.logical(input$build_HTML_report))
        )

      })# end future

      # Monitor the output folder for the newest log file
      output$log_output <- shiny::renderText({
        shiny::invalidateLater(500, session)
        log_files <- list.files(gsub("\\\\", "/", input$output_folder), pattern = "\\.log$", full.names = TRUE)
        if (length(log_files) > 0) {
          newest_log <- log_files[which.max(file.info(log_files)$mtime)]
          log_content <- readLines(newest_log)
          paste(log_content, collapse = "\n")
        } else {
          "No log files found yet."
        }
      })

    }) # end run observerEvent

  } # end server function

  # Run the application
  shiny::shinyApp(ui = ui,
                 server = server,
                 options = list(launch.browser = TRUE))



}
