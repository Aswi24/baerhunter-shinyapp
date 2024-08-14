#' preprocessing UI Function
#'
#' @description user interface that allows input BAM files and annotation files to get information and coverage plot
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
#' '
#'
#'
mod_preprocessing_ui <- function(id) {
  ns <- NS(id)
  list(
    input_prepro_ui = tagList(
      tags$h3("Input:"),
      fileInput(ns("bam"),
                "Choose BAM files",
                multiple = TRUE,
                accept = c('.bam'),
                buttonLabel = "Browse"),
      #actionButton(ns("get_stats"), "Get Stats", class = "btn btn-primary"),
      fileInput(ns("original_annotation_file"),
                "Upload Annotation File",
                accept = c('.gff3'),
                buttonLabel = "Browse"),
      selectInput(ns("original_sRNA_annotation"), "Enter pre-annotated sRNA if known", choices = c("unknown", "Enter Manually")),
      tabsetPanel(
        id = ns("sRNA_biotype"),
        type = "hidden",
        tabPanel("unknown"),
        tabPanel("Enter Manually", textInput(ns("manual_entry"), "Enter data manually:"))
      ),
      textInput(ns("refseq_name"), "Enter name of the reference sequence"),
      selectInput(ns("paired_end_data"), "Paired end data", choices = c(TRUE, FALSE)),
      selectInput(ns("strandedness"), "Select Strandedness", choices = c("stranded", "reversely_stranded", "unstranded")),
      actionButton(ns("submitbutton"), "View coverage plot", class = "btn btn-primary")
    ),
    output_prepro_ui = tagList(
      tableOutput(ns("percentile_table")),
      plotOutput(ns("coverage_plot"))
    )
  )
}
#' preprocessing Server Functions
#'
#' @noRd
mod_preprocessing_server <- function(id){
  moduleServer(id, function(input, output, session){
    ns <- session$ns

    observeEvent(input$original_sRNA_annotation, {
      updateTabsetPanel(inputId = "sRNA_biotype", selected = input$original_sRNA_annotation)
    })

    observeEvent(input$submitbutton, {

      if (input$original_sRNA_annotation == "Enter Manually" && str_trim(input$manual_entry) == "") {
        showNotification("Please enter the data manually as selected.", type = "error")
        return() # Stop further processing
      }

      ncRNA_term <- input$original_sRNA_annotation

      if (ncRNA_term == "Enter Manually") {
        ncRNA_term <- input$manual_entry
      }


      if (is.null(input$bam) || tools::file_ext(input$bam$name) != "bam") {
        showNotification("Please upload a valid BAM file.", type = "error")
        return() # Stop further processing
      }

      if (is.null(input$original_annotation_file) || tools::file_ext(input$original_annotation_file$name) != "gff3") {
        showNotification("Please upload a valid GFF3 file.", type = "error")
        return() # Stop further processing
      }
      bam_files <- input$bam$datapath
      bam_names <- input$bam$name
      paired_end_data <- input$paired_end_data
      original_annotation_file <- input$original_annotation_file$datapath
      original_sRNA_annotation <- ncRNA_term
      refseq_name <- input$refseq_name
      original_annotation_file <- input$original_annotation_file$datapath
      gaps <- igr_strand_specific(original_annotation_file, original_sRNA_annotation )
      gaps_plus <- gaps$gaps_plus
      gaps_minus <- gaps$gaps_minus

      results <- Map(function(bam_file, bam_name) {
        process_bam_file(bam_file = bam_file, bam_name = bam_name, gaps_plus = gaps_plus, gaps_minus = gaps_minus)
      }, bam_files, bam_names)

      # Extract and combine all coverage data frames
      combined_coverage_df <- bind_rows(lapply(results, `[[`, "coverage_df"))

      # Extract and combine all percentiles data frames
      combined_percentiles_df <- bind_rows(lapply(results, `[[`, "percentiles_df"))

      output$coverage_plot <- renderPlot({

        ggplot2::ggplot(combined_coverage_df, ggplot2::aes(x = File, y = log2(Coverage), fill = File)) +
          PupillometryR::geom_flat_violin(position = ggplot2::position_nudge(x = .25, y = 0), adjust = 2, trim = FALSE) +
          ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.3, width = .1, colour = "BLACK") +
          ggplot2::ylab('log2(Coverage)') + ggplot2::xlab('Files') + ggplot2::coord_flip() +
          cowplot::theme_cowplot() + ggplot2::guides(fill = FALSE)
      })

      output$percentile_table <- renderTable({
        combined_percentiles_df
      })


    })
  })
}
## To be copied in the UI
# mod_preprocessing_ui("preprocessing_1")

## To be copied in the server
# mod_preprocessing_server("preprocessing_1")
