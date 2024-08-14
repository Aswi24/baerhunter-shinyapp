#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_server <- function(input, output, session) {
  # Initialize and use the ncRNA annotation module server-side functionality
  results <- mod_annot_ncRNA_server("annot_ncRNA_1")

  # Render the input UI components for ncRNA annotation using a custom module
  output$input_ui <- renderUI({
    mod_annot_ncrna_ui("annot_ncRNA_1")$input_ui
  })

  # Render the output UI components for ncRNA annotation using a custom module
  output$output_ui <- renderUI({
    mod_annot_ncrna_ui("annot_ncRNA_1")$output_ui
  })

  # Initialize the UI for preprocessing BAM files
  preprocessing_ui <- mod_preprocessing_ui("preprocessing_1")

  # Render the input UI components for BAM file preprocessing
  output$preprocessing_input_ui <- renderUI({
    preprocessing_ui$input_prepro_ui
  })

  # Render the output UI components for BAM file preprocessing
  output$preprocessing_output_ui <- renderUI({
    preprocessing_ui$output_prepro_ui
  })

  # Initialize and use the preprocessing module
  mod_preprocessing_server("preprocessing_1")

  # Open a modal dialog with a genome browser when an action button is clicked
  observeEvent(results$action_btn(), {
    showModal(modalDialog(
      title = "Genome Browser",
      size = "l",
      easyClose = TRUE,
      footer = NULL,
      genomeBrowserUI("genome_browser_modal")
    ))
    # Call to server function of the genome browser with table data
    genomeBrowserServer("genome_browser_modal", results$table_df)
  })
}
