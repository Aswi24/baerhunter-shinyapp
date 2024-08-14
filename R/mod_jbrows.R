#' Genome Browser User Interface
#'
#' @description User interface that enables users to view the detected ncRNAs in an interactive
#' genome browser. Requires the path to the directory containing the sequence and annotation files.
#' And requires input of the reference sequence name as mentioned in the FASTA file
#' @param
#'
#'
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#'
#' #'
#' @importFrom shiny NS tagList
#' @noRd




genomeBrowserUI <- function(id) {
  ns <- NS(id)
  tagList(
    titlePanel("Genome browser"),
    textInput(ns("user_path"), "Enter the path to your FASTA files"),
    textInput(ns("ref_seq"), "Enter reference sequence name"),
    actionButton(ns("submit"), "Submit"),
    uiOutput(ns("select_assembly")),
    uiOutput(ns("select_feature")),
    DT::DTOutput(ns("table_browsr")),
    JBrowseR::JBrowseROutput(ns("browser_output"))
  )
}


genomeBrowserServer <- function(id,table_df) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    # Use ns for namespacing IDs
    location_str_rv <- reactiveVal()
    # Use a standard R variable to keep track of the data server
    active_data_server <- NULL

    urls <- reactiveVal(NULL)

    # Function to stop the data server
    stopDataServer <- function() {
      if (!is.null(active_data_server)) {
        active_data_server$stop_server()
        active_data_server <<- NULL
      }
    }

    observeEvent(input$submit, {
      req(input$user_path)

      # Stop existing data server before starting a new one
      stopDataServer()

      # Start the server and generate URLs
      active_data_server <<- serve_data2(input$user_path, 5000)
      url_data <- generate_cli_message(input$user_path, 5000)
      urls(list(
        assembly = setNames(url_data$assembly, basename(url_data$assembly)),
        feature = setNames(url_data$feature, basename(url_data$feature))
      ))
    })

    # Register a callback to stop the data server when the app stops
    onStop(function() {
      stopDataServer()
    })


    jbrowse_location <- reactiveVal()

    output$table_browsr <- DT::renderDT({
      req(table_df())  # Ensure table_df is available
      table_df()  # Return the data frame for rendering
    })

    # Update JBrowseR location based on sRNA table row selection
    observeEvent(input$table_browsr_rows_selected, {
      req(input$table_browsr_rows_selected, table_df())
      selected_data <- table_df()[input$table_browsr_rows_selected, ]

      # Ensure that the selected_data has the necessary columns
      if (ncol(selected_data) > 0 && all(c("start", "stop") %in% names(selected_data))) {
        # Construct location string using the selected row data
        location_string <- paste0(input$ref_seq, ":", selected_data$start[1], "..", selected_data$stop[1])

        # Update the reactive value with the location string
        location_str_rv(location_string)

        # Correctly update jbrowse_location with the value of location_str_rv
        jbrowse_location(location_str_rv())  # Use parentheses to get the value
      }
    })
    # Update JBrowse location
    #location_str_updated <- "AL123456.3:7205..7254"

    output$select_assembly <- renderUI({
      req(urls())  # Ensure that urls are not NULL
      selectInput(ns("selected_assembly"), "Select Assembly", choices = urls()$assembly)
    })

    output$select_feature <- renderUI({
      req(urls())  # Ensure that urls are not NULL
      selectInput(ns("selected_feature"), "Select Feature", choices = urls()$feature)
    })

    # Use the selected assembly and feature in JBrowseR
    output$browser_output <- JBrowseR::renderJBrowseR({
      req(input$selected_assembly, input$selected_feature, jbrowse_location())
      assembly_url <- input$selected_assembly
      feature_url <- input$selected_feature

      assembly <- JBrowseR::assembly(assembly_url, bgzip = TRUE)
      annotations_track <- JBrowseR::track_feature(feature_url, assembly)
      tracks <- JBrowseR::tracks(annotations_track)
      defaultSession <- JBrowseR::default_session(assembly, c(annotations_track))
      JBrowseR::JBrowseR("View", assembly = assembly, tracks = tracks, location = jbrowse_location(), defaultSession = defaultSession, display_assembly = FALSE)
    })
  })

}
