#'  ncRNA Annotation Processing Module
#'
#' @description  User interface that enables users to upload BAM files and GFF3 annotation files,
#' and outputs the summary tables
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList

options(shiny.maxRequestSize = 50 * 1024^3)
mod_annot_ncrna_ui <- function(id) {
  ns <- NS(id)
  list(
    input_ui = tagList(
      tags$h3("Input:"),
      fileInput(ns("bam"),
                "Choose BAM files",
                multiple=TRUE,
                accept=c('.bam'),
                buttonLabel = "Browse"),
      fileInput(ns("original_annotation_file"),
                "Upload annotation file",
                accept=c('.gff3'),
                buttonLabel = "Browse"),
      selectInput(ns("original_sRNA_annotation"), "Enter pre-annotated sRNA if known", choices = c("unknown", "Enter Manually")),
      # panel appears only when unknown is selected
      tabsetPanel(
        id = ns("sRNA_biotype"),
        type = "hidden",
        tabPanelBody("unknown"),
        tabPanelBody("Enter Manually", textInput(ns("manual_entry"), "Enter data manually:"))
      ),
      textInput(ns("refseq_name"), "Enter name of the reference sequence"),
      selectInput(ns("paired_end_data"), "Paired end data", choices = c(TRUE, FALSE)),
      selectInput(ns("strandedness"), "Select Strandedness", choices = c("stranded", "reversely_stranded", "unstranded")),
      numericInput(ns("low_coverage_cutoff"), "Low coverage cutoff", value = 1),
      numericInput(ns("high_coverage_cutoff"), "High coverage cutoff", value = 1),
      numericInput(ns("peak_width"), "Minimum sRNA length", value = 1),
      numericInput(ns("min_UTR_length"), "Minimum UTR length", value = 1),
      actionButton(ns("submitbutton"), "Submit", class = "btn btn-primary")
    ),

    output_ui = tagList(

      tags$label(tags$h3('Output')),
      DT::DTOutput(ns("table1")),
      br(),
      uiOutput(ns("action_button_ui")),
      br(),
      DT::DTOutput(ns("table2")),
      br(),
      tableOutput(ns("table3")),
      br(),
      DT::DTOutput(ns("table4")),
      br(),
      downloadButton(ns("download_gff"), "Download GFF file")
    )
  )
}


mod_annot_ncRNA_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    #reactive table for detected ncRNAs
    table_df <- reactiveVal()

    #view genome browser to be shown only after output of the table
    show_action_button <- reactiveVal(FALSE)


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

      if (!(input$low_coverage_cutoff > 0 &&
            round(input$low_coverage_cutoff) == input$low_coverage_cutoff)) {
        showNotification("Low coverage cutoff should be a positive integer greater than zero.", type = "error")
        return() # Stop further processing
      }

      if (!(input$high_coverage_cutoff > 0 &&
            round(input$high_coverage_cutoff) == input$high_coverage_cutoff)) {
        showNotification("High coverage cutoff should be a positive integer greater than zero.", type = "error")
        return() # Stop further processing
      }

      if (!(input$peak_width > 0 &&
            round(input$peak_width) == input$peak_width)) {
        showNotification("Minimum sRNA length should be a positive integer greater than zero.", type = "error")
        return() # Stop further processing
      }

      if (!(input$peak_width > 0 &&
            round(input$peak_width) == input$peak_width)) {
        showNotification("Minimum sRNA length should be a positive integer greater than zero.", type = "error")
        return() # Stop further processing
      }

      if (!(input$min_UTR_length > 0 &&
            round(input$min_UTR_length) == input$min_UTR_length)) {
        showNotification("Minimum UTR length should be a positive integer greater than zero.", type = "error")
        return() # Stop further processing
      }

      bam_files <- input$bam$datapath
      paired_end_data <- input$paired_end_data
      low_coverage_cutoff <- input$low_coverage_cutoff
      high_coverage_cutoff <- input$high_coverage_cutoff
      peak_width <- input$peak_width
      strandedness <- input$strandedness
      min_UTR_length <- input$min_UTR_length
      original_annotation_file <- input$original_annotation_file$datapath
      original_sRNA_annotation <- ncRNA_term


      refseq_name <- input$refseq_name


      tmpfile <- tempfile(fileext = ".gff3")

      #progress messages while the app is running

      withProgress(message = 'Starting...', value = 0, {
        result <- feature_file_editor(bam_files, original_annotation_file,original_sRNA_annotation, tmpfile,  low_coverage_cutoff, high_coverage_cutoff, peak_width, min_UTR_length, paired_end_data = FALSE, strandedness  = "stranded",
                                      progressCallback = function(value, message) {
                                        setProgress(value = value, message = message)
                                      }
        )

      })


      summary <- ncrna_annot_compare(original_annotation_file, tmpfile, refseq_name)

      # populates the table
      table_df(as.data.frame(summary[[1]]))


      output$download_gff <- downloadHandler(
        filename = function() {
          paste0("New annotation file", ".gff3")
        },
        content = function(file) {
          file.copy(tmpfile, file)
        }
      )

      output$table1 <- DT::renderDT({
        data <- table_df()  # Get the data for table1
        if (!is.null(data)) {
          show_action_button(TRUE)  # Show the action button when table1 has data
          DT::datatable(data)
        }
      })

      output$action_button_ui <- renderUI({
        if (show_action_button()) {
          actionButton(ns("action_button_ui"), "View in browser")
        }
      })

      # Logic to activate the module when the button is clicked
      shinyjs::click(id = "action_button_ui")


      output$table2 <- DT::renderDT({
        summary[[2]]
      })

      output$table3 <- renderTable({
        summary[[3]]
      })

      output$table4 <- DT::renderDT({
        summary[[5]]
      })

    })
    return(list(
      action_btn = reactive(input$action_button_ui),
      table_df = table_df  # Return the reactive expression itself
    ))
  })
}

