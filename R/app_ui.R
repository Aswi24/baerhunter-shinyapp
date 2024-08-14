#' The application User-Interface
#'
#' @param request Internal parameter for `{shiny}`.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_ui <- function(request) {
  tagList(
    # Leave this function for adding external resources
    golem_add_external_resources(),
    # Your application UI logic
    fluidPage(
      shinyjs::useShinyjs(),
      theme = shinythemes::shinytheme("flatly"),
      navbarPage(
        "Baerhunter",
        # First tab: Annotation of ncRNA
        tabPanel("Annotation of ncRNA",
                 sidebarLayout(
                   sidebarPanel(
                     # Dynamic UI for preprocessing inputs, to be rendered by the server
                     uiOutput("input_ui")
                   ),
                   mainPanel(
                     # Dynamic UI for preprocessing outputs, to be rendered by the server
                     uiOutput("output_ui")
                   )
                 )
        ),
        # Second tab: Preprocessing of BAM files
        tabPanel("Preprocessing of BAM files",
                 sidebarLayout(
                   sidebarPanel(
                     uiOutput("preprocessing_input_ui")  # For input elements
                   ),
                   mainPanel(
                     uiOutput("preprocessing_output_ui")  # For output elements
                   )
                 )
        )
      )
    )
  )
}

#' Add external Resources to the Application
#'
#' This function is internally used to add external
#' resources inside the Shiny application.
#'
#' @import shiny
#' @importFrom golem add_resource_path activate_js favicon bundle_resources
#' @noRd
golem_add_external_resources <- function() {
  add_resource_path(
    "www",
    app_sys("app/www")
  )

  tags$head(
    favicon(),
    bundle_resources(
      path = app_sys("app/www"),
      app_title = "baerhunter"
    )
    # Add here other external resources
    # for example, you can add shinyalert::useShinyalert()
  )
}
