#' serve_data2 (I copied these functions from the source code of Jbrowsr from github; I just modified the
#' generate_cli_message code so that the urls can be returned and used in Rshiny)
#'
#' @description A fct function
#'
#' @return The return value, if any, from executing the function.
#'
#' @noRd

strip_gz <- function(track_data) {
  if (stringr::str_ends(track_data, ".gz")) {
    stringr::str_trunc(track_data, nchar(track_data) - 3, "right", "")
  } else {
    track_data
  }
}



serve_data2 <- function(path, port = 5000) {
  path <- normalizePath(path, mustWork = TRUE)

  if (path != ".") {
    # clean up R directory after exit if not current directory
    original_dir <- setwd(path)
    on.exit(setwd(original_dir))
  }

  server <- create_server(path, port)
  app <- list(call = serve_directory(path))
  server$start_server(app)
  invisible(server)
}

create_server <- function(path, port) {
  host <- "127.0.0.1"
  server <- NULL
  list(
    host = host,
    port = as.integer(port),
    interval = 1,
    url = sprintf("http://%s:%d", "127.0.0.1", as.integer(port)),
    start_server = function(app) {
      server_id <- httpuv::startServer(host, port, app)
      server <<- server_id
      invisible(server_id)
    },
    stop_server = function() {
      if (is.null(server)) {
        stop("The server has not been started.")
      } else {
        httpuv::stopServer(server)
      }
    }
  )
}

serve_directory <- function(path) {
  function(req) {
    original_dir <- setwd(path)
    on.exit(setwd(original_dir))
    req_path <- resolve_req_path(req)

    # response content
    if (is.null(req$HTTP_RANGE)) {
      status <- 200L
    } else {
      status <- 206L
      range <- parse_range(req$HTTP_RANGE, req_path)
    }
    type <- guess_type(path)
    body <- generate_body(req, req_path)
    if (is.character(body) && length(body) > 1) {
      body <- stringr::str_c(body, collapse = "")
    }
    headers <- c(
      list(
        "Access-Control-Allow-Origin" = "*",
        "Access-Control-Allow-Headers" = "range, Origin, Content-Type",
        "Access-Control-Expose-Headers" = "Content-Length,Content-Range",
        "Content-Type" = type
      ), if (status == 206L) {
        list(
          "Content-Range" = paste0(
            "bytes ",
            range[2],
            "-",
            range[3],
            "/",
            file.info(req_path)[, "size"]
          )
        )
      },
      "Accept-Ranges" = "bytes"
    )

    list(
      status = status,
      body = body,
      headers = headers
    )
  }
}

resolve_req_path <- function(req) {
  req_path <- httpuv::decodeURIComponent(req$PATH_INFO)
  Encoding(req_path) <- "UTF-8"

  if (grepl("^/", req_path)) {
    stringr::str_c(".", req_path)
  } else if (req_path == "") {
    "."
  }
}

generate_body <- function(req, req_path) {
  range <- req$HTTP_RANGE

  if (is.null(range)) {
    read_raw(req_path)
  } else {
    range <- parse_range(range, req_path)
    byte2 <- as.numeric(range[2])
    byte3 <- as.numeric(range[3])

    if (length(range) < 3 || (range[1] != "bytes") || (byte2 >= byte3)) {
      return(
        list(
          status = 416L, headers = list("Content-Type" = "text/plain"),
          body = "Requested range not satisfiable\r\n"
        )
      )
    }

    connection <- file(req_path, open = "rb", raw = TRUE)
    on.exit(close(connection))
    seek(connection, where = byte2, origin = "start")
    readBin(connection, "raw", byte3 - byte2 + 1)
  }
}

read_raw <- function(req_path) {
  readBin(
    req_path,
    "raw",
    file.info(req_path)[, "size"]
  )
}

guess_type <- function(path) {
  mimetype <- function(...) {
    system2("mimetype", c("-b", shQuote(path)), ...)
  }
  if (Sys.which("mimetype") == "" || mimetype(stdout = NULL) != 0) {
    return(mime::guess_type(path))
  }
  mimetype(stdout = TRUE)
}

parse_range <- function(range, req_path) {
  range <- strsplit(range, split = "(=|-)")[[1]]
  if (length(range) == 2 && range[1] == "bytes") {
    # open-end range request
    range[3] <- file.info(req_path)[, "size"] - 1
  }
  range
}

generate_cli_message <- function(path, port) {
  data_files <- dir(path)

  cli::cli_par()
  cli::cli_text("Serving data from: {path}")
  cli::cli_end()

  cli::cli_par()
  cli::cli_end()

  # Initialize vectors to store the URLs

  assembly <- c()
  alignments <- c()
  feature <- c()
  variant <- c()
  wiggle <- c()

  for (i in seq_along(data_files)) {
    stripped_gz <- strip_gz(data_files[i])

    if (stringr::str_ends(stripped_gz, ".fa") || stringr::str_ends(stripped_gz, ".fasta")) {
      assembly <- c(assembly, stringr::str_glue("http://127.0.0.1:{port}/{data_files[i]}"))
    } else if (stringr::str_ends(stripped_gz, ".bam") || stringr::str_ends(stripped_gz, ".cram")) {
      alignments <- c(alignments, stringr::str_glue("http://127.0.0.1:{port}/{data_files[i]}"))
    } else if (stringr::str_ends(stripped_gz, ".gff") || stringr::str_ends(stripped_gz, ".gff3")) {
      feature <- c(feature, stringr::str_glue("http://127.0.0.1:{port}/{data_files[i]}"))
    } else if (stringr::str_ends(stripped_gz, ".vcf")) {
      variant <- c(variant, stringr::str_glue("http://127.0.0.1:{port}/{data_files[i]}"))
    } else if (stringr::str_ends(stripped_gz, ".bw") || stringr::str_ends(stripped_gz, ".bigWig")) {
      wiggle <- c(wiggle, stringr::str_glue("http://127.0.0.1:{port}/{data_files[i]}"))
    }
  }

  # Compile all URLs into a list
  urls <- list(
    assembly = assembly,
    alignments = alignments,
    feature = feature,
    variant = variant,
    wiggle = wiggle
  )

  cli::cli_par()
  #cli::cli_alert_info("Note: this server is intended for local development and use.")
  #cli::cli_text(
  "For a production deployment, refer to the vignette about creating URLs for a discussion of better options."
  #)
  cli::cli_end()

  #cli::cli_text("To stop the server, use the $stop_server() function.")

  return(urls)
}

log_track_message <- function(title, track_vector, port) {
  if (length(track_vector) > 0) {
    cli::cli_par()
    cli::cli_h2(title)
    cli::cli_ul(track_vector)
    cli::cli_end()
  }

}

