check_if_colnames_null <- function(data_matrix, 
                                   call = rlang::caller_env()){
  no_colnames <- is.null(colnames(data_matrix))
  if (no_colnames) {
    cli::cli_abort(
      message = c("Colnames of {.var data_matrix} must be be specified",
                  "Currently {.code colnames(data_matrix)} returns \\
                  {colnames(data_matrix)}"),
      call = call
    )
  }
}

