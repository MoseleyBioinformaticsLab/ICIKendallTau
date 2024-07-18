setup_missing_matrix = function(data_matrix, global_na)
{
  exclude_loc = matrix(FALSE, nrow = nrow(data_matrix), ncol = ncol(data_matrix))
  rownames(exclude_loc) = rownames(data_matrix)
  colnames(exclude_loc) = colnames(data_matrix)
  if (length(global_na) > 0) {
    if (any(is.na(global_na))) {
      exclude_loc[is.na(data_matrix)] = TRUE
      global_na = global_na[!is.na(global_na)]
    }
    if (any(is.infinite(global_na))) {
      exclude_loc[is.infinite(data_matrix)] = TRUE
      global_na = global_na[!is.infinite(global_na)]
    }
  }
  if (length(global_na) > 0) {
    for (ival in global_na) {
      exclude_loc[data_matrix == ival] = TRUE
    }
  }
  
  return(exclude_loc)
}

check_if_colnames_null = function(data_matrix,
                                  arg = rlang::caller_arg(data_matrix),
                                  call = rlang::caller_env())
{
  if (is.null(colnames(data_matrix))) {
    cli::cli_abort(
      message = c('Colnames of {.arg {arg}} must be be specified.',
                  'x' = 'Currently {.code colnames({arg})} returns \\
                  {colnames(data_matrix)}'),
    call = call
    )  
  }
  
}

check_if_numeric = function(data_matrix,
                            arg = rlang::caller_arg(data_matrix),
                            call = rlang::caller_env())
{
  if (!(is.double(data_matrix) || is.integer(data_matrix))) {
    cli::cli_abort(
      message = c('{.arg {arg}} must be a numeric type.',
                  'x' = 'Currently, {.code class({arg}[1])} returns \\
                  {.obj_type_friendly {data_matrix[1]}}'),
      call = call
    )
  }  
}

transform_to_matrix = function(data_matrix,
                                arg = rlang::caller_arg(data_matrix),
                                call = rlang::caller_env())
{
  if (is.data.frame(data_matrix)) {
    cli::cli_inform(message = c(
      'i' = '{.arg {arg}} is a data.frame, converting to matrix ...'
    ),
    call = call)
    return(as.matrix(data_matrix))
  }
  return(data_matrix)
}

check_furrr = function()
{
  if (requireNamespace("furrr", quietly = TRUE)) {
    ncore = future::nbrOfWorkers()
    names(ncore) = NULL
    split_fun = furrr::future_map
  } else {
    ncore = 1
    split_fun = purrr::map
  }
  return(list(ncore = ncore,
              split_fun = split_fun))
}

check_na_method = function(na_method,
                          na_arg = rlang::caller_arg(na_method),
                          na_call = rlang::caller_env())
{
  if (na_method %in% 'na.or.complete') {
    cli::cli_abort(message = c(
      '\'na.or.complete\' is not a supported value for {.arg use}.',
      'i' = 'Please use one of {.emph all.obs complete.obs pairwise.complete everthing}.'
    ),
  call = na_call)
  }
}

check_x_y = function(x, y,
                     x_arg = rlang::caller_arg(x),
                     y_arg = rlang::caller_arg(y),
                    x_call = rlang::caller_env())
{
  dim_x = dim(x)

  if (is.null(dim_x) && is.null(y)) {
    cli::cli_abort(message = c(
      '{.arg {x_arg}} and {.arg y} should both be provided as vectors, or {.arg {x_arg}} should be matrix-like.',
      'i' = '{.arg {x_arg}} is a single vector, and {.arg y} is `NULL`.'
    ),
    call = x_call)
  }
  
  
  if (!is.null(dim_x)) {
    check_if_colnames_null(x, x_arg, x_call)
  }
  x = transform_to_matrix(x, x_arg, x_call)
  check_if_numeric(x, x_arg, x_call)

  if (!is.null(y)) {
    dim_y = dim(y)

    if (!is.null(dim_x) || !is.null(dim_y)) {
      cli::cli_abort(message = c(
        'Both {.arg {x_arg}} and {.arg {y_arg}} must be vectors.',
        'i' = '{.code ncol({x_arg})} returns {ncol(x)}, {.code ncol({y_arg})} returns {ncol(y)}.'
      ))
    } 
    
    y = transform_to_matrix(y, y_arg, x_call)
    check_if_numeric(y, y_arg, x_call)

    x = cbind(x, y)
    colnames(x) = make.names(c(x_arg, y_arg))
    
  }
  x
}