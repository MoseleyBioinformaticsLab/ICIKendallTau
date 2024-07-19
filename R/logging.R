#' log memory usage
#'
#' Logs the amount of memory being used to a log file if it is available, and
#' generating warnings if the amount of RAM hits zero.
#'
#' @export
#' @return NULL
log_memory = function(){
  if (get("memory", envir = icikt_logger)) {
    linux_memory = system("cat /proc/meminfo", intern = TRUE)
    linux_memory = grep("^MemTotal|^MemAvailable|^Active|^SwapTotal|^SwapFree", linux_memory, value = TRUE)
    linux_memory = grep("anon|active|file", linux_memory, value = TRUE, invert = TRUE)
    
    memory_values = stringr::str_extract(linux_memory, "[[:digit:]].*")
    memory_numbers = as.numeric(stringr::str_extract(memory_values, "[[:digit:]].* "))
    memory_ids = stringr::str_extract(linux_memory, "^[[:alpha:]]+")
    names(memory_numbers) = memory_ids
    memory_string = paste0("Memory: ", paste(paste(c("Total: ", "Available: ", "Active: ", "SwapTotal: ", "SwapFree: "), memory_values, sep = ""), collapse = ", "))
    
    active_to_total = memory_numbers["Active"] / memory_numbers["MemTotal"]
    swapfree_to_swap = memory_numbers["SwapFree"] / memory_numbers["SwapTotal"]
    if (is.nan(swapfree_to_swap)) {
      swapfree_to_swap = 0
    }
    
    swapfree_to_swap = (memory_numbers["SwapTotal"] - memory_numbers["SwapFree"]) / memory_numbers["SwapTotal"]
    if (is.nan(swapfree_to_swap)) {
      swapfree_to_swap = 0
    }
    
    
    if ((active_to_total >= 0.95) || (swapfree_to_swap >= 0.7)) {
      memory_string2 = paste0("HIGH MEMORY USAGE!!! ", memory_string)
      if (get("logger", envir = icikt_logger)) {
        logger::log_warn(memory_string2, namespace = "ICIKendallTau")
      } else {
        warning(memory_string2)
      }
    } else {
      if (get("logger", envir = icikt_logger)) {
        logger::log_info(memory_string, namespace = "ICIKendallTau")
      }
    }
  }
  
}

#' log messages
#'
#' If a log_appender is available, logs the given message at the `info` level.
#'
#' @param message_string the string to put in the message
#'
#' @export
#' @return NULL
log_message = function(message_string){
  if (get("logger", envir = icikt_logger)) {
    logger::log_info(message_string, namespace = "ICIKendallTau")
  }
  if (get("status", envir = icikt_progress)) {
    message(message_string)
  }
}

#' turn logging off
#'
#' There may be good reasons to turn the logging off after it's been turned on. This
#' basically tells the package that the logger isn't available.
#'
#' @export
#' @return NULL
disable_logging = function(){
  assign("logger", FALSE, envir = icikt_logger)
  message("Logging disabled.")
}


#' turn logging on
#'
#' Choose to enable logging, to a specific file if desired.
#'
#' @param log_file the file to log to
#' @param memory provide memory logging too? Only available on Linux and MacOS
#'
#' @details Uses the logger package under the hood, which is suggested in the dependencies.
#'   Having logging enabled is nice to see when things are starting and stopping, and what exactly
#'   has been done, without needing to write messages to the console. It is especially
#'   useful if you are getting errors, but can't really see them, then you can add
#'   "memory" logging to see if you are running out of memory.
#'
#'   Default log file has the pattern:
#'
#'   YYYY.MM.DD.HH.MM.SS_ICIKendallTau_run.log
#'
#' @export
#' @return NULL
enable_logging = function(log_file = NULL, memory = FALSE){
  has_logger = requireNamespace("logger", quietly = TRUE)
  if (!has_logger) {
    stop("logger package is not available. Please install it to enable logging!\ninstall.packages('logger')")
  } else {
    assign("logger", TRUE, envir = icikt_logger)
    # if no log file supplied, and we see an old one, just use it
    if (!is.null(get("log_file", envir = icikt_logger)) && is.null(log_file)) {
      log_file = get("log_file", envir = icikt_logger)
    } else if (is.null(log_file)) {
      log_file = paste0(substring(make.names(Sys.time()), 2), "_ICIKendallTau_run", ".log")
    }
    
    assign("log_file", log_file, envir = icikt_logger)
    if (memory) {
      sys_info = Sys.info()
      if (!grepl("windows", sys_info["sysname"], ignore.case = TRUE)) {
        assign("memory", TRUE, envir = icikt_logger)
      } else {
        message("Memory logging is not available on Windows!\nMemory use will not be logged.")
      }
    }
    logger::log_appender(logger::appender_file(log_file), namespace = "ICIKendallTau")
  }
  
}

#' turn progress on off
#'
#' Allow the user to turn progress messages to the console and off. Default
#' is to provide messages to the console.
#'
#' @param progress logical to have it on or off
#'
#' @export
#' @return NULL
show_progress = function(progress = TRUE){
  assign("status", progress, envir = icikt_progress)
}

icikt_logger = new.env(hash = TRUE)
assign("logger", FALSE, envir = icikt_logger)
assign("memory", FALSE, envir = icikt_logger)
assign("log_file", NULL, envir = icikt_logger)

icikt_progress = new.env(hash = TRUE)
assign("status", FALSE, envir = icikt_progress)
