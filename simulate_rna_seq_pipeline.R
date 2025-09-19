# Updated log_message function

# This function logs messages with proper formatting to fix the sprintf error.
log_message <- function(message) {
    formatted_message <- sprintf("%s - %s", Sys.time(), message)
    cat(formatted_message, "\n")
}