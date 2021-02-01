
#' A wrapper function whether to suppress messages
#'
#' @param verbose Whether to suppress messages
#'
#' @return A function which will suppress messages or not, depending on the verbose parameter
#' @export
#'
#' @examples
verbose_wrapper <- function(verbose){
  return (function(method){
    ifelse(verbose, base::suppressMessages(method), base::identity(method))
  })
}
