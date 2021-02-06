
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
    if(verbose){
      base::suppressMessages(method)
    } else {
      method
    }
  })
}


#' Docker availability check
#'
#' @return A boolean value whether docker is available on the system
#' @export
#'
#' @examples
docker_available <- function(){
  return(system("docker",ignore.stdout = TRUE,ignore.stderr = TRUE)==0)
}
