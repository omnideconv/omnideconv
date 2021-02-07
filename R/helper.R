
#' A wrapper function whether to suppress messages
#'
#' @param verbose Whether to suppress messages
#'
#' @return A function which will suppress messages or not, depending on the verbose parameter
#' @export
#'
verbose_wrapper <- function(verbose){
  return (function(method){
    ifelse(verbose, base::suppressMessages(method), base::identity(method))
  })
}


#' Docker availability check
#'
#' @return A boolean value whether docker is available on the system
#' @export
#'
docker_available <- function(){
  return(system("docker",ignore.stdout = TRUE,ignore.stderr = TRUE)==0)
}

#' Docker connectability check
#'
#' @return A boolean value whether it is possible to connect to docker
#' @export
#'
docker_connectable <- function(){
  return(system("docker ps",ignore.stdout = TRUE,ignore.stderr = TRUE)==0)
}
