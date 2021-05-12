#' A wrapper function whether to suppress messages
#'
#' @param verbose Whether to suppress messages
#'
#' @return A function which will suppress messages or not, depending on the verbose parameter
#'
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
#'
docker_available <- function(){
  return(system("docker",ignore.stdout = TRUE,ignore.stderr = TRUE)==0)
}

#' Docker connectability check
#'
#' @return A boolean value whether it is possible to connect to docker
#'
docker_connectable <- function(){
  return(system("docker ps",ignore.stdout = TRUE,ignore.stderr = TRUE)==0)
}

#' Removes blanks by substituting them with "$_$" which should not be used naturally
#'
#' @param string The string to be escaped
#'
#' @return The String without blanks
#'
escape_blanks <- function(string){
  return(gsub(" ", "$_$", string))
}

#' Removes the substitutions "$_$" and turns them back into blanks
#'
#' @param string The string to be de-escaped
#'
#' @return The String with blanks
#'
deescape_blanks <- function(string){
  return(gsub("\\$_\\$", " ", string))
}
