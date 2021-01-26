call_docker <- function(command){
  return(system(paste("docker ",command)))
}


docker_available <- function(){
  return(system("docker",ignore.stdout = T,ignore.stderr = T)==0)
}


if (docker_available()){
  print("DoStuff")
} else {
  print("Installation of docker can not be found. Please check whether you can
        call 'docker' in the command line and get a help menu")
}

#docker run -v {dir_path}:/src/data -v {dir_path}:/src/outdir cibersortx/fractions
#--refsample Fig2ab-NSCLC_PBMCs_scRNAseq_refsample.txt --mixture Fig2b-WholeBlood_RNAseq.txt --fraction 0 --rmbatchSmode TRUE
create_docker_command <- function(in_dir, out_dir, method = c("create_sig","impute_cell_fractions"),verbose = FALSE, ...){
  base <- paste0("docker run -v ",in_dir,":/src/data -v ",out_dir,":/src/outdir cibersortx/fractions --single_cell TRUE")
  if (verbose){
    base <- paste(base,"--verbose TRUE")
  }
  credentials <- "--username konstantin.pelz@tum.de --token token_obtained_from_CIBERSORTx_website"
  return(paste(base,credentials,get_method_options(method,...)))
}

get_method_options <- function(method,...){
  if (method=="create_sig"){
    return(get_signature_matrix_options(...))
  } else if (method == "impute_cell_fractions"){
    return(get_signature_matrix_options(...))
  } else {
    base::stop(paste("Method",method,"is not valid"))
  }
}

get_signature_matrix_options <- function(refsample, G_min = 300, G_max = 500, q_value = 0.01, filter = FALSE,
                                         k_max = 999, remake = FALSE, replicates = 5, sampling = 0.5, fraction = 0.75){
  return(paste("--refsample",refsample,"--G.min",G_min,"--G.max",G_max,"--q.value",q_value,"--filter",filter,
               "--k.max",k_max,"--remake",remake,"--replicates",replicates,"--sampling",sampling,"--fraction",fraction))
}


get_cell_fractions_options <- function(mixture, sigmatrix){

}

#Working Commands
#docker run -v /c/Users/Konstantin/Desktop/Uni/7Semester/SysBioMed/testIn:/src/data -v /c/Users/Konstantin/Desktop/Uni/7Semester/SysBioMed/testOut:/src/outdir cibersortx/fractions --username mailto:gregor.sturm@cs.tum.edu --token 1860fb7bc414a958e8aa1e91a5229d8a --single_cell TRUE --refsample sampleFileForCibersort.txt





