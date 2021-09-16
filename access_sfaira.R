setup_sfaira <- function(python_path, env_name, basedir){
  tryCatch({
    # check if sfaira is installed in environment
    reticulate::use_python(python_path, required = T)
    reticulate::use_condaenv(env_name, required = T)
    sfaira <- reticulate::import("sfaira")
    
    # create directories for sfaira downloads
    if(!dir.exists(basedir)) dir.create(basedir)
    cachedir <- paste0(basedir, "/cache") 
    if(!dir.exists(cachedir)) dir.create(cachedir)
    rawdir <- paste0(basedir, "/raw") 
    if(!dir.exists(rawdir)) dir.create(rawdir)
    metadir <- paste0(basedir, "/meta") 
    if(!dir.exists(metadir)) dir.create(metadir)
    
    return(list(sfaira=sfaira,
                cachedir=cachedir,
                rawdir=rawdir,
                metadir=metadir))
  }, error=function(e){
    print(e$message)
  })
}

# download from sfaira
# force = force function to return, even though no annotation exists
download_sfaira <- function(setup_list, id, force=F){
  sfaira <- setup_list[["sfaira"]]
  ds <- sfaira$data$Universe(data_path = setup_list[["rawdir"]],
                             meta_path = setup_list[["metadir"]],
                             cache_path = setup_list[["cachedir"]])
  
  # download and load count matrix of given sfaira ID
  tryCatch({
    ds$datasets[[id]]$download()
    ds$datasets[[id]]$load() 
    is_annotated <- ds$datasets[[id]]$annotated
    if(!is_annotated && force){
      warning("The downloaded dataset has no annotation; this might get you into issues down the road.")
    }else if(!is_annotated && !force){
      warning("The downloaded dataset has no annotation; you can force the download regardless with force=TRUE.")
      return(NULL)
    }
    adata <- ds$adata
    return(adata)
  }, error = function(e){
    print(paste0("Could not download dataset for id: ", id))
    print(e$message)
    return(NULL)
  })

}
