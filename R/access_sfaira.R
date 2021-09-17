#' setup the sfaira package
#'
#' @param python_path path to the python executable
#' @param env_name name of the conda environment in which sfaira is installed
#' @param basedir name of the directory, where the raw files will be downloaded into
#'
#' @return list with sfaira package and file directories; must be used for other sfaira functions
#' @export
#'
#' @examples
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
    #streamline features & meta-data
    ds$adata
    ds$datasets[[id]]$streamline_features(match_to_reference=list("human"= "Homo_sapiens.GRCh38.102", "mouse"= "Mus_musculus.GRCm38.102"))
    ds$datasets[[id]]$streamline_metadata(schema="sfaira")
    return(adata)
  }, error = function(e){
    message(paste0("Could not download dataset for id: ", id))
    print(e$message)
    return(NULL)
  })

}

download_sfaira_multiple <- function(setup_list, organisms=NULL, tissues=NULL, assays=NULL, force=F){
  sfaira <- setup_list[["sfaira"]]
  ds <- sfaira$data$Universe(data_path = setup_list[["rawdir"]],
                             meta_path = setup_list[["metadir"]],
                             cache_path = setup_list[["cachedir"]])

  tryCatch({

    # apply filters on sfaira database
    if(all(is.null(c(organism, tissue, assay)))) stop("You must specify at least one filter.", call.=F)
    if(!is.null(organism)) {ds$subset(key="organism", value=organisms)}
    if(!is.null(assays)) {ds$subset(key="assay_sc", value=assays)}
    if(!is.null(tissues)) {ds$subset(key="organ", value=tissues)}
    if(length(ds$datasets) == 0){stop("No datasets found with these filters; please check again", call.=F)}

    # only proceed with datasets with cell-type annotation
    if(!force){
      has_annotation <- lapply(ds$datasets, function(i){i$annotated})
      annotated_sets <- names(has_annotation[which(as.logical(has_annotation))])
      #TODO subset by annotation
    }else{
      warning("Some or all of the downloaded datasets have no annotation; this might get you into issues down the road.")
    }

    #streamline features and meta-data
    ds$streamline_features(match_to_reference = )

  }, error=function(e){
    message(paste0("Could not download all datasets for specified filters."))
    print(e$message)
    return(NULL)
  })
}
