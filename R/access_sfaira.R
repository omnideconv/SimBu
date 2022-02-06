#' setup the sfaira package
#'
#' If you want to download datasets from sfaira, you need to setup the connection to their database first.
#' Sfaira should be installed in its own conda environment; the python-path used inside of this environment
#' as well as the name of the environment have to be given to this function. Additionally you have to give the
#' path to a directory where the sfaira datasets will be downloaded into.
#'
#'
#' @param python_path path to the python executable
#' @param env_name name of the conda environment in which sfaira is installed
#' @param basedir name of the directory, where the raw files will be downloaded into
#'
#' @return list with sfaira package and file directories; must be used for other sfaira functions
#' @export
#'
#' @examples
#' \dontrun{
#' sfaira_list <- setup_sfaira(python_path="/path/to/conda/envs/sfaira/bin/python3",env_name="sfaira", basedir="/home/sfaira_data")
#' }
setup_sfaira <- function(python_path, env_name, basedir){
  tryCatch({
    # check if sfaira is installed in environment
    reticulate::use_python(python_path, required = T)
    reticulate::use_condaenv(env_name, required = T)
    sfaira <- reticulate::import("sfaira")

    # print version of currently used sfaira version
    check_cmd <- "sfaira --version"
    system(check_cmd)

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
#' download a specific dataset from sfaira by an ID
#'
#' @param setup_list the sfaira setup; given by \code{\link{setup_sfaira}}
#' @param id the ID of the datasets
#' @param force logical; TRUE if you want to force the download, even though no cell-type annotation exists for this dataset. Default if FALSE
#' @param synapse_user character; username for synapse portal (https://www.synapse.org)
#' @param synapse_pw character; password for synapse portal (https://www.synapse.org)
#'
#' @return a anndata object, stores counts, metadata and other information on the dataset
#'
download_sfaira <- function(setup_list, id, force=F, synapse_user=NULL, synapse_pw=NULL){
  sfaira <- setup_list[["sfaira"]]
  ds <- sfaira$data$Universe(data_path = setup_list[["rawdir"]],
                             meta_path = setup_list[["metadir"]],
                             cache_path = setup_list[["cachedir"]])

  # download and load count matrix of given sfaira ID
  tryCatch({
    ds$subset(key="id", values=c(id))
    ds$download()
    ds$load()
    is_annotated <- ds$datasets[[id]]$annotated
    if(!is_annotated && force){
      warning("The downloaded dataset has no annotation; this might get you into issues down the road.")
    }else if(!is_annotated && !force){
      warning("The downloaded dataset has no annotation; you can force the download regardless with force=TRUE.")
      return(NULL)
    }
    #streamline features & meta-data
    ds$datasets[[id]]$streamline_features(match_to_release = pkg.globals$sfaira_streamline_feature_version)
    ds$datasets[[id]]$streamline_metadata(schema="sfaira")
    return(ds$datasets[[id]]$adata)
  }, error = function(e){
    message(paste0("Could not download dataset for id: ", id))
    print(e$message)
    return(NULL)
  })

}

#' download multiple datasets from sfaira using filters for organism, tissue and/or assay
#'
#' similar to the filters on the sfaira website (\url{https://theislab.github.io/sfaira-portal/Datasets})
#'
#' @param setup_list the sfaira setup; given by \code{\link{setup_sfaira}}
#' @param organisms list of organisms (only human and mouse available)
#' @param tissues list of tissues
#' @param assays list of assays
#' @param force logical; TRUE if you want to force to download all datasets, otherwise only the ones with cell-type annotation will be returned. Default if FALSE
#'
#' @return
#'
download_sfaira_multiple <- function(setup_list, organisms=NULL, tissues=NULL, assays=NULL, force=F){
  sfaira <- setup_list[["sfaira"]]

  tryCatch({

    if(force){
      warning("Some or all of the downloaded datasets have no annotation; this might get you into issues down the road.")
      ds <- sfaira$data$Universe(data_path = setup_list[["rawdir"]],
                                 meta_path = setup_list[["metadir"]],
                                 cache_path = setup_list[["cachedir"]])
    }else{
      # python code to subset sfaira universe by annotated datasets
      print("Removing datasets without cell-type annotation...")
      reticulate::py_run_string("import sfaira")
      reticulate::py_run_string(paste0("ds = sfaira.data.Universe(data_path=\'",setup_list[["rawdir"]],
                           "\', meta_path=\'",setup_list[["metadir"]],
                           "\', cache_path=\'",setup_list[["cachedir"]],"\')"))
      reticulate::py_run_string("dsg = sfaira.data.DatasetGroup(datasets=dict([(k, v) for k, v in ds.datasets.items() if v.annotated]), collection_id='something')")
      ds <- reticulate::py$dsg
    }

    # apply filters on sfaira database
    if(all(is.null(c(organisms, tissues, assays)))) stop("You must specify at least one filter.", call.=F)
    if(!is.null(organisms)) {ds$subset(key="organism", values=organisms)}
    if(!is.null(assays)) {ds$subset(key="assay_sc", values=assays)}
    if(!is.null(tissues)) {ds$subset(key="organ", values=tissues)}
    if(length(ds$datasets) == 0){stop("No datasets found with these filters; please check again", call.=F)}

    print("Downloading datasets...")
    ds$download()
    ds$load()
    #streamline features and meta-data
    print("Streamlining features & meta-data...")
    ds$streamline_features(match_to_release = pkg.globals$sfaira_streamline_feature_version)
    ds$streamline_metadata(schema = "sfaira")
    return(ds$adata)

  }, error=function(e){
    message(paste0("Could not download all datasets for specified filters."))
    print(e$message)
    return(NULL)
  })
}


#' Gives an overview of the possible datasets you can use from the sfaira database
#'
#' @param setup_list the sfaira setup; given by \code{\link{setup_sfaira}}
#'
#' @return a dataframe with information on each dataset
#' @export
#'
#' @examples
#' \dontrun{all_datasets <- sfaira_overview(setup_list)}
sfaira_overview <- function(setup_list){
  sfaira <- setup_list[["sfaira"]]
  ds <- sfaira$data$Universe(data_path = setup_list[["rawdir"]],
                             meta_path = setup_list[["metadir"]],
                             cache_path = setup_list[["cachedir"]])
  all_datasets <- ds$datasets
  info_list <- lapply(all_datasets, function(x){
    return(list(id=x$id,
                author=x$author,
                doi=x$doi,
                annotated=x$annotated,
                normalization=x$normalization,
                assay=x$assay_sc,
                organ=x$organ,
                organism=x$organism))
  })

  suppressWarnings(out <- data.table::rbindlist(info_list))
  out$annotated[which(is.na(out$annotated))]<-F
  return(out)
}
