#' setup the sfaira package
#'
#' If you want to download datasets from Sfaira, you need to specify a directory where the datasets are saved into.
#' Additionally, when this function is called for the first time, a conda environment will be established and sfaira along
#' all of its dependencies are installed. This can take some time but will be only performed one single time, as the 
#' environment can be re-used.
#'
#'
#' @param basedir name of the directory, where the raw files will be downloaded into
#'
#' @return list with sfaira file directories; must be used as input for other sfaira based functions
#' @export
#'
#' @examples
#' setup_list <- setup_sfaira(basedir=tempdir())
setup_sfaira <- function(basedir){
  tryCatch({
    # create conda environment with sfaira
    proc <- basilisk::basiliskStart(SimBu_env)
    on.exit(basilisk::basiliskStop(proc))
    
    # initialize sfaira environment
    basilisk::basiliskRun(proc, function(){
      sfaira <- reticulate::import("sfaira")
      print('Sucessfully loaded sfaira.')
    })

    # create directories for sfaira downloads
    if(!dir.exists(basedir)) dir.create(basedir)
    cachedir <- paste0(basedir, "/cache")
    if(!dir.exists(cachedir)) dir.create(cachedir)
    rawdir <- paste0(basedir, "/raw")
    if(!dir.exists(rawdir)) dir.create(rawdir)
    metadir <- paste0(basedir, "/meta")
    if(!dir.exists(metadir)) dir.create(metadir)

    return(list(cachedir=cachedir,
                rawdir=rawdir,
                metadir=metadir))
  }, error=function(e){
    print(e$message)
  })
}

# download from sfaira
# force: force the function to return the dataset, even though no annotation exists; this will result in an warning message.
#' download a specific dataset from sfaira by an ID
#'
#' @param setup_list the sfaira setup; given by \code{\link{setup_sfaira}}
#' @param ids the IDs of the datasets
#' @param force logical; TRUE if you want to force the download, even though no cell-type annotation exists for this dataset. Default if FALSE
#' @param synapse_user character; username for synapse portal (https://www.synapse.org)
#' @param synapse_pw character; password for synapse portal (https://www.synapse.org)
#'
#' @return matrix, gene names and cell IDs
#'
download_sfaira <- function(setup_list, ids, force=FALSE, synapse_user=NULL, synapse_pw=NULL){
  # create conda environment with sfaira
  proc <- basilisk::basiliskStart(SimBu_env)
  on.exit(basilisk::basiliskStop(proc))
  
  # download and load count matrix of given sfaira ID
  sfaira_dataset <- basilisk::basiliskRun(proc, function(setup_list, ids, force, feature_version){
    sfaira <- reticulate::import("sfaira")
    
    ds <- sfaira$data$Universe(data_path = setup_list[["rawdir"]],
                               meta_path = setup_list[["metadir"]],
                               cache_path = setup_list[["cachedir"]])
    
    tryCatch({
      
      #check if id(s) is(are) present in sfaira
      if(all(!ids %in% names(ds$datasets))){
        stop('None of the supplied IDs is present in Sfaira. Stopping.', call. = FALSE)
      }
      else if(any(!ids %in% names(ds$datasets))){
        warning('One or more of the IDs you supplied are not availabe in Sfaira and will not be downloaded. Please check for consistent spelling.')
        warning('The missing IDs are:')
        warning(paste0(ids, collapse = '; '))
      }
      ds$subset(key="id", values=c(ids))
      
      #check if all sfaira datasets have accessible files
      accessibility <- check_sfaira_datasets(ds)
      accessible_ids <- names(accessibility$ids_with_data)
      
      accessible_ids <- skip_broken_datasets(accessible_ids)
      
      if(length(accessible_ids) == 0){
        stop('Cannot download data and/or metadata for any of the supplied sfaira IDs. Stopping.',call. = FALSE)
      }
      
      if(length(accessible_ids) != length(ds$ids)){
        missing_ids <- ids[which(!ids %in% accessible_ids)]
        warning('Cannot download data and/or metadata for the following sfaira IDs:')
        warning(paste0(missing_ids, collapse = '; '))
      }
      
      ds$subset(key="id", values=c(accessible_ids))
      
      print("Downloading datasets...")
      ds$download()
      ds$load()
      
      is_annotated <- ds$datasets[[ids]]$annotated
      if(!is_annotated && force){
        warning("The downloaded dataset has no annotation; this might get you into issues down the road.")
      }else if(!is_annotated && !force){
        warning("The downloaded dataset has no annotation; you can force the download regardless with force=TRUE.")
        return(NULL)
      }
      #streamline features & meta-data
      print("Streamlining features & meta-data...")
      ds$datasets[[ids]]$streamline_features(match_to_release = feature_version)
      ds$datasets[[ids]]$streamline_metadata(schema="sfaira")
      
      adata <- ds$datasets[[ids]]$adata
      X <- methods::as(methods::as(adata$X, 'CsparseMatrix'), 'dgCMatrix')
      obs <- adata$obs
      var <- adata$var
      return(list(X=X, obs=obs, var=var))
    }, error = function(e){
      message(paste0("Could not download dataset for id."))
      print(e$message)
      return(NULL)
    })
    
  }, setup_list=setup_list, ids=ids, force=force, feature_version = pkg.globals$sfaira_streamline_feature_version)
  
  return(sfaira_dataset)

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
#' @return annotated data object, contains count matrix and annotation
#'
download_sfaira_multiple <- function(setup_list, organisms=NULL, tissues=NULL, assays=NULL, force=FALSE){
  # create conda environment with sfaira
  proc <- basilisk::basiliskStart(SimBu_env)
  on.exit(basilisk::basiliskStop(proc))
  
  
  sfaira_dataset <- basilisk::basiliskRun(proc, function(setup_list, organisms, tissues, assays, force, feature_version){
    
    sfaira <- reticulate::import("sfaira")
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
      if(all(is.null(c(organisms, tissues, assays)))) stop("You must specify at least one filter.", call.=FALSE)
      if(!is.null(organisms)) {ds$subset(key="organism", values=organisms)}
      if(!is.null(assays)) {ds$subset(key="assay_sc", values=assays)}
      if(!is.null(tissues)) {ds$subset(key="organ", values=tissues)}
      if(length(ds$datasets) == 0){stop("No datasets found with these filters; please check again", call.=FALSE)}
      
      #check if all sfaira datasets have accessible files
      accessibility <- check_sfaira_datasets(ds)
      accessible_ids <- names(accessibility$ids_with_data)
      
      accessible_ids <- skip_broken_datasets(accessible_ids)
      
      if(length(accessible_ids) == 0){
        stop('Cannot download data and/or metadata for any of the supplied sfaira IDs. Stopping.',call. = FALSE)
      }
      
      if(length(accessible_ids) != length(ds$ids)){
        missing_ids <- ds$ids[which(!ds$ids %in% accessible_ids)]
        warning('Cannot download data and/or metadata for the following sfaira IDs:')
        warning(paste0(missing_ids, collapse = '; '))
      }
      
      ds$subset(key="id", values=c(accessible_ids))
      
      print("Downloading datasets...")
      ds$download()
      ds$load()
      
      #streamline features and meta-data
      print("Streamlining features & meta-data...")
      ds$streamline_features(match_to_release = feature_version)
      ds$streamline_metadata(schema = "sfaira")
      
      adata <- ds$adata
      X <- methods::as(methods::as(adata$X, 'CsparseMatrix'), 'dgCMatrix')
      obs <- adata$obs
      var <- adata$var
      return(list(X=X, obs=obs, var=var))
      
    }, error=function(e){
      message(paste0("Could not download all datasets for specified filters."))
      print(e$message)
      return(NULL)
    })
  }, setup_list=setup_list, organisms=organisms, assays=assays, tissues=tissues, force=force, feature_version = pkg.globals$sfaira_streamline_feature_version)

  return(sfaira_dataset)
}


#' Gives an overview of the possible datasets you can use from the sfaira database
#'
#' @param setup_list the sfaira setup; given by \code{\link{setup_sfaira}}
#'
#' @return a dataframe with information on each dataset
#' @export
#'
#' @examples
#' setup_list <- setup_sfaira(basedir=tempdir())
#' all_datasets <- sfaira_overview(setup_list)
sfaira_overview <- function(setup_list){
  
  # create conda environment with sfaira
  proc <- basilisk::basiliskStart(SimBu_env)
  on.exit(basilisk::basiliskStop(proc))
  
  info_list <- basilisk::basiliskRun(proc, function(setup_list){
    sfaira <- reticulate::import("sfaira")
    ds <- sfaira$data$Universe(data_path = setup_list[["rawdir"]],
                               meta_path = setup_list[["metadir"]],
                               cache_path = setup_list[["cachedir"]])
    
    all_datasets <- ds$datasets
    info_list <- lapply(all_datasets, function(x){
      return(list(id=x$id,
                  author=x$author,
                  doi=x$doi,
                  annotated=x$annotated,
                  assay=x$assay_sc,
                  organ=x$organ,
                  organism=x$organism))
    })
    return(info_list)
  }, setup_list=setup_list)

  out <- suppressWarnings(data.table::rbindlist(info_list))
  out$annotated[which(is.na(out$annotated))]<-FALSE
  return(out)
}


# check for given sfaira dataset if all entries have accessible urls to data and metadata 
# Return for all IDs if they can be downloaded
check_sfaira_datasets <- function(dataset){
  
  # create conda environment in which to access sfaira
  proc <- basilisk::basiliskStart(SimBu_env)
  on.exit(basilisk::basiliskStop(proc))
  
  ids_checked <- basilisk::basiliskRun(proc, function(ds){

    ids_list <- lapply(ds$datasets, function(x){
      message(paste('Checking file accessibility for id: ', x$id))
      data_url <- x$download_url_data[[1]]
      meta_url <- x$download_url_meta[[1]]
      
      data_url_endings <- lapply(data_url, function(y){unlist(lapply(strsplit(y,'\\.'), utils::tail, 1))})
      if(!is.null(unlist(data_url))){
        data_url_private <- lapply(data_url, function(y){any(startsWith(y, c('private','manual','syn')))})
      }
      if(!is.null(unlist(meta_url))) {
        meta_url_private <- lapply(meta_url, function(y){any(startsWith(y, c('private','manual','syn')))})
      }
      
      if(is.null(unlist(meta_url)) & all(data_url_endings %in% c('h5','h5ad'))){
        # meta data is contained in h5 obs layer, so no metadata file needs to be present
        meta_accessible <- TRUE
      }else if(is.null(unlist(meta_url))){
        meta_accessible <- FALSE
      }else if(any(meta_url_private)){
        meta_accessible <- FALSE
      }else{
        meta_accessible <- all(RCurl::url.exists(meta_url))
      }
      
      if(is.null(unlist(data_url))){
        data_accessible <- FALSE
      }else if(any(data_url_private)){
        data_accessible <- FALSE
      }else{
        data_accessible <- all(RCurl::url.exists(data_url))
      }
      
      dataset_accessible <- (x$annotated & meta_accessible & data_accessible)
      
      
      return(list(id=x$id,
                  #doi=x$doi,
                  is_annotated=x$annotated,
                  #data_url=data_url,
                  data_accessible=data_accessible,
                  #meta_url=meta_url,
                  meta_accessible=meta_accessible,
                  dataset_accessible=dataset_accessible))
    })
    return(ids_list)
  }, ds=dataset)
  
  out_accessible <- data.table::rbindlist(ids_checked)
  
  ids_with_data <- out_accessible[out_accessible$dataset_accessible == TRUE,]$id
  names(ids_with_data) <- ids_with_data
  
  return(list(out_accessible=out_accessible,
              ids_with_data=ids_with_data))
  
}


# some datasets are currently not working in sfaira and will probably be fixed in the next update
skip_broken_datasets <- function(ids){
  
  print('Checking if IDs contain datasets which are currently known to be broken.')
  
  if('homosapiens_blood_2020_10x3v2_canogamez_001_10.1038/s41467-020-15543-y' %in% ids){
    warning('Removing homosapiens_blood_2020_10x3v2_canogamez_001_10.1038/s41467-020-15543-y from ID list as it is not working in current sfaira release.')
    ids <- ids[ids != 'homosapiens_blood_2020_10x3v2_canogamez_001_10.1038/s41467-020-15543-y']
  }
  
  return(ids)
  
}
