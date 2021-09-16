import os
import sfaira

ds = None

def sfaira_setup(basedir):
    datadir = os.path.join(basedir, 'raw')
    metadir = os.path.join(basedir, 'meta')
    cachedir = os.path.join(basedir, 'cache')
    
    ds = sfaira.data.Universe(data_path=datadir, meta_path=metadir, cache_path=cachedir)
    if(ds != None):
        return 0, ds  
    else:
        return 1, ds


def sfaira_download_dataset(idx, ds):
    try:
        ds.datasets[idx].download()
        adata = ds.datasets[idx].load() 
        return adata
    except:
        print("Could not download dataset with id {idx}.")
    
    
    
