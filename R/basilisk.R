
# dependencies for sfaira and anndata0.7.6 (taken from zellkonverter package)
SimBu_env <- basilisk::BasiliskEnvironment(envname="SimBu_env_0_99_2",
                                            pkgname="SimBu",
                                            packages=c('python==3.9', 
                                                       "anndata==0.7.6",
                                                       "h5py==3.2.1",
                                                       "hdf5==1.10.6",
                                                       "natsort==7.1.1",
                                                       "numpy==1.20.2",
                                                       "packaging==20.9",
                                                       "pandas==1.2.4",
                                                       "scipy==1.6.3",   
                                                       "sqlite==3.35.5"),
                                            pip = c('sfaira==0.3.12', 
                                                    'tables==3.7.0')
)

