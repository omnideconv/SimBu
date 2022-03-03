# "Classical checks"

These comments and suggestions are tailored to the output currently obtained when running good ol' `R CMD check`

- importFrom: some might be missing (aes for ggplot2?, median from stats)
- no visible binding -> might be resolved if you use aes_string instead of aes
- undocumented arguments in documentation object 'calc_scaling_vector' -> check it is documented or the name is matching?
- vignette builder engine must be specified
- the Rproj file is not desired in the final repo on Bioc, just keep it local and put it under .Rbuildignore
- one of the tests is failing?

Some general things on top:
- some info messages are better `message`d than `print`ed (finished simulation)


# BiocCheck specific ones

Output obtained after running `BiocCheck::BiocCheck()` in the folder of the package - some suggestions are next to each item.

```
BiocCheck::BiocCheck()

...

BiocCheck FAILED.
$error
[1] "Packages providing 1 object(s) used in this package should be imported in the NAMESPACE file, otherwise packages importing this package may fail."
--> do as it says and it should be good!
[2] "No biocViews terms found."                                                                              
---> put the biocViews in the DESCRIPTION
[3] "VignetteEngine specified but not in DESCRIPTION. Add the VignetteEngine from the following files to DESCRIPTION:"                            
as above
[4] "At least 80% of man pages documenting exported objects must have runnable examples. The following pages do not:"
not too much to worry, but do put some more?
[5] "Maintainer must register at the support site; visit https://support.bioconductor.org/accounts/signup/ ."                                          
-> easy

$warning
[1] "The following files are over 5MB in size: '.git/objects/pack/pack-1a2a06e11bc2fc33a35436a588f1f8eeefca8577.pack'"
-> clean up some subfolders in git?
[2] "Import scales, matrixStats, testit in NAMESPACE as well as DESCRIPTION."                                         
-> do as it says :)
[3] "Import stats in DESCRIPTION as well as NAMESPACE."                                                               
-> here as well
[4] "Evaluate more vignette chunks."                                                                                  
-> it is pretty much reco
[5] " Avoid T/F variables; If logical, use TRUE/FALSE"                                                                
[6] " Avoid T/F variables; If logical, use TRUE/FALSE (found 38 times)"                                               
-> do that :) 
[7] "Add non-empty \\value sections to the following man pages: man/save_simulation.Rd"                               
-> better put something in

$note
 [1] "Consider clarifying how 3 object(s) are initialized. Maybe they are part of a data set loaded with data(), or perhaps part of an object referenced in with() or within()."      
 
 -> solved w aes_string I guess
 [2] "'LazyData:' in the 'DESCRIPTION' should be set to false or removed"                                                                                                               
 [3] "The Description field in the DESCRIPTION is made up by less than 3 sentences. Please consider\nexpanding this field, and structure it as a full paragraph"                        
 -> say something more
 [4] " 'sessionInfo' not found in vignette(s)"                                                                                     
 -> put it in at the end, as a unnumbered section {-}
 [5] " Avoid sapply(); use vapply()"    
 -> would be better, it is more robust
 [6] " Avoid 1:...; use seq_len() or seq_along()"                                                             -> would be better, it is more robust                                                                          
 [7] " Avoid the use of 'paste' in condition signals" 
 -> it is not really needed
 [8] " Avoid redundant 'stop' and 'warn*' in signal conditions"                                                                                                                         
 [9] "Avoid '<<-' if possible (found 3 times)"                                                                                          -> can you avoid this?                                                
[10] "Avoid 'suppressWarnings'/'*Messages' if possible (found 1 times)"                                                                                                                 
[11] "Recommended function length <= 50 lines."                                                                                                                                         
[12] "Usage of dontrun{} / donttest{} found in man page examples."                                                                                                                      
[13] "Use donttest{} instead of dontrun{}."                                                                                                                                             
[14] "Consider adding a NEWS file, so your package news will be included in Bioconductor release announcements."                                  
-> that is a nice thing to have actually. NEWS.md renders also nicely
[15] "Consider shorter lines; 506 lines (17%) are > 80 characters long."                                                                                                                
[16] "Consider multiples of 4 spaces for line indents, 812 lines(28%) are not."   
-> not so deal-breaker, dont worry here
[17] "Cannot determine whether maintainer is subscribed to the bioc-devel mailing list (requires\nadmin credentials).  Subscribe here: https://stat.ethz.ch/mailman/listinfo/bioc-devel"
```

Running this in your REPL it would give also pointers on where these issues show up - just go ahead and see which spots need some fixes :wink:

