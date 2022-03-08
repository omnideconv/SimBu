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

[4] "At least 80% of man pages documenting exported objects must have runnable examples. The following pages do not:"
not too much to worry, but do put some more?


$warning
[1] "The following files are over 5MB in size: '.git/objects/pack/pack-1a2a06e11bc2fc33a35436a588f1f8eeefca8577.pack'"
-> clean up some subfolders in git?
[4] "Evaluate more vignette chunks."                                                                                  
-> it is pretty much reco

$note

 [5] " Avoid sapply(); use vapply()"    
 -> would be better, it is more robust
                                                                    
 [7] " Avoid the use of 'paste' in condition signals" 
 -> it is not really needed
                                          
 [9] "Avoid '<<-' if possible (found 3 times)"                                                                                          -> can you avoid this?                                                

[17] "Cannot determine whether maintainer is subscribed to the bioc-devel mailing list (requires\nadmin credentials).  Subscribe here: https://stat.ethz.ch/mailman/listinfo/bioc-devel"
```

Running this in your REPL it would give also pointers on where these issues show up - just go ahead and see which spots need some fixes :wink:

# Current issues from BiocCheck:
 ```
$error
[1] "At least 80% of man pages documenting exported objects must have runnable examples. The following pages do not:"
--> all remaining functions require special file/data input (for example h5ad files). I guess a runnable example is too much effort in this case..

$warning
[1] "The following files are over 5MB in size: '.git/objects/pack/pack-cbd871970bfe6032cfa700142fbd42c42e9cae78.pack'"

[2] "Evaluate more vignette chunks."
--> how 'bad' is this?

$note
[1] " Avoid sapply(); use vapply()"    
--> kind of tricky in some cases to replace sapply unfortunantely..

[2] " Avoid the use of 'paste' in condition signals"                                                                                                                                   
[3] "Avoid '<<-' if possible (found in 1 files)"  
--> cannot replace without breaking something

[4] "Recommended function length <= 50 lines."                                                                                                                                         
[5] "Usage of dontrun{} / donttest{} found in man page examples."                                                                                                                      
[6] "Use donttest{} instead of dontrun{}."                                                                                                                                             
[7] "Consider shorter lines; 508 lines (15%) are > 80 characters long."                                                                                                                
[8] "Consider multiples of 4 spaces for line indents, 854 lines(26%) are not."                                                                                                         
[9] "Cannot determine whether maintainer is subscribed to the bioc-devel mailing list (requires admin credentials).  Subscribe\nhere: https://stat.ethz.ch/mailman/listinfo/bioc-devel"
```