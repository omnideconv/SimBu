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

> it is something you might need to justify. Not a killer thing, but it is good to have runnable examples in general. Maybe a small set, to be inserted as a data component?

$warning
[1] "The following files are over 5MB in size: '.git/objects/pack/pack-cbd871970bfe6032cfa700142fbd42c42e9cae78.pack'"

[2] "Evaluate more vignette chunks."
--> how 'bad' is this?

> it depends. If you have one vignette, and 80% are un-eval'd, well, that's bad. If (as in your case) you have 4 vignettes, and some have uneval'd chunks, it is less problematic.
Probably evaluated on a case by case basis, I'd say.

$note
[1] " Avoid sapply(); use vapply()"    
--> kind of tricky in some cases to replace sapply unfortunantely..

> I know- do try, you "just" need to pass the *expected value type...*

[2] " Avoid the use of 'paste' in condition signals"                                                                                                                                   
[3] "Avoid '<<-' if possible (found in 1 files)"  
--> cannot replace without breaking something

> Ok, got it. "It's a note", and it says "if possible". If it really cannot be avoided, well, it will stay :)

[4] "Recommended function length <= 50 lines."                                                                                                                                         
[5] "Usage of dontrun{} / donttest{} found in man page examples."                                                                                                                      
[6] "Use donttest{} instead of dontrun{}."                                                                                                                                             
[7] "Consider shorter lines; 508 lines (15%) are > 80 characters long."                                                                                                                
[8] "Consider multiples of 4 spaces for line indents, 854 lines(26%) are not."                                                                                                         
[9] "Cannot determine whether maintainer is subscribed to the bioc-devel mailing list (requires admin credentials).  Subscribe\nhere: https://stat.ethz.ch/mailman/listinfo/bioc-devel"
```
