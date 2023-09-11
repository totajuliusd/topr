## R CMD check results

0 errors ✓ | 0 warnings ✓ | 0 notes ✓

There were no ERRORs or WARNINGs. 


## R CMD check_rhub results

0 errors ✔ | 0 warnings ✔ | 4 notes ✖

> checking CRAN incoming feasibility ... [11s] NOTE

> * checking examples ... NOTE
Examples with CPU (user + system) or elapsed time > 5s
     user system elapsed
topr 4.84   0.26    5.16

❯ checking for detritus in the temp directory ... NOTE
  Found the following files/directories:
    'lastMiKTeXException'

I got the note lastMiKteXExeption for previous releases of topr too, and still do not know what causes it

❯ checking for non-standard things in the check directory ... NOTE
  Found the following files/directories:
    ''NULL''
    
I havent had this note before and do not understand what the NULL directory or file is.
