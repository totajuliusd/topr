## R CMD check results

0 errors ✓ | 0 warnings ✓ | 0 notes ✓

There were no ERRORs or WARNINGs. 


## R CMD check_rhub results

0 errors ✔ | 0 warnings ✔ | 2 notes ✖

❯ checking for detritus in the temp directory ... NOTE
  Found the following files/directories:
    'lastMiKTeXException'

I got the note lastMiKteXExeption for previous release of topr too, and still do not know what causes it

❯ checking examples ... [8s/21s] NOTE
Examples with CPU (user + system) or elapsed time > 5s
      user system elapsed
topr 4.403  0.135  11.694

I do not understand why the time eplapsed is > 5s. I have tried putting donttest and dontrun after all examples, and have run out of ideas on how to solve this?  