if(!"devtools" %in% installed.packages()[,1]) install.packages("devtools") # installs pacakge devtools if needed
library(devtools) #loads package devtools

if(!"npsurvSS" %in% installed.packages()[,1]) install.packages("npsurvSS") # installs pacakge devtools if needed
library(npsurvSS) #loads package devtools

if(!"SSRMST" %in% installed.packages()[,1]) install.packages("SSRMST") # installs pacakge devtools if needed
library(SSRMST) #loads package devtools

install_github("fridtjofharder/RMSTpackage")
library(RMSTpackage)

# access function help pages
?convert_contrast
?compare_sample_size


