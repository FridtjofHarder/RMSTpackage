# installs RMSTpackage from github, loads, and opens documentation file
if(!"devtools" %in% installed.packages()[,1]) install.packages("devtools") # installs pacakge devtools if needed
library(devtools) #loads package devtools
install_github("fridtjofharder/RMSTpackage") # installs RMST package from devtools
library(RMSTpackage) # loads pacakge RMSTpackage
# ?RMSTpackage # opens help file on RMSTpacakge. NOT WORKING???
?convert_contrast

