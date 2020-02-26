# Add packages that will be used
usethis::use_package("mclust")
usethis::use_package("ellipse")
usethis::use_package("flexclust")
usethis::use_package("methods")
usethis::use_package("mvtnorm")
usethis::use_package("label.switching")

# Add the Pipe
usethis::use_pipe()

# Add the License
usethis::use_gpl3_license()

# https://r-forge.r-project.org/scm/viewvc.php/pkg/twopiece/R/ptp4.R?root=twopiece&view=log

usethis::use_build_ignore("init.R")
usethis::use_build_ignore("Makefile")

# Future Vignettes
# https://rpubs.com/FJRubio/DTP
# https://rpubs.com/FJRubio/twopiece
