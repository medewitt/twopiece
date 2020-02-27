# Add packages that will be used
usethis::use_package("mclust")
usethis::use_package("ellipse")
usethis::use_package("flexclust")
usethis::use_package("methods")
usethis::use_package("mvtnorm")
usethis::use_package("label.switching")

# add package R
usethis::use_package_doc()

# Add the Pipe
usethis::use_pipe()

# Add the License
usethis::use_gpl3_license()

# Build Website

# https://r-forge.r-project.org/scm/viewvc.php/pkg/twopiece/R/ptp4.R?root=twopiece&view=log

usethis::use_build_ignore("init.R")
usethis::use_build_ignore("Makefile")
usethis::use_build_ignore("docs")
usethis::use_build_ignore("Makefile")
usethis::use_build_ignore(".github")

# CI and Package Architecture
usethis::use_appveyor()
usethis::use_travis()
usethis::use_coverage()
usethis::use_readme_rmd()
# Future Vignettes
# https://rpubs.com/FJRubio/DTP
# https://rpubs.com/FJRubio/twopiece
