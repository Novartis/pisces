getScriptPath <- function() {
  cmd.args <- commandArgs()
  m <- regexpr("(?<=^--file=).+", cmd.args, pe = TRUE)
  script.dir <- dirname(regmatches(cmd.args, m))
  if (length(script.dir) == 0) 
    stop("can't determine script dir: please call the script with Rscript")
  if (length(script.dir) > 1) 
    stop("can't determine script dir: more than one '--file' argument detected")
  return(script.dir)
}
script.dir <- getScriptPath()
setwd(script.dir)
if (!requireNamespace("renv", quietly = TRUE)) { 
    install.packages("renv",  repos='http://cran.us.r-project.org') }
renv::settings$use.cache(FALSE, persist = TRUE)
renv::consent(provided = TRUE)
renv::activate()
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager",  repos='http://cran.us.r-project.org') }
options(repos = BiocManager::repositories())
renv::hydrate(packages=renv::dependencies()$Package)
renv::snapshot()