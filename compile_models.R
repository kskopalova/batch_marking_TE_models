#===============================================================================
# compile_models.R  -  compile the C++ files once, before running any simulation
#
# Compiled files are specific to the machine they were built on: run this script
# on the machine where the simulations will run. Set force <- FALSE to compile
# only the files that are missing or older than their source.
#===============================================================================
root  <- "."           # folder that contains cpp/
force <- TRUE

cpp.files <- c("random_emigration", "markovian_atbt", "markovian_atbc",
               "markovian_acbt", "markovian_acbc", "phitpt")

library(TMB)
old <- setwd(file.path(root, "cpp")); on.exit(setwd(old))

status <- data.frame(file = cpp.files, compiled = FALSE, loaded = FALSE, message = "")
for (i in seq_along(cpp.files)) {
  f <- cpp.files[i]
  cat("\n---- ", f, " ----\n")
  res <- tryCatch({
    if (!file.exists(paste0(f, ".cpp"))) stop("file not found in cpp/")
    if (force || !file.exists(dynlib(f)) || file.mtime(dynlib(f)) < file.mtime(paste0(f, ".cpp")))
      compile(paste0(f, ".cpp"), "-O1 -g", DLLFLAGS = "")
    status$compiled[i] <- file.exists(dynlib(f))
    dyn.load(dynlib(f))
    status$loaded[i] <- f %in% names(getLoadedDLLs())
    "ok"
  }, error = function(e) conditionMessage(e))
  status$message[i] <- res
}
cat("\n==== Summary ====\n"); print(status, row.names = FALSE)
