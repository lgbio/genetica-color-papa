# Set renv path manually before activating
Sys.setenv(RENV_PATHS_ROOT = "../renv")
Sys.setenv(RENV_PATHS_LIBRARY = "../library")

# Then activate renv
renv::activate("../renv")
