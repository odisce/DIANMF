.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Welcome to my package DIANMF :) ")
}

# .onLoad <- function(libname, pkgname) {
#   user_permission <- utils::askYesNo("Install miniconda? downloads 50MB and takes time")
#
#   if (isTRUE(user_permission)) {
#     reticulate::install_miniconda()
#   } else {
#     message("You should run `reticulate::install_miniconda()` before using this package")
#   }
# }



# .onLoad <- function(libname, pkgname) {
#   reticulate::use_python('C:/Users/DK273056/AppData/Local/anaconda3/python.exe', required = TRUE)
#
#   # Reset Python environment
#   reticulate::py_run_string("import sys; sys.modules.clear()")
#
#   script_path <- system.file("python/rank_determination_using_concordance.py", package = pkgname)
#
#   if (file.exists(script_path)) {
#     message("Sourcing Python script: ", script_path)
#     reticulate::source_python(script_path)
#     main_env <- reticulate::import_main()
#
#     # List all objects in the main Python environment to verify the functions are loaded
#     python_objects <- ls(main_env)
#     message("Contents of main Python environment: ", paste(python_objects, collapse = ", "))
#
#     # Verify that the specific function is available
#     if (!("perform_rank_determination" %in% python_objects)) {
#       stop("Function 'perform_rank_determination' not found in Python script")
#     }
#   } else {
#     stop("Python script not found: ", script_path)
#   }
# }


# reticulate::py_run_file(system.file("python", "rank_determination_using_concordance.py", package = "DIANMF"))



