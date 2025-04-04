# 1. Install reticulate if not already installed
if (!requireNamespace("reticulate", quietly = TRUE)) {
  install.packages("reticulate")
}

# Load reticulate
library(reticulate)

# 2. Ensure pip is installed in the system
system("python -m ensurepip --default-pip")

# 3. Ensure pandas is installed globally (outside virtualenv)
system("python -m pip install --upgrade pip pandas")

# 4. Create the virtual environment if it doesn't exist
if (!dir.exists("~/mi_env")) {
  system("python -m venv ~/mi_env")
}

# 5. Use the virtual environment
reticulate::use_virtualenv("~/mi_env", required = TRUE)

# 6. Install Python packages in the virtual environment
reticulate::py_install(c("pandas", "scikit-learn"))

# 7. Check Python configuration
reticulate::py_config()

# 8. Test if the libraries are correctly installed
reticulate::py_run_string("import pandas; import sklearn; print('Success!')")
