setup_py_env <- function(py_env_name, py_location) {
  
  Sys.setenv(PYTORCH_ROCM_ARCH = "gfx1201")
  options(reticulate.conda_binary = py_location)
  
  envs <- reticulate::conda_list()$name
  
  if (!(py_env_name %in% envs)) {
    
    # Create env
    reticulate::conda_create(
      envname = py_env_name,
      python_version = "3.12",
      packages = c("pip", "umap-learn")
    )
    
    # Activate
    reticulate::use_condaenv(py_env_name, required = TRUE)
    
    # Get pip path for this environment
    conda_envs <- reticulate::conda_list()
    env_info <- conda_envs[conda_envs$name == py_env_name, ]
    pip_path <- file.path(dirname(dirname(env_info$python)), "bin", "pip")
    
    # Install PyTorch ROCm using index URL
    system(paste(pip_path, "install torch torchvision --index-url https://download.pytorch.org/whl/rocm6.4"))
    
    # Remove bundled HSA runtime
    torch_lib_path <- reticulate::py_capture_output({
      reticulate::py_run_string("import torch; import os; print(os.path.join(os.path.dirname(torch.__file__), 'lib'))")
    })
    torch_lib_path <- trimws(torch_lib_path)
    system(paste("rm -f", file.path(torch_lib_path, "libhsa-runtime64.so*")))
    
    # Install CellBender
    system(paste(pip_path, "install git+https://github.com/broadinstitute/CellBender.git@refs/pull/420/head"))
    
    # Verify
    torch <- reticulate::import("torch")
    
    message("PyTorch version: ", torch$`__version__`)
    
    gpu_available <- torch$cuda$is_available()
    if (gpu_available) {
      message("GPU detected: ", torch$cuda$get_device_name(0L))
    } else {
      warning("GPU not detected")
    }
    
    return(gpu_available)
    
  } else {
    reticulate::use_condaenv(py_env_name, required = TRUE)
    torch <- reticulate::import("torch")
    gpu_available <- torch$cuda$is_available()
    
    if (gpu_available) {
      message("Using existing env with GPU: ", torch$cuda$get_device_name(0L))
    } else {
      warning("Existing env - GPU not detected")
    }
    
    return(gpu_available)
  }
}