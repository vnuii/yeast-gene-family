# This script automates the CAFE analysis by running multiple iterations of the CAFE program with different parameters.
# It takes command line arguments from a Python script and performs the following steps:
# 1. Reads the command line arguments and sets the necessary file paths.
# 2. Finds the directory names in the specified directory and sorts them using the mixedsort function from the gtools package.
# 3. Defines a function to get the likelihood value from the results.txt file in each directory.
# 4. Finds the directory with the maximum likelihood value and assigns it to maxHood_k_number.
# 5. Creates a new directory for storing intermediate files.
# 6. Defines a function to determine the path of the OGcounts_filtered file based on the maximum likelihood k value.
# 7. Checks if the OGcounts_filtered file exists and if not, iterates through the k values until a valid file is found or the k value reaches 0.
# 8. Returns the path of the OGcounts_filtered file.
# 9. Checks if the error_model.txt file exists in the specified directory.
# 10. Runs the CAFE program 10 times with the specified arguments, creating a new directory for each run.
# 11. Checks if the results.txt file exists in each directory and increments the real_number variable if it does.
# 12. Deletes the directory if the results.txt file does not exist.
# 13. Finds the directory names in the forAL_path directory and sorts them using the mixedsort function.
# 14. Finds the directory with the maximum likelihood value and assigns it to maxHood_order.
# 15. Changes the working directory to the directory with the maximum likelihood value.
# 16. Finds the path of the final results.txt file.
# 17. Reads the final results.txt file and assigns the alpha and lambda values.
# 18. Creates a new directory for storing the final output files.
# 19. Runs the CAFE program 10 times with the specified arguments, using the alpha and lambda values from step 17.
# 20. Checks if the final output directory is empty and re-runs the CAFE program if it is.
library(data.table)
library(dplyr)
library(tibble)
library(gtools)

args <- commandArgs(trailingOnly = TRUE)

# get the arguments from the python script
kout_dir_path <- args[1]
OGcounts_filtered_path <- args[2]
tree_path <- args[3]
erroModel_path <- args[4]
forAL_path <- args[5]
OGcounts_path <- args[6]

final_outPath <- args[7]
command <- args[8]
parament_file_path <- args[9]

kout_dir_names <- list.dirs(kout_dir_path) %>%
  c() %>%
  .[-1] %>%
  gtools::mixedsort()

get_all_GammaResults_Likelihood <- function(dir_name) {
  setwd(dir_name)
  GammaResults_FilePath <- list.files(pattern = "results.txt") %>% tail(1)
  if (length(GammaResults_FilePath) != 0) {
    GammaResults_likelihood <- fread(GammaResults_FilePath, header = FALSE, fill = TRUE)[1, 6] %>%
      unname() %>%
      unlist()
    return(GammaResults_likelihood)
  }
}

maxHood_k_number <- lapply(kout_dir_names, get_all_GammaResults_Likelihood) %>%
  unlist() %>%
  c() %>%
  which.max()

# run the CAFE program with the maximum likelihood k
system(paste("mkdir", forAL_path))

test_OGcounts_filtered_path <- function() {
  k_tmp_number <- maxHood_k_number
  if (file.exists(paste0(OGcounts_filtered_path, "/k", k_tmp_number + 1, ".tmp"))) {
    OGcounts_filtered_path_file_r <- paste0(OGcounts_filtered_path, "/k", k_tmp_number + 1, ".tmp")
  } else {
    while (k_tmp_number > 0) {
      k_tmp_number <- k_tmp_number - 1
      if (file.exists(paste0(OGcounts_filtered_path, "/k", k_tmp_number + 1, ".tmp"))) {
        OGcounts_filtered_path_file_r <- paste0(OGcounts_filtered_path, "/k", k_tmp_number + 1, ".tmp")
        break
      }
      if (k_tmp_number == 0) {
        OGcounts_filtered_path_file_r <- parament_file_path
        break
      }
    }
  }
  return(OGcounts_filtered_path_file_r)
}

OGcounts_filtered_path_file <- test_OGcounts_filtered_path()

# used to determine whether the previous task has been executed successfully.
# If the previous task is successful, the script proceeds to the next task.
# If the previous task is not successful, the script deletes the task.
erroModel_path_file <- list.files(erroModel_path, pattern = "error_model.txt", full.names = T)
real_number <- 1
while (real_number <= 10) {
  args <- c(
    "-c", 80,
    "-i", OGcounts_filtered_path_file,
    "-t", tree_path,
    paste0("-e", erroModel_path_file),
    "-p",
    "-k", maxHood_k_number + 1,
    "-o", paste0(forAL_path, "/forAL_", real_number)
  )
  system2(command, args = args)
  files_number <- list.files(paste0(forAL_path, "/forAL_", real_number),
    pattern = "results.txt"
  ) %>%
    length()
  if (files_number == 3) {
    real_number <- real_number + 1
  } else {
    system(paste("rm -r", paste0(forAL_path, "/forAL_", real_number)))
  }
}
# choose the alpha and lambda with the highest likelihood
forAL_dirs <- list.dirs(forAL_path) %>%
  c() %>%
  .[-1] %>%
  .[order(as.numeric(gsub("[^0-9]", "", .)))]

maxHood_order <- lapply(forAL_dirs, get_all_GammaResults_Likelihood) %>%
  unlist() %>%
  c() %>%
  which.max()

setwd(forAL_dirs[maxHood_order])
final_path <- list.files(pattern = "results.txt") %>% tail(1)

alpha <- fread(final_path, header = FALSE, fill = TRUE) %>% .[nrow(.), 2]
lambda <- fread(final_path, header = FALSE, fill = TRUE) %>% .[2, 2]

# -z version
system(paste("mkdir", final_outPath))

for (i in seq_len(10)) {
  args <- c(
    "-c", 80,
    "-i", OGcounts_path,
    "-t", tree_path,
    paste0("-e", erroModel_path_file),
    "-p",
    "-k", maxHood_k_number + 1,
    "-l", lambda,
    "-a", alpha,
    "-o", paste0(final_outPath, "/final_out_", i),
    "-z"
  )
  system2(command, args = args)
  # check if the output folder is empty and reruns the process if it is.
  files_num <- list.files(paste0(final_outPath, "/final_out_", i)) %>% length()
  while (files_num == 0) {
    system(paste("rm -r", paste0(final_outPath, "/final_out_", i)))
    system2(command, args = args)
    files_num <- list.files(paste0(final_outPath, "/final_out_", i)) %>% length()
  }
}
