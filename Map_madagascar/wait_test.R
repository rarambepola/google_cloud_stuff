path_output <- Sys.getenv("path_output")
path_input <- Sys.getenv("path_input")

print(paste0("path input: ", path_input))
print(paste0("path output: ", path_output))

wait_start <- Sys.time()
print(paste0("starting waiting at: ", wait_start))
Sys.sleep(17454)
wait_end <- Sys.time()
print(paste0("stopped waiting waiting at: ", wait_end))
print(paste0("waited for ", wait_end - wait_start))

print(paste0("path input: ", path_input))
print(paste0("path output: ", path_output))

print("getting environment variables again")
path_output <- Sys.getenv("path_output")
path_input <- Sys.getenv("path_input")

print(paste0("path input: ", path_input))
print(paste0("path output: ", path_output))
