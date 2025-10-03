# Load required libraries
library(tidyverse)
library(secr)
library(dplyr)

# Load data
rodent <- read.csv("data/mam_pertrapnight.csv")
plot_rodent <- read.csv("data/mam_perplotnight.csv")

# Select necessary columns
plot_rodent <- plot_rodent |> 
  dplyr::select(nightuid, eventID)

# Merge data based on nightuid
rodent <- inner_join(rodent, plot_rodent, by = "nightuid")

# Convert collectDate to POSIXct and extract year, month, day
rodent$collectDate <- as.POSIXct(rodent$collectDate, format = "%Y-%m-%d")
rodent$Year <- format(rodent$collectDate, "%Y")
rodent$Month <- format(rodent$collectDate, "%m")
rodent$Day <- format(rodent$collectDate, "%d")

All_rodent <- rodent |> 
  dplyr::select(nlcdClass, sex, lifeStage, siteID, namedLocation, domainID, trapCoordinate, decimalLatitude, decimalLongitude, Year, Month, Day, eventID, tagID) 

All_rodent_ID <- All_rodent %>%
  mutate(dif = as.factor(abs(round(decimalLatitude - decimalLongitude, 4))),
         TRAPID = paste(trapCoordinate, dif, sep = "_")) |> 
  mutate(TRAPID = as.integer(as.factor(TRAPID))) |> 
  filter(!grepl("X", trapCoordinate))

summary(All_rodent_ID)

# 1. Define the output directory for the trap files
trap_output <- "data/Trap/Trap_all_rodents/"

# Create the directory if it doesn't already exist to prevent errors
dir.create(trap_output, showWarnings = FALSE, recursive = TRUE)

# 2. Define the parameters for coordinate adjustment, as in your original script
# This "jitters" the coordinates to ensure traps on a grid don't overlap perfectly.
coords <- LETTERS[1:10] # Corresponds to grid rows 'A' through 'J'
nums <- 1:10            # Corresponds to grid columns '1' through '10'
offsets <- c(0.0004, 0.0003, 0.0002, 0.0001, 0, -0.0001, -0.0002, -0.0003, -0.0004, -0.0005)

# 3. Get a list of unique site IDs to loop through
sites <- unique(All_rodent_ID$siteID)

# 4. Loop through each site to process and write its trap layout file
for (current_site in sites) {
  
  cat("Processing traps for site:", current_site, "\n")
  
  # First, get a unique list of traps for the current site.
  # We use distinct() to ensure we only have one row per TRAPID.
  site_traps <- All_rodent_ID %>%
    filter(siteID == current_site) %>%
    # This is a robust way to get one unique entry for each trap
    distinct(TRAPID, .keep_all = TRUE) %>%
    select(TRAPID, decimalLatitude, decimalLongitude, trapCoordinate, siteID)
  
  # Proceed only if there are traps for the site
  if (nrow(site_traps) > 0) {
    
    # 5. Apply the coordinate jittering logic from your original script.
    # This loop adjusts the lat/lon for each trap based on its grid position.
    for (j in seq_along(coords)) {
      # Adjust latitude based on the letter prefix (e.g., 'A' in 'A1', 'A2', etc.)
      lat_indices <- site_traps$trapCoordinate %in% paste0(coords[j], nums)
      site_traps$decimalLatitude[lat_indices] <- site_traps$decimalLatitude[lat_indices] + offsets[j]
      
      # Adjust longitude based on the number suffix (e.g., '1' in 'A1', 'B1', etc.)
      lon_indices <- site_traps$trapCoordinate %in% paste0(coords, nums[j])
      site_traps$decimalLongitude[lon_indices] <- site_traps$decimalLongitude[lon_indices] - offsets[j]
    }
    
    # 6. Convert adjusted Lat/Lon to a local Cartesian (x, y) coordinate system.
    # This method simply scales the coordinates as you did in your script.
    # Finally, select the columns in the required SECR format: TrapID, x, y
    trap_data_for_secr <- site_traps %>%
      mutate(
        x = abs(decimalLongitude * 100000),
        y = abs(decimalLatitude * 100000)
      ) %>%
      select(TRAPID, x, y)
    
    # 7. Define the full path for the output file
    output_file <- file.path(trap_output, paste0("trap_data_all_rodents_", current_site, ".txt"))
    
    # 8. Write the processed trap data to a text file
    # The format is critical for SECR's read.capthist function:
    # - Tab-separated (sep = "\t")
    # - No column names (col.names = FALSE)
    # - No row names (row.names = FALSE)
    write.table(
      trap_data_for_secr,
      file = output_file,
      sep = "\t",
      row.names = FALSE,
      col.names = FALSE,
      quote = FALSE
    )
  } else {
    cat("Skipping site:", current_site, "(no trap data)\n")
  }
}

# 1. Define the output directory for the capture history files
cap_output <- "data/Caphist/Caphist_year_occasion/"

# Create the directory if it doesn't exist to avoid errors
dir.create(cap_output, showWarnings = FALSE, recursive = TRUE)

# 2. Get a list of unique site IDs to loop through
sites <- unique(All_rodent_ID$siteID)

# 3. Loop through each site to process and write its data
for (current_site in sites) {
  
  cat("Processing site:", current_site, "\n")
  
  # Prepare the capture history for the current site
  caphist_for_secr <- All_rodent_ID %>%
    
    # Filter for data from the current site only
    filter(siteID == current_site) %>%
    
    # Remove any records that don't have a valid animal identifier (tagID)
    filter(!is.na(tagID) & tagID != "") %>%
    
    # THE KEY STEP: Create the Occasion column by mapping Year to integers.
    # factor(Year) creates a factor with levels sorted chronologically (e.g., 2014, 2015, ...).
    # as.numeric() then converts these levels to integers (1, 2, ...).
    mutate(Occasion = as.numeric(factor(Year))) %>%
    
    # CORRECTED SYNTAX: Create the new 'Session' column using mutate().
    # This is the standard way to add new columns in modern dplyr.
    mutate(Session = 1) %>%
    
    # Select and rename the columns into the standard SECR format.
    select(
      Session,
      ID = tagID,
      Occasion,              # Use the new Occasion column we just created
      TrapID = TRAPID
    ) %>%
    
    # Remove duplicate rows. This is important because the 'occasion' is now a full year.
    # This ensures that an animal caught multiple times in the same trap
    # within the same year is only recorded once for that occasion.
    distinct() %>%
    
    # Arrange the data for readability.
    arrange(ID, Occasion)
  
  # 4. Proceed only if there is data to write
  if (nrow(caphist_for_secr) > 0) {
    # Define the full path for the output file
    output_file <- file.path(cap_output, paste0("caphist_year_occasion_", current_site, ".txt"))
    
    # 5. Write the processed data to a text file in the required SECR format
    write.table(
      caphist_for_secr,
      file = output_file,
      sep = " ",
      row.names = FALSE,
      col.names = FALSE,
      quote = FALSE
    )
  } else {
    cat("Skipping site:", current_site, "(no valid capture data)\n")
  }
}

All_rodent_secr <- list()

for (site in sites) {
  
  # Construct file paths for the current site
  capthist_file <- paste0(cap_output, "caphist_year_occasion_", site, ".txt")
  trap_file <- paste0(trap_output, "trap_data_all_rodents_", site, ".txt")
  
  # Check if both files exist before trying to read them
  if (file.exists(capthist_file) && file.exists(trap_file)) {
    All_rodent_secr[[site]] <- tryCatch(
      read.capthist(
        captfile = capthist_file,
        trapfile = trap_file,
        detector = "multi"
      ),
      error = function(e) {
        # Print a message if an error occurs for a specific site
        cat("Error reading data for site:", site, "\nError was:", conditionMessage(e), "\n")
        return(NULL) # Return NULL on error
      }
    )
  } else {
    cat("Skipping site:", site, "(capture or trap file missing)\n")
  }
}
All_rodent_secr <- All_rodent_secr[!sapply(All_rodent_secr, is.null)]

initialsigma <- list()
for (i in seq_along(All_rodent_secr)) {
  initialsigma[[i]] <- tryCatch(RPSV(All_rodent_secr[[i]], CC = TRUE), 
                                error = function(e) e)
  names(initialsigma)[[i]] <- names(All_rodent_secr)[[i]]
}

initialsigma <- rapply(initialsigma, f = function(x) ifelse(is.nan(x),1,x), how = "replace" )
# Calculate the mean for each list
mean_values <- lapply(initialsigma, function(x) {
  mean(unlist(x), na.rm = TRUE)
})

masks <- list()
for (i in seq_along(All_rodent_secr)) {
  masks[[i]] <- tryCatch(make.mask(traps(All_rodent_secr[[i]]), buffer = initialsigma[[i]], nx = 100, type = 'trapbuffer'),
                         error = function(e) e)
}

fit.HN <- list()
fit.EX <- list()
fit.HR <- list()

# --- Initialization ---
# Get the total number of datasets to process
n_total <- length(All_rodent_secr)

# Pre-allocate lists to store the results
fit.HN <- vector("list", n_total)
fit.EX <- vector("list", n_total)
fit.HR <- vector("list", n_total)

# Start the main timer
start_time <- Sys.time()

# --- Main Loop ---
for (i in seq_along(All_rodent_secr)) {
  
  # Get the name of the current dataset for clear messaging
  dataset_name <- names(All_rodent_secr)[[i]]
  
  # --- 1. Print Progress Header ---
  cat(paste0("\n----------------------------------------------------------\n",
             "Processing ", i, " of ", n_total, ": '", dataset_name, "'\n",
             "----------------------------------------------------------\n"))
  
  # --- 2. Fit Half-Normal (HN) model ---
  cat(" -> Fitting Half-Normal (HN) model... ")
  fit.HN[[i]] <- tryCatch({
    # The 'suppressMessages' keeps the console cleaner by hiding intermediate secr output
    suppressMessages(
      secr.fit(All_rodent_secr[[1]],
               mask = masks[[1]],
               model = D~session,
               details = list(Dlambda = TRUE),
               trace = FALSE,
               detectfn = 'HN')
    )
  }, error = function(e) {
    # Return the error object itself if the fit fails
    return(e) 
  })
  
  # Check if the result is an error and print a status message
  if (inherits(fit.HN[[i]], "error")) {
    cat("FAILED.\n")
  } else {
    cat("DONE.\n")
  }
  names(fit.HN)[[i]] <- dataset_name
  
  
  # --- 3. Fit Exponential (EX) model ---
  cat(" -> Fitting Exponential (EX) model... ")
  fit.EX[[i]] <- tryCatch({
    suppressMessages(
      secr.fit(All_rodent_secr[[i]],
               mask = masks[[i]],
               trace = FALSE,
               detectfn = 'EX')
    )
  }, error = function(e) {
    return(e)
  })
  
  if (inherits(fit.EX[[i]], "error")) {
    cat("FAILED.\n")
  } else {
    cat("DONE.\n")
  }
  names(fit.EX)[[i]] <- dataset_name
        
        
        # --- 4. Fit Hazard-Rate (HR) model ---
        cat(" -> Fitting Hazard-Rate (HR) model... ")
        fit.HR[[i]] <- tryCatch({
          suppressMessages(
            secr.fit(All_rodent_secr[[i]],
                     mask = masks[[i]],
                     trace = FALSE,
                     detectfn = 'HR')
          )
        }, error = function(e) {
          return(e)
        })
        
        if (inherits(fit.HR[[i]], "error")) {
          cat("FAILED.\n")
        } else {
          cat("DONE.\n")
        }
        names(fit.HR)[[i]] <- dataset_name
        
        
        # --- 5. Calculate and Display Time ---
        elapsed_time <- Sys.time() - start_time
        avg_time_per_iter <- elapsed_time / i
        time_remaining <- avg_time_per_iter * (n_total - i)
        
        # Format the time for better readability
        elapsed_formatted <- format(round(elapsed_time, 2))
        remaining_formatted <- format(round(time_remaining, 2))
        
        cat(paste0("\n   Status: Elapsed Time = ", elapsed_formatted, 
                   ", Estimated Time Remaining = ", remaining_formatted, "\n"))
}

# Combine fitted models into a list
fits <- list()
for (i in seq_along(All_rodent_secr)) {
  fits[[i]] <- tryCatch(secrlist(HN = fit.HN[[i]], EX = fit.EX[[i]], HR = fit.HR[[i]]),
                        error = function(e) e)
  names(fits)[[i]] <- names(All_rodent_secr)[[i]]
}

# Calculate AIC for model selection
AIC_fits <- list()
for (i in seq_along(All_rodent_secr)) {
  AIC_fits[[i]] <- tryCatch(AIC(fits[[i]]),
                            error = function(e) e)
  names(AIC_fits)[[i]] <- names(All_rodent_secr)[[i]]
  AIC_fits[[i]]$site <- names(All_rodent_secr)[[i]]
}
AIC_fits
AIC_fits_dat <- Filter(function(x) !inherits(x, "error"), AIC_fits)
AIC_fits_dat <- dplyr::bind_rows(AIC_fits_dat) 

# Predict population density from the best model
fits_predict <- list()
best_mod <- list()
for (i in seq_along(All_rodent_secr)) {
  try(best_mod[[i]] <- names(fits[[i]]) %>% 
        str_detect(rownames(AIC_fits[[i]])) %>%
        keep(fits[[i]], .), silent = TRUE)
}

for (i in seq_along(All_rodent_secr)) {
  fits_predict[[i]] <- tryCatch(as.data.frame(predict(best_mod[[i]][[1]])), error = function(e) e)
  fits_predict[[i]]$siteID <- names(All_rodent_secr)[[i]]
  fits_predict[[i]]$Detection_fn <- names(best_mod[[i]])
  fits_predict[[i]] <- tryCatch(tibble::rownames_to_column(fits_predict[[i]]), error = function(e) e)
}

fits_predict <- Filter(function(x) !inherits(x, "error"), fits_predict)

# Combine predictions into a single data frame and write to CSV
dat_model <- dplyr::bind_rows(fits_predict) 
write.csv(dat_model, "output/population_All.csv")

