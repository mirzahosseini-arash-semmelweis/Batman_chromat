smart_fread <- function(file, encoding) {
  # Description
  # Helper function required for reading files in batch.eval functions.
  
  # Arguments
  # file      The .csv file to be read.
  # encoding  The recognized encoding of the file.
  
  # Value
  # list      The recognized separator and decimal of the file.
  
  raw_lines <- readLines(con <- file(file, encoding = encoding), n = 10, warn = FALSE)
  sep_counts <- sapply(c("," = ",", ";" = ";", "\t" = "\t", "|" = "|"), function(sep) {
    length(strsplit(raw_lines[1], sep, fixed = TRUE)[[1]]) - 1
  })
  guessed_sep <- names(which.max(sep_counts))
  all_fields <- unlist(strsplit(paste(raw_lines, collapse = guessed_sep), split = guessed_sep, fixed = TRUE))
  number_like <- grep("^[0-9]+[.,][0-9]+$", all_fields, value = TRUE)
  dot_count <- sum(grepl("\\.", number_like))
  comma_count <- sum(grepl(",", number_like))
  guessed_dec <- if (comma_count > dot_count) "," else "."
  list(sep = guessed_sep, dec = guessed_dec)
}

plot.runs <- function(path) {
  # Description
  # This function plots all chromatograms into new subdirectory as .png files.
  
  # Arguments
  # path      The path where files are found.
  
  # Value
  # plots     Folder is created in path with plots of chromatograms.
  
  # load required packages
  require(stringr)
  require(tools)
  require(readr)
  require(ggplot2)
  require(fs)
  require(data.table)
  
  # reading files
  setwd(path)
  csv_files <- list.files(pattern = "(?i)\\.csv$", full.names = TRUE)
  csv_files <- csv_files[!tolower(basename(csv_files)) %in% c("summary_data.csv", "peak_table.csv")]
  if (length(csv_files) == 0) {stop("No .csv files found in path directory.")}
  
  # creating output initialization
  output_dir <- "plots"
  dir_create(output_dir)
  
  # evaluating files
  for (file in csv_files) {
    base_name <- tools::file_path_sans_ext(basename(file))
    base_name_parts <- str_split(base_name, "_", simplify = TRUE)[1:4]
    comp_name <- base_name_parts[1]
    col_name <- base_name_parts[2]
    flow <- as.numeric(str_replace(base_name_parts[3], "^(\\d)(\\d+)$", "\\1.\\2"))
    temp <- as.numeric(base_name_parts[4])
    title <- paste(paste(comp_name, col_name, flow, temp, sep = ", "), "(comp, col, flow ml/min, °C)")
    encoding <- readr::guess_encoding(file)$encoding
    sep <- smart_fread(file, encoding = encoding)$sep
    dec <- smart_fread(file, encoding = encoding)$dec
    df <- read.csv(file, header = FALSE, sep = sep, dec = dec, fileEncoding = encoding)
    if (!is.numeric(df[, 1])) {
      df <- read.csv(file, header = TRUE, sep = sep, dec = dec, fileEncoding = encoding)
    }
    
    data <- data.frame(df[, 1:2])
    colnames(data) <- c("t", "A")
    
    plot <- ggplot(data = data) +
      geom_line(aes(x = t, y = A)) +
      labs(title = title, x = "Time (min)", y = "Intensity") +
      scale_x_continuous(breaks = seq(round(min(data$t), 0), round(max(data$t), 0), by = 1),
                         guide = guide_axis(check.overlap = TRUE)) +
      theme_bw()
    out_file <- file.path(output_dir, paste0(base_name, ".png"))
    ggsave(out_file, plot = plot, width = 6, height = 4, dpi = 300)
  }
}

fit.batch.kue <- function(data, threshold = 0.3, minSNR = 10, minpeakdist = 1,
                          peak_table = FALSE, file, t0_given = FALSE,
                          microrev = FALSE, A0 = 0.5) {
  # Description
  # Helper function to fit chromatography interconversion kinetic constant from
  # unified equation using a two-step iterative minimization.
  
  # Arguments
  # data        A data frame, data table or matrix with at least 2 columns, in
  #             which the first contains numeric indicating time, the second
  #             column indicating measured chromatography signal. All other
  #             columns will be disregarded.
  # threshold   The minimum peak threshold passed to the `pracma::findpeaks` algorithm.
  #             Default is 0.3. Adjust lower value for noisier or broader peaks.
  # minSNR      The minimum SNR a peak must have to be included as peak. The
  #             baseline standard deviation is estimated from the beginning and
  #             portions of the chromatogram. Default value is 10. Adjust value
  #             noisier chromatograms.
  #             Default is 50.
  # minpeakdist The minimum peak distance enforced, expressed as percentage point
  #             of total run time. Default is 1%. Warning, during the algorithm
  #             closer than `minpeakdistance` peaks will be excluded from lowest
  #             intensity first.
  # peak_table  Logical, indicating whether the retention times of the dead peak
  #             and two Batman peaks are given in a guide. This guide must be in
  #             the same path named "peak_table.csv" and must contain the following
  #             columns: file_name, t_M, t_A, t_B. The file_name must be the vector
  #             of file names as appearing in the output "summary_data.csv".
  #             The t_M, t_A, t_B vectors contain retention time estimates where
  #             the algorithm will search for the peaks. If t_A, t_B values are
  #             empty in the "peak_table.csv" file (because of coalesced peaks),
  #             the algorithm will skip unified equation calculation for that
  #             chromatogram. Providing t_M for these coalesced peaks is still
  #             advised, since accurate dead peak identification is required for
  #             stochastic modeling. Default is FALSE.
  # file        The file of the chromatogram being processed. Passed on from the
  #             batch.eval.kue() function.
  # t0_given    Logical, indicating whether dead peaks are included in the chromatograms.
  #             If TRUE, the algorithm will not look for dead peaks and take t_M
  #             values in "peak_table.csv" as given. Default is FALSE.
  # microrev    Logical, indicating whether the principle of microscopic reversibility
  #             applies. Default is FALSE, which takes into account that the rate
  #             constants of conversion of A and B are rendered differently.
  # A0          Fraction of the first eluted isomer in the injected sample.
  #             Default is 0.5.
  
  # Value
  # result      The fitted peak picking parameters and kinetic constants (kue_f and
  #             kue_r for the forward and reverse direction) are returned.
  
  # References
  # Trapp, O. (2006) Unified Equation for Access to Rate Constants of First-Order
  # Reactions in Dynamic and On-Column Reaction Chromatography. Anal. Chem. 78, 189-198.
  
  # load required packages
  require(pracma)
  require(dplyr)
  require(stringr)
  
  # error handling
  if (missing(data)) stop("No data provided.")
  if (!is.numeric(threshold)) stop("Argument 'threshold' must be numeric.")
  if (!is.numeric(minSNR)) stop("Argument 'minSNR' must be numeric.")
  if (!is.numeric(minpeakdist)) stop("Argument 'minpeakdist' must be numeric.")
  if (!is.logical(microrev)) stop("Argument 'microrev' must be boolean")
  if (!is.numeric(A0)) stop("Argument 'A0' must be numeric.")
  if (A0 > 1 | A0 < 0) stop("Argument 'A0' must be between 0 and 1.")
  if (peak_table == FALSE & t0_given == TRUE) stop("Peak table must be provided to fetch t0 values.")
  
  # use first 2 columns of dataframe
  data <- data.frame(data[, 1:2])
  colnames(data) <- c("t", "A")
  
  # baseline correction
  all_peaks <- findpeaks(data$A, threshold = 0.3) # find all peaks (rough estimate)
  first_peak_start <- min(min(all_peaks[, 3]), 0.02*nrow(data))
  if (first_peak_start == 1) {message("Peak found at first time point; baseline might be unreliable.")}
  last_peak_ends <- max(max(all_peaks[, 4]), 0.98*nrow(data))
  if (last_peak_ends == nrow(data)) {message("Peak found at last time point; baseline might be unreliable.")}
  baseline_sample <- data[c(1:first_peak_start, last_peak_ends:nrow(data)), ]
  baseline_sd <- sd(baseline_sample$A)
  baseline <- data[abs(data$A - mean(baseline_sample$A)) < 3*baseline_sd, ]
  fit <- lm(A ~ poly(t, 3, raw = TRUE), data = baseline)
  baseline_pred <- predict(fit, newdata = data, type = "response")
  data$A <- data$A - baseline_pred
  
  # peak find
  peaks <- findpeaks(data$A,
                     threshold = threshold,
                     minpeakheight = minSNR*baseline_sd,
                     minpeakdistance = (minpeakdist/100)*nrow(data))
  peaks <- data.frame(peaks)
  if (nrow(peaks) == 0) {
    stop("No peaks found.")
  }
  colnames(peaks) <- c("A", "ti", "ti_start", "ti_end")
  peaks <- peaks %>% arrange(ti)
  peaks$t <- data$t[peaks$ti]
  # find peaks on chromatogram
  if (!peak_table) {
    t0peak <- peaks[1, ]
    t_M <- t0peak$t
    peaks <- peaks[-1, ]
    npeaks <- nrow(peaks)
    if (npeaks == 0) {
      stop("No peaks found other than dead time.")
    }
    if (npeaks == 1) {
      stop("Only one peak found; perform stochastic modeling.")
    }
    peaks <- peaks %>% filter(row_number() %in% c(1, n()))
  }
  # find peaks based on peak_table
  if (peak_table) {
    peak_tab <- read.csv("peak_table.csv", header = TRUE, sep = ",", dec = ".")
    peak_values <- peak_tab %>% filter(file_name == basename(file))
    if (anyNA(peak_values)) {stop("Peak values not given (completely) in peak_tab.")}
    # either t0 is given or needs to be found
    if (t0_given) {
      if (nrow(peaks) < 2) {stop("Not enough peaks found; perform stochastic modeling.")}
      t_M <- peak_values$t_M
      idx <- c(
        which.min(abs(peaks$t - peak_values$t_A)),
        which.min(abs(peaks$t - peak_values$t_B))
      )
      peaks <- peaks[idx,]
      t0peak <- data.frame(A = NA, ti = NA, ti_start = NA, ti_end = NA, t = NA)
    } else {
      if (nrow(peaks) < 3) {stop("Not enough peaks found; perform stochastic modeling.")}
      t0peak <- peaks %>% slice_min(abs(peaks$t - peak_values$t_M), n = 1, with_ties = FALSE)
      t_M <- t0peak$t
      idx <- c(
        which.min(abs(peaks$t - peak_values$t_A)),
        which.min(abs(peaks$t - peak_values$t_B))
      )
      peaks <- peaks[idx,]
    }
  }
  
  # evaluate chromatographic parameters
  if (!anyNA(t0peak)) {
    t_Mi <- t0peak$ti
    h_M <- t0peak$A
    cutoff_idx <- max(10, round(0.5*nrow(data)/max(data$t))) # 0.5 min with security minimum
    idx_left <- (t_Mi - cutoff_idx):t_Mi
    idx_right <- t_Mi:(t_Mi + cutoff_idx)
    w_M.i <- which.min(abs(data$A[idx_left] - h_M/2)) + t_Mi - cutoff_idx - 1
    w_M.i2 <- which.min(abs(data$A[idx_right] - h_M/2)) + t_Mi - 1
    w_M <- data$t[w_M.i2] - data$t[w_M.i]
    s_M <- w_M/sqrt(8*log(2))
  } else {
    s_M <- NA
  }
  thresh <- max(abs(diff(data$A)))
  t_Ai <- peaks[1, 2]
  t_Bi <- peaks[2, 2]
  t_A <- data$t[t_Ai]
  t_B <- data$t[t_Bi]
  t_As <- t_A*60
  t_Bs <- t_B*60
  h_A <- peaks[1, 1]
  h_B <- peaks[2, 1]
  h_mid <- data$A[mean(c(t_Ai, t_Bi))]
  if (h_A >= h_B) {
    h_p <- 100*h_mid/h_A
  } else {
    h_p <- 100*h_mid/h_B
  }
  idx_leftA <- 1:t_Ai
  idx_rightA <- t_Ai:nrow(data)
  w_A.i <- idx_leftA[max(which(abs(data$A[idx_leftA] - h_A/2) <= thresh))]
  if (data$A[w_A.i] > h_A/2) {
    while(data$A[w_A.i] > h_A/2) {w_A.i <- w_A.i - 1}
    w_A.i <- ifelse(abs(data$A[w_A.i] - h_A/2) > abs(data$A[w_A.i + 1] - h_A/2), w_A.i + 1, w_A.i)
  } else {
    while(data$A[w_A.i] < h_A/2) {w_A.i <- w_A.i + 1}
    w_A.i <- ifelse(abs(data$A[w_A.i] - h_A/2) > abs(data$A[w_A.i - 1] - h_A/2), w_A.i - 1, w_A.i)
  }
  w_A.i2 <- idx_rightA[min(which(abs(data$A[idx_rightA] - h_A/2) <= thresh))]
  if (h_mid < h_A/2) {
    if (data$A[w_A.i2] > h_A/2) {
      while(data$A[w_A.i2] > h_A/2) {w_A.i2 <- w_A.i2 + 1}
      w_A.i2 <- ifelse(abs(data$A[w_A.i2] - h_A/2) > abs(data$A[w_A.i2 - 1] - h_A/2), w_A.i2 - 1, w_A.i2)
    } else {
      while(data$A[w_A.i2] < h_A/2) {w_A.i2 <- w_A.i2 - 1}
      w_A.i2 <- ifelse(abs(data$A[w_A.i2] - h_A/2) > abs(data$A[w_A.i2 + 1] - h_A/2), w_A.i2 + 1, w_A.i2)
    }
  }
  w_A <- ifelse(h_mid >= h_A/2, 2*abs(data$t[w_A.i] - t_A), data$t[w_A.i2] - data$t[w_A.i])
  w_As <- w_A*60
  w_A.1 <- data$t[w_A.i]
  w_A.2 <- ifelse(h_mid >= h_A/2, w_A + w_A.1, data$t[w_A.i2])
  idx_leftB <- 1:t_Bi
  idx_rightB <- t_Bi:nrow(data)
  w_B.i <- idx_leftB[max(which(abs(data$A[idx_leftB] - h_B/2) <= thresh))]
  if (h_mid < h_B/2) {
    if (data$A[w_B.i] > h_B/2) {
      while(data$A[w_B.i] > h_B/2) {w_B.i <- w_B.i - 1}
      w_B.i <- ifelse(abs(data$A[w_B.i] - h_B/2) > abs(data$A[w_B.i + 1] - h_B/2), w_B.i + 1, w_B.i)
    } else {
      while(data$A[w_B.i] < h_B/2) {w_B.i <- w_B.i + 1}
      w_B.i <- ifelse(abs(data$A[w_B.i] - h_B/2) > abs(data$A[w_B.i - 1] - h_B/2), w_B.i - 1, w_B.i)
    }
  }
  w_B.i2 <- idx_rightB[min(which(abs(data$A[idx_rightB] - h_B/2) <= thresh))]
  if (data$A[w_B.i2] > h_B/2) {
    while(data$A[w_B.i2] > h_B/2) {w_B.i2 <- w_B.i2 + 1}
    w_B.i2 <- ifelse(abs(data$A[w_B.i2] - h_B/2) > abs(data$A[w_B.i2 - 1] - h_B/2), w_B.i2 - 1, w_B.i2)
  } else {
    while(data$A[w_B.i2] < h_B/2) {w_B.i2 <- w_B.i2 - 1}
    w_B.i2 <- ifelse(abs(data$A[w_B.i2] - h_B/2) > abs(data$A[w_B.i2 + 1] - h_B/2), w_B.i2 + 1, w_B.i2)
  }
  w_B <- ifelse(h_mid >= h_B/2, 2*abs(data$t[w_B.i2] - t_B), data$t[w_B.i2] - data$t[w_B.i])
  w_Bs <- w_B*60
  w_B.2 <- data$t[w_B.i2]
  w_B.1 <- ifelse(h_mid >= h_B/2, w_B.2 - w_B, data$t[w_B.i])
  s_As <- w_As/sqrt(8*log(2))
  s_Bs <- w_Bs/sqrt(8*log(2))
  A_spl <- spline(x = data$t[peaks[1, 3]:mean(c(t_Ai, t_Bi))],
                  y = data$A[peaks[1, 3]:mean(c(t_Ai, t_Bi))],
                  n = 2000)
  A <- sum(A_spl$y)
  B_spl <- spline(x = data$t[mean(c(t_Ai, t_Bi)):peaks[2, 4]],
                  y = data$A[mean(c(t_Ai, t_Bi)):peaks[2, 4]],
                  n = 2000)
  B <- sum(B_spl$y)
  K <- A/B
  N <- mean(5.545*(c(t_A, t_B)/c(w_A, w_B))^2)
  
  # calculate kue
  calc.kue1 <- function(t_As, t_Bs, s_As, s_Bs, N, h_p, A0) {
    -(1/t_As)*(
      log(
        (100*(1 - A0) + A0*(100 - h_p*(1 + sqrt(2/(pi*N)))))/(t_Bs - t_As)
      ) - log(
        (1 - A0)*(
          (h_p*exp(-(t_Bs - t_As)^2/(2*s_Bs^2)) - 100*exp(-(t_Bs - t_As)^2/(8*s_Bs^2)))/(s_Bs*sqrt(2*pi)) + 100/(t_Bs - t_As)
        ) - A0*(
          (100*exp(-(t_Bs - t_As)^2/(8*s_As^2)) - h_p)/(s_As*sqrt(2*pi)) + (h_p*(1 + sqrt(2/(pi*N))) - 100)/(t_Bs - t_As)
        )
      )
    )
  }
  calc.kue2 <- function(t_As, t_Bs, s_As, s_Bs, N, h_p, A0) {
    -(1/t_As)*(
      log(
        (100*A0 + (1 - A0)*(100 - h_p*(1 - sqrt(2/(pi*N)))))/(t_Bs - t_As)
      ) - log(
        A0*(
          (h_p*exp(-(t_Bs - t_As)^2/(2*s_As^2)) - 100*exp(-(t_Bs - t_As)^2/(8*s_As^2)))/(s_As*sqrt(2*pi)) + 100/(t_Bs - t_As)
        ) - (1 - A0)*(
          (100*exp(-(t_Bs - t_As)^2/(8*s_Bs^2)) - h_p)/(s_Bs*sqrt(2*pi)) + (h_p*(1 - sqrt(2/(pi*N))) - 100)/(t_Bs - t_As)
        )
      )
    )
  }
  
  if (h_A >= h_B) {
    kue_f <- calc.kue1(t_As, t_Bs, s_As, s_Bs, N, h_p, A0)
  } else {
    kue_f <- calc.kue2(t_As, t_Bs, s_As, s_Bs, N, h_p, A0)
  }
  if (microrev) {
    kue_r <- kue_f*K*t_As/t_Bs
  } else {kue_r <- kue_f*K}
  
  # export results
  result <- list(t_M = t_M,
                 s_M = s_M,
                 t_Ai = t_Ai,
                 t_Bi = t_Bi,
                 t_A = t_A,
                 t_B = t_B,
                 t_As = t_As,
                 t_Bs = t_Bs,
                 h_A = h_A,
                 h_B = h_B,
                 h_mid = h_mid,
                 h_p = h_p,
                 w_A.i = w_A.i,
                 w_A.i2 = w_A.i2,
                 w_A = w_A,
                 w_As = w_As,
                 w_A.1 = w_A.1,
                 w_A.2 = w_A.2,
                 w_B.i = w_B.i,
                 w_B.i2 = w_B.i2,
                 w_B = w_B,
                 w_Bs = w_Bs,
                 w_B.1 = w_B.1,
                 w_B.2 = w_B.2,
                 s_As = s_As,
                 s_Bs = s_Bs,
                 A = A,
                 B = B,
                 K = K,
                 N = N,
                 kue_f = kue_f,
                 kue_r = kue_r,
                 A_bc = data$A)
  return(result)
}

batch.eval.kue <- function(path, threshold = 0.3, minSNR = 10, minpeakdist = 1,
                           peak_table = FALSE, t0_given = FALSE,
                           microrev = FALSE, A0 = 0.5) {
  # Description
  # This function performs batch evaluation of the unified equation on chromatography
  # files provided in .csv format. It uses the algorithm of fit.batch.kue().
  
  # Arguments
  # path  The path where files are found.
  # All other arguments have the same defaults and are passed to fit.batch.kue().
  
  # Value
  # The function call generates a new folder in the path and saves the chromatography
  # plots as .png files. A summary_log.txt is generated in the path with Success/Error
  # status and error messages. A summary_data.csv is also generated in the path
  # with fitted peak parameters, kinetic constants and calculated thermodynamic parameters.
  
  # References
  # Trapp, O. (2006) Unified Equation for Access to Rate Constants of First-Order
  # Reactions in Dynamic and On-Column Reaction Chromatography. Anal. Chem. 78, 189-198.
  #
  # Trapp, O. (2008) A novel software tool for high throughput measurements of
  # interconversion barriers: DCXplorer. Journal of Chromatography B, 875, 42-47.
  
  # load required packages
  require(fs)
  require(readr)
  require(data.table)
  require(tools)
  require(ggplot2)
  require(stringr)
  require(dplyr)
  require(broom)
  
  # reading files
  setwd(path)
  csv_files <- list.files(pattern = "(?i)\\.csv$", full.names = TRUE)
  csv_files <- csv_files[!tolower(basename(csv_files)) %in% c("summary_data.csv", "peak_table.csv")]
  if (length(csv_files) == 0) {stop("No .csv files found in path directory.")}
  
  # creating output initialization
  output_dir <- "plots_kue"
  dir_create(output_dir)
  log_entries <- list()
  summary_rows <- list()
  
  # evaluating files
  for (file in csv_files) {
    if (t0_given) {
      peak_tab <- read.csv("peak_table.csv", header = TRUE, sep = ",", dec = ".")
      peak_values <- peak_tab %>% filter(file_name == basename(file))
    }
    base_name <- tools::file_path_sans_ext(basename(file))
    base_name_parts <- str_split(base_name, "_", simplify = TRUE)[1:4]
    comp_name <- base_name_parts[1]
    col_name <- base_name_parts[2]
    flow <- as.numeric(str_replace(base_name_parts[3], "^(\\d)(\\d+)$", "\\1.\\2"))
    Temp <- as.numeric(base_name_parts[4]) + 273.15
    title <- paste(paste(comp_name, col_name, flow, Temp - 273.15, sep = ", "), "(comp, col, flow ml/min, °C)")
    summary_row <- data.frame(file = basename(file),
                              comp_name = comp_name,
                              col_name = col_name,
                              flow = flow,
                              Temp = Temp)
    log_entry <- list(
      file = basename(file),
      kue_status = "Success",
      kue_error = NA
    )
    tryCatch({
      encoding <- readr::guess_encoding(file)$encoding
      sep <- smart_fread(file, encoding = encoding)$sep
      dec <- smart_fread(file, encoding = encoding)$dec
      df <- read.csv(file, header = FALSE, sep = sep, dec = dec, fileEncoding = encoding)
      if (!is.numeric(df[, 1])) {
        df <- read.csv(file, header = TRUE, sep = sep, dec = dec, fileEncoding = encoding)
      }
      colnames(df) <- c("V1", "V2")
      # result table
      result <- fit.batch.kue(data = df, threshold, minSNR, minpeakdist,
                              peak_table, file, t0_given, microrev, A0)
      df$V2 <- result$A_bc
      result <- result[-length(result)]
      Gue_f <- -log(result$kue_f/(1.380662e-23*Temp*0.5/6.626176e-34))*8.31441*Temp
      Gue_r <- -log(result$kue_r/(1.380662e-23*Temp*0.5/6.626176e-34))*8.31441*Temp
      result$Gue_f <- Gue_f
      result$Gue_r <- Gue_r
      table <- data.frame(result)
      table <- table %>% select(t_M,
                                s_M,
                                t_A,
                                t_B,
                                h_A,
                                h_B,
                                h_p,
                                w_A,
                                w_B,
                                A,
                                B,
                                K,
                                N,
                                kue_f,
                                kue_r,
                                Gue_f,
                                Gue_r)
      summary_row <- cbind(
        summary_row,
        table
      )
      # chromatography plot
      h_max <- max(result$h_A, result$h_B)
      t_max <- result$t_B
      plot <- ggplot(data = df) +
        geom_line(aes(x = V1, y = V2)) +
        labs(title = title, x = "Time (min)", y = "Intensity") +
        annotate("segment", x = result$t_A, y = result$h_A - 0.05*h_max,
                 xend = result$t_A, yend = result$h_A + 0.05*h_max,
                 col = "red3") +
        annotate("segment", x = result$t_B, y = result$h_B - 0.05*h_max,
                 xend = result$t_B, yend = result$h_B + 0.05*h_max,
                 col = "red3") +
        annotate("segment", x = result$t_M, y = -0.03*h_max,
                 xend = result$t_M, yend = 0,
                 col = "red3") +
        annotate("text", x = (result$t_A + 0.1*t_max), y = result$h_A,
                 label = paste(round(result$t_A, 1), " min"), col = "red3") +
        annotate("text", x = (result$t_B + 0.1*t_max), y = result$h_B,
                 label = paste(round(result$t_B, 1), " min"), col = "red3") +
        annotate("text", x = result$t_M, y = -0.05*h_max,
                 label = paste("italic(t)[M]"), col = "red3", parse = TRUE) +
        annotate("segment", x = mean(c(result$t_A, result$t_B)), y = 0,
                 xend = mean(c(result$t_A, result$t_B)), yend = result$h_mid,
                 col = "dodgerblue3") +
        annotate("text", x = mean(c(result$t_A, result$t_B)), y = result$h_mid + 0.05*h_max,
                 label = round(result$h_mid, 2), col = "dodgerblue3") +
        annotate("segment", x = result$w_A.1, y = result$h_A/2,
                 xend = result$w_A.2, yend = result$h_A/2, col = "dodgerblue3") +
        annotate("text", x = result$w_A.2 + 0.1*t_max, y = result$h_A/2, col = "dodgerblue3",
                 label = paste(round(result$w_A, 1), " min")) +
        annotate("segment", x = result$w_B.1, y = result$h_B/2,
                 xend = result$w_B.2, yend = result$h_B/2, col = "dodgerblue3") +
        annotate("text", x = result$w_B.2 + 0.1*t_max, y = result$h_B/2, col = "dodgerblue3",
                 label = paste(round(result$w_B, 1), " min")) +
        annotate("text", x = mean(c(result$t_A, result$t_B)), y = 0,
                 label = paste("italic(k)[1]^ue ==", signif(result$kue_f, digits = 3)),
                 col = "black", parse = TRUE) +
        theme_bw()
      out_file <- file.path(output_dir, paste0(base_name, ".png"))
      ggsave(out_file, plot = plot, width = 6, height = 4, dpi = 300)
    }, error = function(e) {
      log_entry$kue_status <<- "Error"
      log_entry$kue_error <<- conditionMessage(e)
      summary_row <<- cbind(summary_row,
                            data.frame(
                              t_M = ifelse(t0_given, peak_values$t_M, NA),
                              s_M = NA,
                              t_A = NA,
                              t_B = NA,
                              h_A = NA,
                              h_B = NA,
                              h_p = NA,
                              w_A = NA,
                              w_B = NA,
                              A = NA,
                              B = NA,
                              K = NA,
                              N = NA,
                              kue_f = NA,
                              kue_r = NA,
                              Gue_f = NA,
                              Gue_r = NA
                            ))
    })
    summary_rows[[length(summary_rows) + 1]] <- summary_row
    log_entries[[length(log_entries) + 1]] <- log_entry
  }
  
  # create log and summary files
  log_df <- rbindlist(log_entries, fill = TRUE)
  summary_df <- rbindlist(summary_rows, fill = TRUE)
  
  #add linreg Eyring-Polanyi
  summary_df <- summary_df %>%
    group_by(comp_name, col_name) %>%
    group_modify(~ {
      valid_data <- .x %>% filter(!is.na(kue_f) & !is.na(Temp)) %>%
        group_by(Temp) %>%
        summarise(mean_kue_f = mean(kue_f, na.rm = TRUE),
                  mean_kue_r = mean(kue_r, na.rm = TRUE),
                  .groups = "drop")
      EP_f_intercept <- NA_real_
      EP_f_slope <- NA_real_
      EP_f_intercept_se <- NA_real_
      EP_f_slope_se <- NA_real_
      EP_f_r2 <- NA_real_
      EP_r_intercept <- NA_real_
      EP_r_slope <- NA_real_
      EP_r_intercept_se <- NA_real_
      EP_r_slope_se <- NA_real_
      EP_r_r2 <- NA_real_
      if (nrow(valid_data) >= 2) {
        EP_f_model <- lm(I(log(mean_kue_f/Temp)) ~ I(1/Temp), data = valid_data)
        EP_f_tidy_model <- tidy(EP_f_model)
        EP_f_glance_model <- glance(EP_f_model)
        EP_f_intercept <- EP_f_tidy_model$estimate[EP_f_tidy_model$term == "(Intercept)"]
        EP_f_slope <- EP_f_tidy_model$estimate[EP_f_tidy_model$term == "I(1/Temp)"]
        EP_f_intercept_se <- EP_f_tidy_model$std.error[EP_f_tidy_model$term == "(Intercept)"]
        EP_f_slope_se <- EP_f_tidy_model$std.error[EP_f_tidy_model$term == "I(1/Temp)"]
        EP_f_r2 <- EP_f_glance_model$r.squared
        EP_r_model <- lm(I(log(mean_kue_r/Temp)) ~ I(1/Temp), data = valid_data)
        EP_r_tidy_model <- tidy(EP_r_model)
        EP_r_glance_model <- glance(EP_r_model)
        EP_r_intercept <- EP_r_tidy_model$estimate[EP_r_tidy_model$term == "(Intercept)"]
        EP_r_slope <- EP_r_tidy_model$estimate[EP_r_tidy_model$term == "I(1/Temp)"]
        EP_r_intercept_se <- EP_r_tidy_model$std.error[EP_r_tidy_model$term == "(Intercept)"]
        EP_r_slope_se <- EP_r_tidy_model$std.error[EP_r_tidy_model$term == "I(1/Temp)"]
        EP_r_r2 <- EP_r_glance_model$r.squared
      }
      .x %>%
        mutate(dHue_f = -EP_f_slope*8.314,
               dHue_f_se = EP_f_slope_se*8.314,
               dSue_f = (EP_f_intercept - log(2.084e10))*8.314,
               dSue_f_se = EP_f_intercept_se*8.314,
               EPue_f_r2 = EP_f_r2,
               dHue_r = -EP_r_slope*8.314,
               dHue_r_se = EP_r_slope_se*8.314,
               dSue_r = (EP_r_intercept - log(2.084e10))*8.314,
               dSue_r_se = EP_r_intercept_se*8.314,
               EPue_r_r2 = EP_r_r2)
    }) %>%
    ungroup()
  
  # export log and summary files
  fwrite(log_df, "summary_log.txt", sep = "\t", na = "")
  fwrite(summary_df, "summary_data.csv", sep = ",", dec = ".")
  message("Processing complete.")
}

Batman <- function(np, t_run, n_A, n_B, tau_A, tau_B, a, b, alpha) {
  # Description
  # Helper function to simulate Batman peaks given a set of chromatography and
  # interconversion parameters.
  
  # Arguments
  # np            Number of data points to generate.
  # t_run         Total run time in minutes.
  # n_A, n_B      Average number of the adsorption–desorption steps for the two peaks.
  # tau_A, tau_B  Sojourn times for the two peaks.
  # a, b          Number of forward and reverse direction interconversions during the run.
  # alpha         Fraction of the first eluted isomer in the injected sample.
  
  # Value
  # result        List of peak data frames (run time, peak signals before and after convolution),
  #               and state probabilities.
  
  # References
  # Felinger, A. (2008) Molecular dynamic theories in chromatography. Journal of
  # Chromatography A, 1184, 20-41.
  #
  # Sepsey, A., Németh, D. R., Németh, G., Felinger, A. (2018) Rate constant determination of 
  # interconverting enantiomers by chiral chromatography using a stochastic model.
  # Journal of Chromatography A, 1564, 155-162.
  
  # load required packages
  require(pracma)
  
  np <- np - 1 # ensures number of points in exported data is np
  
  # ensure peak A elutes first
  peakA <- n_A*tau_A
  peakB <- n_B*tau_B
  if (peakA <= peakB) {
    n1 <- n_A
    n2 <- n_B
    tau1 <- tau_A
    tau2 <- tau_B
  } else {
    n1 <- n_B
    n2 <- n_A
    tau1 <- tau_B
    tau2 <- tau_A
  }
  n_A <- n1
  n_B <- n2
  tau_A <- tau1
  tau_B <- tau2
  
  beta <- 1 - alpha
  
  # initialize range of characteristic function
  omega_max <- np*pi/t_run
  domega <- 2*omega_max/np
  omega <- seq(from = -omega_max, to = omega_max, by = domega)
  
  # calculate characteristic function of peaks and inverse Fourier transform
  CF_A <- exp(n_A*(1/(1 - 1i*tau_A*omega) - 1))
  CF_B <- exp(n_B*(1/(1 - 1i*tau_B*omega) - 1))
  f_A <- rev(Mod(ifft(CF_A)))
  f_B <- rev(Mod(ifft(CF_B)))
  t <- (1:(np + 1))*(pi/omega_max) # ensures time domain is correct
  f <- alpha*f_A + beta*f_B
  
  # fetch peak parameters
  tA.i <- which.max(f_A)
  tB.i <- which.max(f_B)
  tA_B.i <- sort(c(tA.i:tB.i))
  
  # calculate probability distributions of two states
  x <- seq(from = 1, to = 0, length = length(tA_B.i))
  P_AA <- sqrt(a*b*(1 - x)/x)*
    exp(-a*(1 - x) - b*x)*
    besselI(sqrt(4*a*b*x*(1 - x)), 1)
  P_AB <- a*
    exp(-a*(1 - x) - b*x)*
    besselI(sqrt(4*a*b*x*(1 - x)), 0)
  P_A0 <- exp(-a)
  P_BB <- sqrt(a*b*x/(1 - x))*
    exp(-a*(1 - x) - b*x)*
    besselI(sqrt(4*a*b*x*(1 - x)), 1)
  P_BA <- b*
    exp(-a*(1 - x) - b*x)*
    besselI(sqrt(4*a*b*x*(1 - x)), 0)
  P_B0 <- exp(-b)
  P_AA_norm <- (1 - P_A0)*P_AA/sum(P_AA + P_AB, na.rm = T)
  P_AB_norm <- (1 - P_A0)*P_AB/sum(P_AA + P_AB, na.rm = T)
  P_BB_norm <- (1 - P_B0)*P_BB/sum(P_BB + P_BA, na.rm = T)
  P_BA_norm <- (1 - P_B0)*P_BA/sum(P_BB + P_BA, na.rm = T)
  P_A <- alpha*P_AA_norm + beta*P_BA_norm
  P_A[1] <- alpha*P_A0
  P_A <- c(rep(0, tA.i - 1), P_A, rep(0, np - tB.i + 1))
  P_B <- beta*P_BB_norm + alpha*P_AB_norm
  P_B[length(P_B)] <- beta*P_B0
  P_B <- c(rep(0, tA.i - 1), P_B, rep(0, np - tB.i + 1))
  P_A[is.nan(P_A)] <- 0
  P_B[is.nan(P_B)] <- 0
  
  # convolution of peaks with state probabilities
  f_A_conv <- convolve(f_A, rev(P_A), type = "o")[tA.i:(tA.i + np)]
  f_B_conv <- convolve(f_B, rev(P_B), type = "o")[tB.i:(tB.i + np)]
  f_conv <- f_A_conv + f_B_conv
  
  # return result
  return(list(t = t, f_A = f_A, f_B = f_B, f = f,
              f_A_conv = f_A_conv, f_B_conv = f_B_conv, f_conv = f_conv,
              P_A = P_A, P_B = P_B))
}

preproc.Batman1 <- function(data, summary_df, i) {
  # Description
  # Helper function to preprocess chromatography data for batch.eval.stoch().
  # Used for data sets that were successfully evaluated with unified equation.
  
  # Arguments
  # data        Chromatography data to be processed.
  # summary_df  The summary data frame generated by batch.eval.kue()
  # i           The index of data (in summary_df) being processed.
  
  # Value
  # list        List of processed chromatography data together with characteristic
  #             function parameters calculated.
  
  # load required packages
  require(pracma)
  require(dplyr)
  require(minpack.lm)
  require(gamlss.dist)
  
  # baseline correction
  colnames(data) <- c("t", "f")
  all_peaks <- findpeaks(data$f, threshold = 0.3)
  first_peak_start <- min(min(all_peaks[, 3]), 0.02*nrow(data))
  if (first_peak_start == 1) {message("Peak found at first time point; baseline might be unreliable.")}
  last_peak_ends <- max(max(all_peaks[, 4]), 0.98*nrow(data))
  if (last_peak_ends == nrow(data)) {message("Peak found at last time point; baseline might be unreliable.")}
  baseline_sample <- data[c(1:first_peak_start, last_peak_ends:nrow(data)), ]
  baseline_sd <- sd(baseline_sample$f)
  baseline <- data[abs(data$f - mean(baseline_sample$f)) < 3*baseline_sd, ]
  fit <- lm(f ~ poly(t, 3, raw = TRUE), data = baseline)
  baseline_pred <- predict(fit, newdata = data, type = "response")
  data$f <- data$f - baseline_pred
  
  # fetch peak parameters and calculate characteristic function parameters
  t_M <- summary_df$t_M[i]
  s_M <- summary_df$s_M[i]
  t_A <- summary_df$t_A[i]
  t_B <- summary_df$t_B[i]
  tR_A <- t_A - t_M
  tR_B <- t_B - t_M
  w_A <- summary_df$w_A[i]
  w_B <- summary_df$w_B[i]
  sd_A <- w_A/sqrt(8*log(2))
  sd_B <- w_B/sqrt(8*log(2))
  var_A <- (sd_A)^2
  var_B <- (sd_B)^2
  tau_A <- (var_A/tR_A)/2
  tau_B <- (var_B/tR_B)/2
  n_A <- tR_A/tau_A
  n_B <- tR_B/tau_B
  
  # deconvolute chromatography data from dead time signal
  if (is.na(s_M)) {
    data$t <- data$t - t_M
    data <- data %>% filter(t > 0)
  } else {
    if (s_M < 1e-3) {s_M <- 0.01} # sanity check
    tM_data <- subset(data, t > (t_M - 3*s_M) & t < (t_M + 3*s_M))
    tM.LM <- function(parms) {
      pred <- parms[4]*dexGAUS(tM_data$t, parms[1], exp(parms[2]), exp(parms[3]))
      obs <- tM_data$f
      res <- (pred - obs)
      return(res)
    }
    dx <- diff(tM_data$t)
    norm <- sum(0.5*(tM_data$f[-1] + tM_data$f[-length(tM_data$f)])*dx)
    exp_val <- sum(0.5*(tM_data$t[-1]*tM_data$f[-1] +
                          tM_data$t[-length(tM_data$t)]*tM_data$f[-length(tM_data$f)])*dx)/norm
    var <- sum(0.5*(((tM_data$t[-1] - exp_val)^2*tM_data$f[-1]) +
                      ((tM_data$t[-length(tM_data$t)] - exp_val)^2*tM_data$f[-length(tM_data$f)]))*dx)/norm
    m3 <- sum(0.5*(((tM_data$t[-1] - exp_val)^3*tM_data$f[-1]) +
                     ((tM_data$t[-length(tM_data$t)] - exp_val)^3*tM_data$f[-length(tM_data$f)]))*dx)/norm
    skew <- ifelse(m3 > 0, (m3/2)^(1/3), 1e-6)
    sd <- ifelse(var > 0, sqrt(var), 1e-6)
    parms <- c(exp_val, log(sd), log(skew), norm)
    fitval <- nls.lm(par = parms, fn = tM.LM)
    fit_parms <- fitval$par
    tM_pred <- fit_parms[4]*dexGAUS(data$t, fit_parms[1], exp(fit_parms[2]), exp(fit_parms[3]))
    deconv_A <- data$f - tM_pred
    deconv_t <- data$t - t_M
    data <- data.frame(t = deconv_t, f = deconv_A)
    data <- data %>% filter(t > 0)
  }
  
  # normalization
  data$f <- data$f/sum(data$f)
  
  # fetch run time parameters
  np <- nrow(data)
  t_run <- max(data$t)
  
  return(list(t = data$t, f = data$f, np = np, t_run = t_run,
              n_A = n_A, n_B = n_B, tau_A = tau_A, tau_B = tau_B))
}

preproc.Batman2 <- function(data, summary_df, i, threshold = 0.3, minSNR = 10) {
  # Description
  # Helper function to preprocess chromatography data for batch.eval.stoch().
  # Used for data sets that could not be evaluated with unified equation.
  
  # Arguments
  # data        Chromatography data to be processed.
  # summary_df  The summary data frame generated by batch.eval.kue()
  # i           The index of data (in summary_df) being processed.
  # All other arguments have the same defaults and are passed to fit.batch.kue().
  
  # Value
  # list        List of processed chromatography data together with characteristic
  #             function parameters extrapolated.
  
  # load required packages
  require(pracma)
  require(dplyr)
  require(minpack.lm)
  require(gamlss.dist)
  
  # grouping of data in summary_df with the same metadata
  group_data <- summary_df %>%
    filter(comp_name == summary_df$comp_name[i],
           col_name == summary_df$col_name[i],
           Temp == summary_df$Temp[i])
  if (sum(!is.na(group_data$t_A)) > 0) {
    group_data <- group_data %>% mutate(
      tR_A = t_A - t_M,
      tR_B = t_B - t_M,
      sd_A = w_A/sqrt(8*log(2)),
      sd_B = w_B/sqrt(8*log(2)),
      var_A = (sd_A)^2,
      var_B = (sd_B)^2,
      tau_A = (var_A/tR_A)/2,
      tau_B = (var_B/tR_B)/2,
      n_A = tR_A/tau_A,
      n_B = tR_B/tau_B
    )
    model_n_A <- lm(n ~ flow, data = data.frame(n = group_data$n_A, flow = group_data$flow))
    n_A <- predict(model_n_A, newdata = data.frame(flow = summary_df$flow[i]))
    if (sum(!is.na(group_data$t_A)) == 1) {
      n_A <- summary_df$flow[i]*na.omit(group_data$n_A)[[1]]/na.omit(group_data$flow)[[1]]
    }
    model_n_B <- lm(n ~ flow, data = data.frame(n = group_data$n_B, flow = group_data$flow))
    n_B <- predict(model_n_B, newdata = data.frame(flow = summary_df$flow[i]))
    if (sum(!is.na(group_data$t_A)) == 1) {
      n_B <- summary_df$flow[i]*na.omit(group_data$n_B)[[1]]/na.omit(group_data$flow)[[1]]
    }
    tau_A <- mean(group_data$tau_A, na.rm = TRUE)
    tau_B <- mean(group_data$tau_B, na.rm = TRUE)
  } else {
    lower_temps <- summary_df %>%
      filter(comp_name == summary_df$comp_name[i],
             col_name == summary_df$col_name[i],
             Temp < summary_df$Temp[i],
             !is.na(t_A))
    if (nrow(lower_temps) == 0) {stop("No data available with different flow rates, and no lower temperature data available for extrapolation.")}
    nearest_lower_temp <- max(lower_temps$Temp, na.rm = TRUE)
    group_data <- summary_df %>%
      filter(comp_name == summary_df$comp_name[i],
             col_name == summary_df$col_name[i],
             Temp == nearest_lower_temp)
    calc_vals <- group_data %>%
      mutate(
        tau_A = (t_A - t_M)/2*(w_A/sqrt(8*log(2)))^2,
        tau_B = (t_B - t_M)/2*(w_B/sqrt(8*log(2)))^2,
        n_A = (t_A - t_M)/tau_A,
        n_B = (t_B - t_M)/tau_B
      ) %>%
      select(tau_A, tau_B, n_A, n_B)
    mean_vals <- calc_vals %>%
      summarize(
        mean_n_A = mean(n_A, na.rm = TRUE),
        mean_n_B = mean(n_B, na.rm = TRUE),
        mean_tau_A = mean(tau_A, na.rm = TRUE),
        mean_tau_B = mean(tau_B, na.rm = TRUE)
      )
    n_A <- mean_vals$mean_n_A
    n_B <- mean_vals$mean_n_B
    tau_A <- mean_vals$mean_tau_A
    tau_B <- mean_vals$mean_tau_B
  }
  
  # baseline correction
  colnames(data) <- c("t", "f")
  all_peaks <- findpeaks(data$f, threshold = 0.3)
  first_peak_start <- min(min(all_peaks[, 3]), 0.02*nrow(data))
  if (first_peak_start == 1) {message("Peak found at first time point; baseline might be unreliable.")}
  last_peak_ends <- max(max(all_peaks[, 4]), 0.98*nrow(data))
  if (last_peak_ends == nrow(data)) {message("Peak found at last time point; baseline might be unreliable.")}
  baseline_sample <- data[c(1:first_peak_start, last_peak_ends:nrow(data)), ]
  baseline_sd <- sd(baseline_sample$f)
  baseline <- data[abs(data$f - mean(baseline_sample$f)) < 3*baseline_sd, ]
  fit <- lm(f ~ poly(t, 3, raw = TRUE), data = baseline)
  baseline_pred <- predict(fit, newdata = data, type = "response")
  data$f <- data$f - baseline_pred
  
  # finding dead time
  t_M <- summary_df$t_M[i]
  s_M <- summary_df$s_M[i]
  peaks <- findpeaks(data$f,
                     threshold = threshold,
                     minpeakheight = minSNR*baseline_sd,
                     minpeakdistance = (1/100)*nrow(data))
  peaks <- data.frame(peaks)
  colnames(peaks) <- c("f", "ti", "ti_start", "ti_end")
  peaks <- peaks %>% arrange(ti)
  peaks$t <- data$t[peaks$ti]
  if (is.na(t_M)) {
    peak_tab <- tryCatch(
      read.csv("peak_table.csv", header = TRUE, sep = ",", dec = "."),
      error = function(e) NULL
    )
    if (!is.null(peak_tab)) {
      peak_values <- peak_tab %>% filter(file_name == summary_df$file[i])
      if (!is.na(peak_values$t_M)) {
        t0peak <- peaks %>% slice_min(abs(peaks$t - peak_values$t_M), n = 1, with_ties = FALSE)
        t_M <- t0peak$t
        t_Mi <- t0peak$ti
        h_M <- t0peak$f
        cutoff_idx <- max(10, round(0.5*nrow(data)/max(data$t))) # 0.5 min with security minimum
        idx_left <- (t_Mi - cutoff_idx):t_Mi
        idx_right <- t_Mi:(t_Mi + cutoff_idx)
        w_M.i <- which.min(abs(data$f[idx_left] - h_M/2)) + t_Mi - cutoff_idx - 1
        w_M.i2 <- which.min(abs(data$f[idx_right] - h_M/2)) + t_Mi - 1
        w_M <- data$t[w_M.i2] - data$t[w_M.i]
        s_M <- w_M/sqrt(8*log(2))
      }
    } else {
      t0peak <- peaks[1, ]
      t_M <- t0peak$t
      t_Mi <- t0peak$ti
      h_M <- t0peak$f
      cutoff_idx <- max(10, round(0.5*nrow(data)/max(data$t))) # 0.5 min with security minimum
      idx_left <- (t_Mi - cutoff_idx):t_Mi
      idx_right <- t_Mi:(t_Mi + cutoff_idx)
      w_M.i <- which.min(abs(data$f[idx_left] - h_M/2)) + t_Mi - cutoff_idx - 1
      w_M.i2 <- which.min(abs(data$f[idx_right] - h_M/2)) + t_Mi - 1
      w_M <- data$t[w_M.i2] - data$t[w_M.i]
      s_M <- w_M/sqrt(8*log(2))
    }
  }
  
  # sanity check n and tau start values
  if (!exists("tau_A") || is.na(tau_A) || tau_A < 1e-6) {tau_A <- 1e-5}
  if (!exists("tau_B") || is.na(tau_B) || tau_B < 1e-6) {tau_B <- 1e-5}
  if (!exists("n_A") || is.na(n_A) || n_A == 0) {n_A <- 1e2}
  if (!exists("n_B") || is.na(n_B) || n_B == 0) {n_B <- 1e2}
  colaesced_peak <- filter(peaks, t > t_M)$t[which.max(filter(peaks, t > t_M)$f)]
  if (abs(n_A*tau_A - colaesced_peak + t_M) > 0.5) {
    n_A <- (colaesced_peak - t_M)/tau_A
  }
  if (abs(n_B*tau_B - colaesced_peak + t_M) > 0.5) {
    n_B <- (colaesced_peak - t_M)/tau_B
  }
  
  # deconvolute chromatography data from dead time signal
  if (is.na(s_M)) {
    data$t <- data$t - t_M
    data <- data %>% filter(t > 0)
  } else {
    if (s_M < 1e-3) {s_M <- 0.01} # sanity check
    tM_data <- subset(data, t > (t_M - 3*s_M) & t < (t_M + 3*s_M))
    tM.LM <- function(parms) {
      pred <- parms[4]*dexGAUS(tM_data$t, parms[1], exp(parms[2]), exp(parms[3]))
      obs <- tM_data$f
      res <- (pred - obs)
      return(res)
    }
    dx <- diff(tM_data$t)
    norm <- sum(0.5*(tM_data$f[-1] + tM_data$f[-length(tM_data$f)])*dx)
    exp_val <- sum(0.5*(tM_data$t[-1]*tM_data$f[-1] +
                          tM_data$t[-length(tM_data$t)]*tM_data$f[-length(tM_data$f)])*dx)/norm
    var <- sum(0.5*(((tM_data$t[-1] - exp_val)^2*tM_data$f[-1]) +
                        ((tM_data$t[-length(tM_data$t)] - exp_val)^2*tM_data$f[-length(tM_data$f)]))*dx)/norm
    m3 <- sum(0.5*(((tM_data$t[-1] - exp_val)^3*tM_data$f[-1]) +
                        ((tM_data$t[-length(tM_data$t)] - exp_val)^3*tM_data$f[-length(tM_data$f)]))*dx)/norm
    skew <- ifelse(m3 > 0, (m3/2)^(1/3), 1e-6)
    sd <- ifelse(var > 0, sqrt(var), 1e-6)
    parms <- c(exp_val, log(sd), log(skew), norm)
    fitval <- nls.lm(par = parms, fn = tM.LM)
    fit_parms <- fitval$par
    tM_pred <- fit_parms[4]*dexGAUS(data$t, fit_parms[1], exp(fit_parms[2]), exp(fit_parms[3]))
    deconv_A <- data$f - tM_pred
    deconv_t <- data$t - t_M
    data <- data.frame(t = deconv_t, f = deconv_A)
    data <- data %>% filter(t > 0)
  }
  
  # normalization
  data$f <- data$f/sum(data$f)
  
  # fetch run time parameters
  np <- nrow(data)
  t_run <- max(data$t)
  
  return(list(t = data$t, f = data$f, np = np, t_run = t_run, t_M = t_M,
              n_A = n_A, n_B = n_B, tau_A = tau_A, tau_B = tau_B))
}

batch.eval.stoch <- function(path, alpha = 0.5, threshold = 0.3, minSNR = 10, max_conv = 10) {
  # Description
  # This function performs batch evaluation of the stochastic model on chromatography
  # files provided in .csv format if unified equation modeling has been performed.
  
  # Arguments
  # path      The path where the files are found.
  # alpha     Fraction of the first eluted isomer in the injected sample. Default is 0.5.
  # max_conv  The maximum number of interconversions during the run, passed on as a starting parameter to DE algorithm. Defaults to 10.
  # All other arguments have the same defaults and are passed to fit.batch.kue().
  
  # Value
  # The function call generates a new folder in the path and saves the chromatography
  # plots and fitted stochastic models as .png files. The summary_log.txt and summary_data.csv
  # are appended with fitted characteristic parameters, interconversion rates and
  # fitted thermodynamic parameters.
  
  # References
  # Felinger, A. (2008) Molecular dynamic theories in chromatography. Journal of
  # Chromatography A, 1184, 20-41.
  #
  # Sepsey, A., Németh, D. R., Németh, G., Felinger, A. (2018) Rate constant determination of 
  # interconverting enantiomers by chiral chromatography using a stochastic model.
  # Journal of Chromatography A, 1564, 155-162.
  
  # load required packages
  require(data.table)
  require(readr)
  require(fs)
  require(ggplot2)
  require(stringr)
  require(tools)
  require(dplyr)
  require(broom)
  require(DEoptim)
  require(minpack.lm)
  require(parallel)
  require(parallelly)
  
  # set path as working directory and create directory for plots
  setwd(path)
  output_dir <- "plots_stoch"
  dir_create(output_dir)
  
  # read summary and log files from kue modeling, update files
  summary_df <- tryCatch({
    fread("summary_data.csv")
  }, error = function(e) {
    stop("The file 'summary_data.csv' was not found. Perform unified equation modeling first.")
  })
  summary_df <- summary_df[, lapply(.SD, function(x) {
    if (is.logical(x)) as.numeric(x) else x
  })]
  summary_df <- summary_df %>%
    arrange(Temp, desc(flow))
  log_df <- fread("summary_log.txt", na.strings = "")
  summary_df$n_A <- as.numeric(NA)
  summary_df$n_B <- as.numeric(NA)
  summary_df$tau_A <- as.numeric(NA)
  summary_df$tau_B <- as.numeric(NA)
  summary_df$a <- as.numeric(NA)
  summary_df$b <- as.numeric(NA)
  log_df$stoch_status <- "Success"
  log_df$stoch_error <- rep(NA_character_, nrow(log_df))
  
  # evaluate files with stochastic model
  for (i in 1:nrow(summary_df)) {
    message(sprintf("Evaluating file %d out of %d", i, nrow(summary_df)))
    
    # extract file and name
    file <- summary_df$file[i]
    base_name <- tools::file_path_sans_ext(basename(file))
    base_name_parts <- str_split(base_name, "_", simplify = TRUE)[1:4]
    comp_name <- base_name_parts[1]
    col_name <- base_name_parts[2]
    flow <- as.numeric(str_replace(base_name_parts[3], "^(\\d)(\\d+)$", "\\1.\\2"))
    Temp <- as.numeric(base_name_parts[4]) + 273.15
    title <- paste(paste(comp_name, col_name, flow, Temp - 273.15, sep = ", "), "(comp, col, flow ml/min, °C)")
    
    # read chromatograghy files
    encoding <- readr::guess_encoding(file)$encoding
    sep <- smart_fread(file, encoding = encoding)$sep
    dec <- smart_fread(file, encoding = encoding)$dec
    df <- read.csv(file, header = FALSE, sep = sep, dec = dec, fileEncoding = encoding)
    if (!is.numeric(df[, 1])) {
      df <- read.csv(file, header = TRUE, sep = sep, dec = dec, fileEncoding = encoding)
    }
    colnames(df) <- c("V1", "V2")
    
    tryCatch({
      # preprocess chromatography data
      if (!is.na(summary_df$t_A[i])) {
        preproc_data <- preproc.Batman1(df, summary_df, i)
        overwrite_t <- FALSE
      } else {
        preproc_data <- preproc.Batman2(df, summary_df, i, threshold, minSNR)
        summary_df$t_M[i] <- preproc_data$t_M
        overwrite_t <- TRUE
      }
      
      #required optim functions
      Batman.DEoptim <- function(parms) {
        # Differential Evolution Optimization genetic algorithm
        parms1 <- parms[1]
        parms2 <- parms[2]
        parms3 <- parms[3]
        parms4 <- parms[4]
        parms5 <- parms[5]
        sum((Batman(preproc_data$np, preproc_data$t_run,
                    parms1, parms2, parms3,
                    parms4, parms5, parms5,
                    alpha)$f_conv -
               preproc_data$f)^2, na.rm = TRUE)
      }
      Batman.LM <- function(parms) {
        # Levenberg-Marquadt algorithm
        pred <- Batman(preproc_data$np, preproc_data$t_run,
                       parms[1], parms[2], parms[3],
                       parms[4], parms[5], parms[6],
                       alpha)$f_conv
        obs <- preproc_data$f
        res <- (pred - obs)
        return(res)
      }
      
      # set bounds and control for genetic algorithm
      lower <- c(preproc_data$n_A/2, preproc_data$n_B/2,
                 preproc_data$tau_A/2, preproc_data$tau_B/2, 0)
      upper <- c(2*preproc_data$n_A, 2*preproc_data$n_B,
                 2*preproc_data$tau_A, 2*preproc_data$tau_B, max_conv)
      for (start_i in 1:2) {
        if (lower[start_i] > upper[start_i]) {
          lower[start_i] <- 0
          upper[start_i] <- 1e4
        }
      }
      for (start_i in 3:4) {
        if (lower[start_i] > upper[start_i]) {
          lower[start_i] <- 0
          upper[start_i] <- 1e-1
        }
      }
      
      DEoptim_control <- DEoptim.control(itermax = 200,
                                         trace = TRUE,
                                         parallelType = 1,
                                         parVar = c("alpha", "Batman"))
      DEoptim_res <- DEoptim(Batman.DEoptim, lower, upper, control = DEoptim_control)
      
      # Levenberg-Marquadt algorithm with warm start
      parms <- c(DEoptim_res$optim$bestmem, DEoptim_res$optim$bestmem[5])
      names(parms) <- NULL
      fitval <- nls.lm(par = parms, fn = Batman.LM, control = nls.lm.control(maxiter = 100))
      fit_parms <- fitval$par
      pred <- Batman(preproc_data$np, preproc_data$t_run,
                     fit_parms[1], fit_parms[2], fit_parms[3],
                     fit_parms[4], fit_parms[5], fit_parms[6],
                     alpha)
      
      # write fitted parameters
      summary_df$n_A[i] <- fit_parms[1]
      summary_df$n_B[i] <- fit_parms[2]
      summary_df$tau_A[i] <- fit_parms[3]
      summary_df$tau_B[i] <- fit_parms[4]
      summary_df$a[i] <- fit_parms[5]
      summary_df$b[i] <- fit_parms[6]
      if (overwrite_t) {
        summary_df$t_A[i] <- min(summary_df$n_A[i]*summary_df$tau_A[i],
                                 summary_df$n_B[i]*summary_df$tau_B[i]) + summary_df$t_M[i]
        summary_df$t_B[i] <- max(summary_df$n_A[i]*summary_df$tau_A[i],
                                 summary_df$n_B[i]*summary_df$tau_B[i]) + summary_df$t_M[i]
        summary_df$w_A[i] <- ifelse(summary_df$t_A[i] < summary_df$t_B[i],
                                    sqrt(16*log(2)*summary_df$n_A[i]*summary_df$tau_A[i]^2),
                                    sqrt(16*log(2)*summary_df$n_B[i]*summary_df$tau_B[i]^2))
        summary_df$w_B[i] <- ifelse(summary_df$t_A[i] > summary_df$t_B[i],
                                    sqrt(16*log(2)*summary_df$n_A[i]*summary_df$tau_A[i]^2),
                                    sqrt(16*log(2)*summary_df$n_B[i]*summary_df$tau_B[i]^2))
      }
      
      #export plots
      plot <- ggplot() +
        geom_line(aes(x = preproc_data$t, y = preproc_data$f)) +
        geom_line(aes(x = pred$t, y = pred$f_conv), col = "red", linetype = 2) +
        labs(title = title,
             subtitle = "Stochastic modeling",
             x = "Time (min)",
             y = "Intensity") +
        geom_segment(aes(x = 0, y = max(pred$f_conv),
                         xend = max(pred$t)/20, yend = max(pred$f_conv)),
                     col = "red", linetype = 2) +
        annotate("text", x = max(pred$t)/15, y = max(pred$f_conv), label = "fit") +
        theme_bw()
      out_file <- file.path(output_dir, paste0(base_name, ".png"))
      ggsave(out_file, plot = plot, width = 6, height = 4, dpi = 300)
    }, error = function(e) {
      log_df$stoch_status[log_df$file == file] <<- "Error"
      log_df$stoch_error[log_df$file == file] <<- conditionMessage(e)
    })
  }
  
  #add linreg a, b ~ t_A, t_B
  summary_df <- summary_df %>%
    group_by(comp_name, col_name, Temp) %>%
    group_modify(~ {
      valid_data <- .x %>% filter(!is.na(a) & !is.na(b) & !is.na(t_A) & !is.na(t_B))
      f_intercept <- NA_real_
      f_slope <- NA_real_
      f_intercept_se <- NA_real_
      f_slope_se <- NA_real_
      f_r2 <- NA_real_
      r_intercept <- NA_real_
      r_slope <- NA_real_
      r_intercept_se <- NA_real_
      r_slope_se <- NA_real_
      r_r2 <- NA_real_
      EP_intercept <- NA_real_
      EP_slope <- NA_real_
      EP_intercept_se <- NA_real_
      EP_slope_se <- NA_real_
      EP_r2 <- NA_real_
      if (nrow(valid_data) >= 2) {
        f_model <- lm(a ~ I(60*t_A), data = valid_data)
        f_tidy_model <- tidy(f_model)
        f_glance_model <- glance(f_model)
        f_intercept <- f_tidy_model$estimate[f_tidy_model$term == "(Intercept)"]
        f_slope <- f_tidy_model$estimate[f_tidy_model$term == "I(60 * t_A)"]
        f_intercept_se <- f_tidy_model$std.error[f_tidy_model$term == "(Intercept)"]
        f_slope_se <- f_tidy_model$std.error[f_tidy_model$term == "I(60 * t_A)"]
        f_r2 <- f_glance_model$r.squared
        r_model <- lm(b ~ I(60*t_B), data = valid_data)
        r_tidy_model <- tidy(r_model)
        r_glance_model <- glance(r_model)
        r_intercept <- r_tidy_model$estimate[r_tidy_model$term == "(Intercept)"]
        r_slope <- r_tidy_model$estimate[r_tidy_model$term == "I(60 * t_B)"]
        r_intercept_se <- r_tidy_model$std.error[r_tidy_model$term == "(Intercept)"]
        r_slope_se <- r_tidy_model$std.error[r_tidy_model$term == "I(60 * t_B)"]
        r_r2 <- r_glance_model$r.squared
      }
      .x %>%
        mutate(kstoch_f_tA = f_slope,
               kstoch_f_tA_se = f_slope_se,
               kstoch_f_tA_intcpt = f_intercept,
               kstoch_f_tA_intcpt_se = f_intercept_se,
               kstoch_f_tA_r2 = f_r2,
               kstoch_r_tB = r_slope,
               kstoch_r_tB_se = r_slope_se,
               kstoch_r_tB_intcpt = r_intercept,
               kstoch_r_tB_intcpt_se = r_intercept_se,
               kstoch_r_tB_r2 = r_r2)
    }) %>%
    ungroup()
  
  #add linreg Eyring-Polanyi
  summary_df <- summary_df %>%
    group_by(comp_name, col_name) %>%
    group_modify(~ {
      valid_data <- .x %>% filter(!is.na(kstoch_f_tA) & !is.na(Temp)) %>%
        group_by(Temp) %>%
        summarise(mean_kstoch_f = mean(kstoch_f_tA, na.rm = TRUE),
                  mean_kstoch_r = mean(kstoch_r_tB, na.rm = TRUE),
                  .groups = "drop")
      EP_f_intercept <- NA_real_
      EP_f_slope <- NA_real_
      EP_f_intercept_se <- NA_real_
      EP_f_slope_se <- NA_real_
      EP_f_r2 <- NA_real_
      EP_r_intercept <- NA_real_
      EP_r_slope <- NA_real_
      EP_r_intercept_se <- NA_real_
      EP_r_slope_se <- NA_real_
      EP_r_r2 <- NA_real_
      if (nrow(valid_data) >= 2) {
        EP_f_model <- lm(I(log(mean_kstoch_f/Temp)) ~ I(1/Temp), data = valid_data)
        EP_f_tidy_model <- tidy(EP_f_model)
        EP_f_glance_model <- glance(EP_f_model)
        EP_f_intercept <- EP_f_tidy_model$estimate[EP_f_tidy_model$term == "(Intercept)"]
        EP_f_slope <- EP_f_tidy_model$estimate[EP_f_tidy_model$term == "I(1/Temp)"]
        EP_f_intercept_se <- EP_f_tidy_model$std.error[EP_f_tidy_model$term == "(Intercept)"]
        EP_f_slope_se <- EP_f_tidy_model$std.error[EP_f_tidy_model$term == "I(1/Temp)"]
        EP_f_r2 <- EP_f_glance_model$r.squared
        EP_r_model <- lm(I(log(mean_kstoch_r/Temp)) ~ I(1/Temp), data = valid_data)
        EP_r_tidy_model <- tidy(EP_r_model)
        EP_r_glance_model <- glance(EP_r_model)
        EP_r_intercept <- EP_r_tidy_model$estimate[EP_r_tidy_model$term == "(Intercept)"]
        EP_r_slope <- EP_r_tidy_model$estimate[EP_r_tidy_model$term == "I(1/Temp)"]
        EP_r_intercept_se <- EP_r_tidy_model$std.error[EP_r_tidy_model$term == "(Intercept)"]
        EP_r_slope_se <- EP_r_tidy_model$std.error[EP_r_tidy_model$term == "I(1/Temp)"]
        EP_r_r2 <- EP_r_glance_model$r.squared
      }
      .x %>%
        mutate(dH_f = -EP_f_slope*8.314,
               dH_f_se = EP_f_slope_se*8.314,
               dS_f = (EP_f_intercept - log(2.084e10))*8.314,
               dS_f_se = EP_f_intercept_se*8.314,
               EP_f_r2 = EP_f_r2,
               dH_r = -EP_r_slope*8.314,
               dH_r_se = EP_r_slope_se*8.314,
               dS_r = (EP_r_intercept - log(2.084e10))*8.314,
               dS_r_se = EP_r_intercept_se*8.314,
               EP_r_r2 = EP_r_r2)
    }) %>%
    ungroup()
  
  # rewrite files
  fwrite(log_df, "summary_log.txt", sep = "\t", na = "")
  fwrite(summary_df, "summary_data.csv", sep = ",", dec = ".")
  message("Processing complete.")
}




