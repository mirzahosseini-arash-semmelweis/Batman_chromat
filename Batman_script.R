#TODO
#check calculation of kue_r
#check if linreg kue vs t_A is necessary
#simulate se
#rescale A and B AUC (?)

#fit function for batch analysis
fit.batch.kue <- function(data, t0 = "run", threshold = 0.3, minSNR = 50, minpeakdist = 5,
                          enantio = TRUE, microrev = TRUE, A0 = 0.5) {
  require(pracma)
  require(dplyr)
  require(stringr)
  #error handling
  if (missing(data)) stop("No data provided.")
  if (!is.numeric(threshold)) stop("Argument 'threshold' must be numeric.")
  if (!is.numeric(minSNR)) stop("Argument 'minSNR' must be numeric.")
  if (!is.numeric(minpeakdist)) stop("Argument 'minpeakdist' must be numeric.")
  if (!is.logical(enantio)) stop("Argument 'enantio' must be boolean")
  if (!is.logical(microrev)) stop("Argument 'microrev' must be boolean")
  if (!is.numeric(A0)) stop("Argument 'A0' must be numeric.")
  if (A0 > 1 | A0 < 0) stop("Argument 'A0' must be between 0 and 1.")
  #peak find
  data <- data.frame(data[, 1:2])
  colnames(data) <- c("t", "A")
  baseline_sample <- data[1:(0.02*nrow(data)), ]
  baseline_sd <- sd(baseline_sample$A)
  baseline <- data[abs(data$A) < 12*baseline_sd, ]
  base_corr_coef <- summary(lm(A ~ t, baseline))$coef
  base_corr <- data$t*base_corr_coef[2, 1] + base_corr_coef[1, 1]
  data$A <- data$A - base_corr
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
  t0peak <- peaks[1, ]
  peaks <- peaks[-1, ]
  if (is.numeric(t0)) {
    if (abs(data$t[t0peak$ti] - t0) > 1) {
      warning("Warning: Identified dead time farther than 1 min from provided t0; using provided value. \n")
      peaks <- rbind(t0peak, peaks)
      t_M <- t0
    } else {
      t_M <- data$t[t0peak$ti]
    }
  }
  if (t0 == "run") {
    t_M <- data$t[t0peak$ti]
  } else {stop("Possible values for t0 are numeric or 'run'.")}
  npeaks <- nrow(peaks)
  if (npeaks == 0) {
    stop("Something went wrong during peak picking. No peaks found other than dead time.")
  }
  if (npeaks == 1) {
    stop("Only one peak found; perform stochastic modeling.")
  }
  if (npeaks > 2) {
    peaks <- peaks %>% filter(row_number() %in% c(1, n()))
    warning("More than 2 peaks found; peaks at extrema are kept. \n")
  }
  npeaks <- nrow(peaks)
  if (npeaks == 2) {
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
    w_A.i <- which.min(abs(data$A[peaks[1, 3]:t_Ai] - h_A/2)) + peaks[1, 3] - 1
    w_A.i2 <- which.min(abs(data$A[t_Ai:peaks[1, 4]] - h_A/2)) + t_Ai - 1
    w_A <- ifelse(h_mid >= h_A/2, 2*abs(data$t[w_A.i] - t_A), data$t[w_A.i2] - data$t[w_A.i])
    w_As <- w_A*60
    w_A.1 <- data$t[w_A.i]
    w_A.2 <- ifelse(h_mid >= h_A/2, w_A + w_A.1, data$t[w_A.i2])
    w_B.i <- which.min(abs(data$A[peaks[2, 3]:t_Bi] - h_B/2)) + peaks[2, 3] - 1
    w_B.i2 <- which.min(abs(data$A[t_Bi:peaks[2, 4]] - h_B/2)) + t_Bi - 1
    w_B <- ifelse(h_mid >= h_B/2, 2*abs(data$t[w_B.i2] - t_B), data$t[w_B.i2] - data$t[w_B.i])
    w_Bs <- w_B*60
    w_B.1 <- data$t[w_B.i]
    w_B.2 <- ifelse(h_mid >= h_B/2, w_B + w_B.1, data$t[w_B.i2])
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
    N <- mean(5.54*(c(t_A, t_B)/c(w_A, w_B))^2)
  } else {stop("Something went wrong during peak picking.")}
  if (enantio == TRUE) {K <- 1}
  if (microrev == TRUE) {
    t_xs <- t_As
  } else {t_xs <- t_Bs}
  #calculate kue
  calc.kue1 <- function(k, K, t_As, t_Bs, t_xs, s_As, s_Bs, N, h_p, A0) {
    -(1/t_As)*(
      log(
        (1 - A0)*exp(-K*k*t_xs)
        *
          (
            (100*exp(-(t_Bs - t_As)^2/(8*s_Bs^2)) - h_p*exp(-(t_Bs - t_As)^2/(2*s_Bs^2)))
            /
              (s_Bs*sqrt(2*pi))
            -
              100/(t_Bs - t_As)
          )
        +
          (100*(1 - A0) + A0*(100 - h_p*(1 + sqrt(2/(pi*N)))))/(t_Bs - t_As)
      )
      -
        log(
          A0*(
            (h_p - 100*exp(-(t_As - t_Bs)^2/(8*s_As^2)))/(s_As*sqrt(2*pi))
            +
              (100 - h_p*(1 + sqrt(2/(pi*N))))/(t_Bs - t_As)
          )
        )
    )
  }
  calc.kue2 <- function(k, K, t_As, t_Bs, t_xs, s_As, s_Bs, N, h_p, A0) {
    -(1/t_As)*(
      log(
        (1 - A0)*exp(-K*k*t_xs)
        *
          (
            (100*exp(-(t_As - t_Bs)^2/(8*s_Bs^2)) - h_p)
            /
              (s_Bs*sqrt(2*pi))
            +
              (h_p*(1 - sqrt(2/(pi*N))) - 100)/(t_Bs - t_As)
          )
        +
          (100*A0 + (1 - A0)*(100 - h_p*(1 - sqrt(2/(pi*N)))))/(t_Bs - t_As)
      )
      -
        log(
          A0*(
            (h_p*exp(-(t_Bs - t_As)^2/(2*s_As^2)) - 100*exp(-(t_Bs - t_As)^2/(8*s_As^2)))/(s_As*sqrt(2*pi))
            +
              100/(t_Bs - t_As)
          )
        )
    )
  }
  k <- seq(from = 1e-6, to = 1e-2, length = 100)
  if (h_A >= h_B) {
    f_k <- calc.kue1(k, K, t_As, t_Bs, t_xs, s_As, s_Bs, N, h_p, A0)
  } else {
    f_k <- calc.kue2(k, K, t_As, t_Bs, t_xs, s_As, s_Bs, N, h_p, A0)
  }
  root.i <- which.min(abs(k - f_k))
  root <- k[root.i]
  k <- seq(from = root/2, to = root*2, length = 100)
  if (h_A >= h_B) {
    f_k <- calc.kue1(k, K, t_As, t_Bs, t_xs, s_As, s_Bs, N, h_p, A0)
  } else {
    f_k <- calc.kue2(k, K, t_As, t_Bs, t_xs, s_As, s_Bs, N, h_p, A0)
  }
  lm <- summary(lm((k - f_k) ~ k))
  kue_f <- as.numeric(-lm$coefficients[1, 1]/lm$coefficients[2, 1])
  if (t_xs == t_As) {
    kue_r <- kue_f*K*t_As/t_Bs
  } else {kue_r <- kue_f*K}
  #export results
  result <- list(t_M = t_M,
                 t_Ai = t_Ai,
                 t_Bi = t_Bi,
                 t_A = t_A,
                 t_B = t_B,
                 t_As = t_As,
                 t_Bs = t_Bs,
                 t_xs = t_xs,
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
                 kue_r = kue_r)
  return(result)
}

#function required for reading files
smart_fread <- function(file, encoding) {
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

batch.eval.kue <- function(path, t0 = "run", threshold = 0.3, minSNR = 50, minpeakdist = 5,
                           enantio = TRUE, microrev = TRUE, A0 = 0.5) {
  require(fs)
  require(readr)
  require(data.table)
  require(tools)
  require(ggplot2)
  require(stringr)
  require(dplyr)
  require(broom)
  
  #reading files
  setwd(path)
  csv_files <- list.files(pattern = "(?i)\\.csv$", full.names = TRUE)
  csv_files <- csv_files[!tolower(basename(csv_files)) %in% "summary_data.csv"]
  if (length(csv_files) == 0) {stop("No .csv files found in path directory.")}
  output_dir <- "plots_kue"
  dir_create(output_dir)
  log_entries <- list()
  summary_rows <- list()
  for (file in csv_files) {
    base_name <- tools::file_path_sans_ext(basename(file))
    base_name_parts <- str_split(base_name, "_", simplify = TRUE)
    comp_name <- base_name_parts[1]
    col_name <- base_name_parts[2]
    flow <- as.numeric(str_replace(base_name_parts[3], "^(\\d)(\\d+)$", "\\1.\\2"))
    Temp <- str_extract(base_name, "\\d+$") |> as.numeric() + 273.15
    title <- str_replace_all(base_name, "_", ", ")
    title <- str_replace(title, "(\\b)(\\d)(\\d+)(,\\s*\\d+\\s*$)", "\\2.\\3\\4")
    title <- paste(title, "(comp, col, flow ml/min, °C)")
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
      #table
      result <- fit.batch.kue(data = df, t0, threshold, minSNR, minpeakdist,
                              enantio, microrev, A0)
      Gue_f <- -log(result$kue_f/(1.380662e-23*Temp/6.626176e-34))*8.31441*Temp/1000
      Gue_r <- -log(result$kue_r/(1.380662e-23*Temp/6.626176e-34))*8.31441*Temp/1000
      result$Gue_f <- Gue_f
      result$Gue_r <- Gue_r
      table <- data.frame(result)
      table <- table %>% select(t_M,
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
      #plot
      h_max <- max(result$h_A, result$h_B)
      t_max <- result$t_B
      plot <- ggplot(data = df) +
        geom_line(aes(x = V1, y = V2)) +
        labs(title = title, x = "time [min]", y = "intensity") +
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
                 label = round(result$h_mid, 3), col = "dodgerblue3") +
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
                              t_M = NA,
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
  #create log file and table
  log_df <- rbindlist(log_entries, fill = TRUE)
  summary_df <- rbindlist(summary_rows, fill = TRUE)
  #add linreg kue ~ t_A
  summary_df <- summary_df %>%
    group_by(comp_name, col_name, Temp) %>%
    group_modify(~ {
      valid_data <- .x %>% filter(!is.na(kue_f) & !is.na(kue_r) & !is.na(t_A))
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
      if (nrow(valid_data) >= 2) {
        f_model <- lm(kue_f ~ I(60*t_A), data = valid_data)
        f_tidy_model <- tidy(f_model)
        f_glance_model <- glance(f_model)
        f_intercept <- f_tidy_model$estimate[f_tidy_model$term == "(Intercept)"]
        f_slope <- f_tidy_model$estimate[f_tidy_model$term == "I(60 * t_A)"]
        f_intercept_se <- f_tidy_model$std.error[f_tidy_model$term == "(Intercept)"]
        f_slope_se <- f_tidy_model$std.error[f_tidy_model$term == "I(60 * t_A)"]
        f_r2 <- f_glance_model$r.squared
        r_model <- lm(kue_r ~ I(60*t_A), data = valid_data)
        r_tidy_model <- tidy(r_model)
        r_glance_model <- glance(r_model)
        r_intercept <- r_tidy_model$estimate[r_tidy_model$term == "(Intercept)"]
        r_slope <- r_tidy_model$estimate[r_tidy_model$term == "I(60 * t_A)"]
        r_intercept_se <- r_tidy_model$std.error[r_tidy_model$term == "(Intercept)"]
        r_slope_se <- r_tidy_model$std.error[r_tidy_model$term == "I(60 * t_A)"]
        r_r2 <- r_glance_model$r.squared
      }
      .x %>%
        mutate(kue_f_tA = f_slope,
               kue_f_tA_se = f_slope_se,
               kue_f_tA_intcpt = f_intercept,
               kue_f_tA_intcpt_se = f_intercept_se,
               kue_f_tA_r2 = f_r2,
               kue_r_tA = r_slope,
               kue_r_tA_se = r_slope_se,
               kue_r_tA_intcpt = r_intercept,
               kue_r_tA_intcpt_se = r_intercept_se,
               kue_r_tA_r2 = r_r2)
    }) %>%
    ungroup()
  fwrite(log_df, "summary_log.txt", sep = "\t", na = "")
  fwrite(summary_df, "summary_data.csv", sep = ",", dec = ".")
  message("Processing complete.")
}

batch.update.kue <- function(path, drop = "") {
  require(data.table)
  summary_df <- fread("summary_data.csv")
  log_df <- fread("summary_log.txt", na.strings = "")
  #drop kue data of rows in which "file" is in drop = c("filename")
  files_to_drop <- drop
  cols_to_drop <- c("t_M", "t_A", "t_B", "h_A", "h_B", "h_p", "w_A", "w_B",
                    "A", "B", "K", "N", "kue_f", "kue_r", "Gue_f", "Gue_r",
                    "kue_f_tA", "kue_f_tA_se", "kue_f_tA_intcpt", "kue_f_tA_intcpt_se",
                    "kue_f_tA_r2", "kue_r_tA", "kue_r_tA_se", "kue_r_tA_intcpt",
                    "kue_r_tA_intcpt_se", "kue_r_tA_r2")
  for (file_drop in files_to_drop) {
    summary_df[summary_df$file == file_drop, c(cols_to_drop)] <- NA
    log_df[log_df$file == file_drop, "kue_status"] <- "Dropped"
    log_df[log_df$file == file_drop, "kue_error"] <- NA
  }
  fwrite(log_df, "summary_log.txt", sep = "\t", na = "")
  fwrite(summary_df, "summary_data.csv", sep = ",", dec = ".")
  message("Update complete.")
}

Batman <- function(np, t_run, n_A, n_B, tau_A, tau_B, a, b, alpha) {
  require(pracma)
  np <- np - 1 #ensures number of points is np
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
  omega_max <- np*pi/t_run
  domega <- 2*omega_max/np
  omega <- seq(from = -omega_max, to = omega_max, by = domega)
  CF_A <- exp(n_A*(1/(1 - 1i*tau_A*omega) - 1))
  CF_B <- exp(n_B*(1/(1 - 1i*tau_B*omega) - 1))
  f_A <- rev(Mod(ifft(CF_A)))
  f_B <- rev(Mod(ifft(CF_B)))
  t <- (1:(np + 1))*(pi/omega_max)
  f <- alpha*f_A + beta*f_B
  tA.i <- which.max(f_A)
  tB.i <- which.max(f_B)
  tA_B.i <- sort(c(tA.i:tB.i))
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
  P_A <- alpha*(1 - P_A0)*(P_AA + P_BA)/sum(P_AA + P_BA, na.rm = T)
  P_A[1] <- alpha*P_A0
  P_A <- c(rep(0, tA.i - 1), P_A, rep(0, np - tB.i + 1))
  P_B <- beta*(1 - P_B0)*(P_BB + P_AB)/sum(P_BB + P_AB, na.rm = T)
  P_B[length(P_B)] <- beta*P_B0
  P_B <- c(rep(0, tA.i - 1), P_B, rep(0, np - tB.i + 1))
  P_A[is.nan(P_A)] <- 0
  P_B[is.nan(P_B)] <- 0
  f_A_conv <- convolve(f_A, rev(P_A), type = "o")[tA.i:(tA.i + np)]
  f_B_conv <- convolve(f_B, rev(P_B), type = "o")[tB.i:(tB.i + np)]
  f_conv <- f_A_conv + f_B_conv
  return(list(t = t, f_A = f_A, f_B = f_B, f = f,
              f_A_conv = f_A_conv, f_B_conv = f_B_conv, f_conv = f_conv,
              P_A = P_A, P_B = P_B))
}

preproc.Batman1 <- function(data, summary_df, i) {
  colnames(data) <- c("t", "f")
  baseline_sample <- data[1:(0.02*nrow(data)), ]
  baseline_sd <- sd(baseline_sample$f)
  baseline <- data[abs(data$f) < 12*baseline_sd, ]
  base_corr_coef <- summary(lm(f ~ t, baseline))$coef
  base_corr <- data$t*base_corr_coef[2, 1] + base_corr_coef[1, 1]
  data$f <- data$f - base_corr
  
  t_M <- summary_df$t_M[i]
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
  
  data$t <- data$t - t_M
  data <- data %>% filter(t >= 0)
  lower_bound <- 1
  upper_bound <- which.min(abs(data$t - max(0.01*max(data$t), tR_A - 20*sd_A)))
  data$f[lower_bound:upper_bound] <- rnorm(upper_bound, 0, baseline_sd)
  data$f <- data$f/sum(data$f)
  np <- nrow(data)
  t_run <- max(data$t)
  
  return(list(t = data$t, f = data$f, np = np, t_run = t_run,
              n_A = n_A, n_B = n_B, tau_A = tau_A, tau_B = tau_B))
}

preproc.Batman2 <- function(data, summary_df, i) {
  group_data <- summary_df %>%
    filter(comp_name == summary_df$comp_name[i],
           col_name == summary_df$col_name[i],
           Temp == summary_df$Temp[i])
  if(nrow(group_data) < 2) {stop("No data available with different flow rates to extrapolate n and tau.")}
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
  model_n_B <- lm(n ~ flow, data = data.frame(n = group_data$n_B, flow = group_data$flow))
  n_B <- predict(model_n_B, newdata = data.frame(flow = summary_df$flow[i]))
  tau_A <- mean(group_data$tau_A, na.rm = TRUE)
  tau_B <- mean(group_data$tau_B, na.rm = TRUE)
  
  colnames(data) <- c("t", "f")
  baseline_sample <- data[1:(0.02*nrow(data)), ]
  baseline_sd <- sd(baseline_sample$f)
  baseline <- data[abs(data$f) < 12*baseline_sd, ]
  base_corr_coef <- summary(lm(f ~ t, baseline))$coef
  base_corr <- data$t*base_corr_coef[2, 1] + base_corr_coef[1, 1]
  data$f <- data$f - base_corr
  
  peaks <- findpeaks(data$f,
                     threshold = 0.3,
                     minpeakheight = 50*baseline_sd,
                     minpeakdistance = (5/100)*nrow(data))
  peaks <- data.frame(peaks)
  colnames(peaks) <- c("f", "ti", "ti_start", "ti_end")
  peaks <- peaks %>% arrange(ti)
  t0peak <- peaks[1, ]
  t_M <- data$t[t0peak$ti]
  
  data$t <- data$t - t_M
  data <- data %>% filter(t >= 0)
  lower_bound <- 1
  upper_bound <- which.min(abs(data$t - max(0.01*max(data$t), n_A*tau_A - 40*n_A*tau_A^2)))
  data$f[lower_bound:upper_bound] <- rnorm(upper_bound, 0, baseline_sd)
  data$f <- data$f/sum(data$f)
  np <- nrow(data)
  t_run <- max(data$t)
  
  return(list(t = data$t, f = data$f, np = np, t_run = t_run, t_M = t_M,
              n_A = n_A, n_B = n_B, tau_A = tau_A, tau_B = tau_B))
}

batch.eval.stoch <- function(path, alpha = 0.5) {
  require(data.table)
  require(readr)
  require(ggplot2)
  require(stringr)
  require(tools)
  require(dplyr)
  require(DEoptim)
  require(minpack.lm)
  require(parallel)
  
  setwd(path)
  output_dir <- "plots_stoch"
  dir_create(output_dir)
  
  #read in summary and log files from kue modeling
  summary_df <- tryCatch({
    fread("summary_data.csv")
  }, error = function(e) {
    stop("The file 'summary_data.csv' was not found. Perform unified equation modeling first.")
  })
  log_df <- fread("summary_log.txt", na.strings = "")
  summary_df$n_A <- as.numeric(NA)
  summary_df$n_B <- as.numeric(NA)
  summary_df$tau_A <- as.numeric(NA)
  summary_df$tau_B <- as.numeric(NA)
  summary_df$a <- as.numeric(NA)
  summary_df$b <- as.numeric(NA)
  log_df$stoch_status <- "Success"
  log_df$stoch_error <- rep(NA_character_, nrow(log_df))
  
  #loop for stoch modeling
  for (i in 1:nrow(summary_df)) {
    file <- summary_df$file[i]
    base_name <- tools::file_path_sans_ext(basename(file))
    title <- str_replace_all(base_name, "_", ", ")
    title <- str_replace(title, "(\\b)(\\d)(\\d+)(,\\s*\\d+\\s*$)", "\\2.\\3\\4")
    title <- paste(title, "(comp, col, flow ml/min, °C)")
    
    encoding <- readr::guess_encoding(file)$encoding
    sep <- smart_fread(file, encoding = encoding)$sep
    dec <- smart_fread(file, encoding = encoding)$dec
    df <- read.csv(file, header = FALSE, sep = sep, dec = dec, fileEncoding = encoding)
    if (!is.numeric(df[, 1])) {
      df <- read.csv(file, header = TRUE, sep = sep, dec = dec, fileEncoding = encoding)
    }
    
    tryCatch({
      if (!is.na(summary_df$t_A[i])) {
        preproc_data <- preproc.Batman1(df, summary_df, i)
      } else {
        preproc_data <- preproc.Batman2(df, summary_df, i)
        summary_df$t_M[i] <- preproc_data$t_M
      }
      
      # np <- preproc_data$np
      # t_run <- preproc_data$t_run
      # alpha <- alpha
      
      #required optim functions
      Batman.DEoptim <- function(parms) {
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
        pred <- Batman(preproc_data$np, preproc_data$t_run,
                       parms[1], parms[2], parms[3],
                       parms[4], parms[5], parms[6],
                       alpha)$f_conv
        obs <- preproc_data$f
        res <- (pred - obs)
        return(res)
      }
      
      lower <- c(preproc_data$n_A/2, preproc_data$n_B/2,
                 preproc_data$tau_A/2, preproc_data$tau_B/2, 0)
      upper <- c(2*preproc_data$n_A, 2*preproc_data$n_B,
                 2*preproc_data$tau_A, 2*preproc_data$tau_B, 6)
      
      DEoptim_control <- DEoptim.control(itermax = 200,
                                         trace = TRUE,
                                         parallelType = 1,
                                         parVar = c("alpha", "Batman"))
      DEoptim_res <- DEoptim(Batman.DEoptim, lower, upper, control = DEoptim_control)
      
      parms <- c(DEoptim_res$optim$bestmem, DEoptim_res$optim$bestmem[5])
      fitval <- nls.lm(par = parms, fn = Batman.LM, control = nls.lm.control(maxiter = 100))
      fit_parms <- fitval$par
      pred <- Batman(preproc_data$np, preproc_data$t_run,
                     fit_parms[1], fit_parms[2], fit_parms[3],
                     fit_parms[4], fit_parms[5], fit_parms[6],
                     alpha)
      
      # write fit_parms
      summary_df$n_A[i] <- fit_parms[1]
      summary_df$n_B[i] <- fit_parms[2]
      summary_df$tau_A[i] <- fit_parms[3]
      summary_df$tau_B[i] <- fit_parms[4]
      summary_df$a[i] <- fit_parms[5]
      summary_df$b[i] <- fit_parms[6]
      
      #export plots
      plot <- ggplot() +
        geom_line(aes(x = preproc_data$t, y = preproc_data$f)) +
        geom_line(aes(x = pred$t, y = pred$f_conv), col = "red", linetype = 2) +
        labs(title = title,
             subtitle = "Stochastic modeling",
             x = "time [min]",
             y = "intensity") +
        geom_segment(aes(x = 0, y = max(pred$f_conv),
                         xend = max(pred$t)/20, yend = max(pred$f_conv)),
                     col = "red", linetype = 2) +
        annotate("text", x = max(pred$t)/15, y = max(pred$f_conv), label = "fit") +
        theme_bw()
      out_file <- file.path(output_dir, paste0(base_name, ".png"))
      ggsave(out_file, plot = plot, width = 6, height = 4, dpi = 300)
    }, error = function(e) {
      log_df$stoch_status[i] <<- "Error"
      log_df$stoch_error[i] <<- conditionMessage(e)
    })
  }
  
  #add linreg a/b ~ t_A
  summary_df <- summary_df %>%
    group_by(comp_name, col_name, Temp) %>%
    group_modify(~ {
      valid_data <- .x %>% filter(!is.na(a) & !is.na(b) & !is.na(t_M))
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
        f_model <- lm(a ~ I(60*(n_A*tau_A + t_M)), data = valid_data)
        f_tidy_model <- tidy(f_model)
        f_glance_model <- glance(f_model)
        f_intercept <- f_tidy_model$estimate[f_tidy_model$term == "(Intercept)"]
        f_slope <- f_tidy_model$estimate[f_tidy_model$term == "I(60 * (n_A * tau_A + t_M))"]
        f_intercept_se <- f_tidy_model$std.error[f_tidy_model$term == "(Intercept)"]
        f_slope_se <- f_tidy_model$std.error[f_tidy_model$term == "I(60 * (n_A * tau_A + t_M))"]
        f_r2 <- f_glance_model$r.squared
        r_model <- lm(b ~ I(60*(n_A*tau_A + t_M)), data = valid_data)
        r_tidy_model <- tidy(r_model)
        r_glance_model <- glance(r_model)
        r_intercept <- r_tidy_model$estimate[r_tidy_model$term == "(Intercept)"]
        r_slope <- r_tidy_model$estimate[r_tidy_model$term == "I(60 * (n_A * tau_A + t_M))"]
        r_intercept_se <- r_tidy_model$std.error[r_tidy_model$term == "(Intercept)"]
        r_slope_se <- r_tidy_model$std.error[r_tidy_model$term == "I(60 * (n_A * tau_A + t_M))"]
        r_r2 <- r_glance_model$r.squared
      }
      .x %>%
        mutate(kstoch_f_tA = f_slope,
               kstoch_f_tA_se = f_slope_se,
               kstoch_f_tA_intcpt = f_intercept,
               kstoch_f_tA_intcpt_se = f_intercept_se,
               kstoch_f_tA_r2 = f_r2,
               kstoch_r_tA = r_slope,
               kstoch_r_tA_se = r_slope_se,
               kstoch_r_tA_intcpt = r_intercept,
               kstoch_r_tA_intcpt_se = r_intercept_se,
               kstoch_r_tA_r2 = r_r2)
    }) %>%
    ungroup()
  
  #add linreg Eyring-Polanyi
  summary_df <- summary_df %>%
    group_by(comp_name, col_name) %>%
    group_modify(~ {
      valid_data <- .x %>% filter(!is.na(kstoch_f_tA) & !is.na(Temp)) %>%
        group_by(Temp) %>%
        summarise(mean_kstoch = mean(kstoch_f_tA, na.rm = TRUE), .groups = "drop")
      EP_intercept <- NA_real_
      EP_slope <- NA_real_
      EP_intercept_se <- NA_real_
      EP_slope_se <- NA_real_
      EP_r2 <- NA_real_
      if (nrow(valid_data) >= 2) {
        EP_model <- lm(I(log(mean_kstoch/Temp)) ~ I(1/Temp), data = valid_data)
        EP_tidy_model <- tidy(EP_model)
        EP_glance_model <- glance(EP_model)
        EP_intercept <- EP_tidy_model$estimate[EP_tidy_model$term == "(Intercept)"]
        EP_slope <- EP_tidy_model$estimate[EP_tidy_model$term == "I(1/Temp)"]
        EP_intercept_se <- EP_tidy_model$std.error[EP_tidy_model$term == "(Intercept)"]
        EP_slope_se <- EP_tidy_model$std.error[EP_tidy_model$term == "I(1/Temp)"]
        EP_r2 <- EP_glance_model$r.squared
      }
      .x %>%
        mutate(dH = -EP_slope*8.314,
               dH_se = EP_slope_se*8.314,
               dS = (EP_intercept - log(2.084e10))*8.314,
               dS_de = EP_intercept_se*8.314,
               EP_r2 = EP_r2)
    }) %>%
    ungroup()
  
  fwrite(log_df, "summary_log.txt", sep = "\t", na = "")
  fwrite(summary_df, "summary_data.csv", sep = ",", dec = ".")
  message("Processing complete.")
}

################################################################################
#test batch import

path <- "~/01_Research/03_Statistics/Batman/sample data/lorabatmancell2"
batch.eval.kue(path = path)
batch.eval.stoch(path = path)
