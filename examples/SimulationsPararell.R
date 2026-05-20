# ---------- setup ----------
library(tidyverse)
library(future.apply)
library(ALDEx2)            # asegúrate de tener todos los paquetes instalados
library(coda4microbiome)   # idem
library(codabiocom)

setwd("~/2026/Irene/CodaPaper")
dir.create("LRRelev_sim_results", showWarnings = FALSE)
dir.create("LRRelev_metrics", showWarnings = FALSE)

# ---------- escenarios (idéntico a tu grid) ----------
G_vals <- c(2,3,5,10)
n_vals <- c(10,50,100,1000)
m_vals <- c(10,50,100,1000)
k_vals <- c(2,4,8,16,32)
replicates <- 50

scenarios <- expand_grid(
  m = m_vals,
  n = n_vals,
  k = k_vals,
  G = G_vals,
  rep = 1:replicates
) %>%
  filter(k <= m) %>%
  arrange(m, n, k, G, rep) %>%
  mutate(sim_id = row_number())

# ---------- define run_one_sim (adapté tu bloque) ----------
run_one_sim <- function(scenario_row) {
  # scenario_row is a single-row tibble/data.frame
  m <- scenario_row$m
  n <- scenario_row$n
  k <- scenario_row$k
  G <- scenario_row$G
  rep <- scenario_row$rep
  sim_id <- scenario_row$sim_id



  generate_fold_changes <- function(k, G, prop_all_groups = 0.5,
                                    fc_down = c(0.1,0.5), fc_up = c(2,10)) {
    fold_change_list <- vector("list", k)

    for(otu_idx in 1:k){
      if(runif(1) < prop_all_groups){
        # cambia en todos los grupos
        fold_change_list[[otu_idx]] <- sapply(1:G, function(g) {
          if(runif(1) < 0.5) runif(1, fc_down[1], fc_down[2])
          else runif(1, fc_up[1], fc_up[2])
        })
      } else {
        # cambia solo en un grupo
        g_change <- sample(1:G, 1)
        fold_vec <- rep(1, G)   # no cambia en otros grupos
        fold_vec[g_change] <- if(runif(1) < 0.5) runif(1, fc_down[1], fc_down[2])
        else runif(1, fc_up[1], fc_up[2])
        fold_change_list[[otu_idx]] <- fold_vec
      }
    }

    return(fold_change_list)
  }

  simulate_OTU_data_counts <- function(n, m, G, k_relevant,
                                       N_total = runif(1,10^4, 10^8),
                                       prop_all_groups = 0.5,
                                       fc_down = c(0.3,0.5),
                                       fc_up = c(2,5)) {

    # OTUs relevantes
    rel_otus <- 1:k_relevant

    # generar fold-changes heterogéneos
    fold_change_list <- generate_fold_changes(
      k = k_relevant,
      G = G,
      prop_all_groups = prop_all_groups,
      fc_down = fc_down,
      fc_up = fc_up
    )

    # asignar grupos
    group <- sample(1:G, n, replace = TRUE)

    # matriz de conteos
    X_counts <- matrix(0, nrow = n, ncol = m)
    colnames(X_counts) <- paste0("OTU", 1:m)

    # generar cada muestra
    for(i in 1:n){

      # alpha base NUEVA por muestra
      alpha_i <- rep(1, m)

      # aplicar fold-change a los OTUs relevantes
      for(otu_idx in 1:k_relevant){
        fc <- fold_change_list[[otu_idx]][ group[i] ]
        alpha_i[ rel_otus[otu_idx] ] <- alpha_i[ rel_otus[otu_idx] ] * fc
      }

      # convertir a probabilidades
      probs_i <- alpha_i / sum(alpha_i)

      # simular recuentos
      X_counts[i, ] <- rmultinom(1, size = N_total, prob = probs_i)
    }

    # salida final
    data.frame(
      SampleID = paste0("S", 1:n),
      Group = factor(group),
      X_counts
    )
  }


  file_name <- file.path("LRRelev_sim_results",
                         paste0("sim_", sim_id, "_G", G, "_n", n,
                                "_m", m, "_k", k, "_rep", rep, ".rds"))

  metrics_file <- file.path("LRRelev_metrics",
                            paste0("metrics_", sim_id, "_G", G, "_n", n,
                                   "_m", m, "_k", k, "_rep", rep, ".rds"))

  # skip if metrics already exist
  if (file.exists(metrics_file)) {
    message("Skipping sim ", sim_id, " (metrics exist)")
    return(NULL)
  }

  # --- simulate (or reuse existing sim file if present) ---
  if (file.exists(file_name)) {
    sim_data <- readRDS(file_name)
  } else {
    sim_data <- simulate_OTU_data_counts(n = n, m = m, G = G, k_relevant = k)
    saveRDS(sim_data, file_name)
  }

  OTU_names <- paste0("OTU", 1:m)
  relevant_otus <- OTU_names[1:k]

  # ---- LRRelev ----
  t0 <- Sys.time()
  res <- tryCatch({
    LRRelev(data = sim_data[, c(-1, -2)],
            sample = sim_data$SampleID,
            group  = sim_data$Group,
            taxa   = OTU_names,
            otus   = OTU_names,
            threshold = 2,
            X = NULL,
            method = "hanley")
  }, error = function(e) {
    message("Error in LRRelev sim ", sim_id, ": ", e$message)
    return(NULL)
  })
  t1 <- Sys.time()
  t1Mio <- as.numeric(difftime(t1, t0, units = "secs"))

  # ---- ALDEx2 ----
  t0 <- Sys.time()
  resALDEx <- tryCatch({
    ALDEx2::aldex(
      reads = t(sim_data[, -c(1,2)]),
      conditions = as.character(sim_data$Group),
      denom = "all",
      mc.samples = 128,
      test = "kw",
      effect = TRUE
    )
  }, error = function(e) {
    message("Error in ALDEx sim ", sim_id, ": ", e$message)
    NULL
  })
  t1 <- Sys.time()
  t1Aldex <- as.numeric(difftime(t1, t0, units = "secs"))

  metricsALDEx <- NULL
  if(!is.null(resALDEx)) {
    detected_otusALDEx <- rownames(resALDEx)[resALDEx$kw.eBH < 0.05]
    TPALDEx <- sum(detected_otusALDEx %in% relevant_otus)
    FPALDEx <- sum(!detected_otusALDEx %in% relevant_otus)
    FNALDEx <- k - TPALDEx
    TPRALDEx <- TPALDEx / k
    FDRALDEx <- ifelse(TPALDEx + FPALDEx > 0, FPALDEx / (TPALDEx + FPALDEx), NA)

    metricsALDEx <- list(
      TP = TPALDEx, FP = FPALDEx, FN = FNALDEx,
      TPR = TPRALDEx, FDR = FDRALDEx,
      n_detected = length(detected_otusALDEx),
      relevant_otus = relevant_otus,
      detected_otus = detected_otusALDEx,
      sim_k = k, sim_m = m, sim_n = n, sim_G = G,
      Time = t1Aldex
    )
  }

  # ---- coda_glmnet (MALU) only if G == 2 ----
  metricsMALU <- NULL
  if (G == 2) {
    t0 <- Sys.time()
    resMALU <- tryCatch({
      coda4microbiome::coda_glmnet(sim_data[, c(-1,-2)], sim_data$Group, showPlots = FALSE)
    }, error = function(e) {
      message("Error in MALU sim ", sim_id, ": ", e$message)
      return(NULL)
    })
    t1 <- Sys.time()
    t1MALU <- as.numeric(difftime(t1, t0, units = "secs"))

    if (!is.null(resMALU)) {
      detected_otusMALU <- resMALU$taxa.name
      TPMALU <- sum(detected_otusMALU %in% relevant_otus)
      FPMALU <- sum(!detected_otusMALU %in% relevant_otus)
      FNMALU <- k - TPMALU
      TPRMALU <- TPMALU / k
      FDRMALU <- ifelse(TPMALU + FPMALU > 0, FPMALU / (TPMALU + FPMALU), NA)
      k_selected <- length(resMALU$taxa.num)
      auc_relev <- resMALU$`mean cv-AUC`

      metricsMALU <- list(
        TP = TPMALU, FP = FPMALU, FN = FNMALU, TPR = TPRMALU, FDR = FDRMALU,
        n_detected = length(detected_otusMALU),
        relevant_otus = relevant_otus,
        detected_otus = detected_otusMALU,
        k_selected = k_selected,
        sim_k = k, sim_m = m, sim_n = n, sim_G = G,
        mean_auc_relev = auc_relev,
        Time = t1MALU
      )
    }
  }

  # ---- Collect LRRelev metrics (if res exists) ----
  metricsLR <- NULL
  if (!is.null(res)) {
    detected_otus <- res$OTUSRelev
    OTUs_relev <- 1:k
    TP <- sum(detected_otus %in% relevant_otus)
    FP <- sum(!detected_otus %in% relevant_otus)
    FN <- k - TP
    TPR <- TP / k
    FDR <- ifelse(TP + FP > 0, FP / (TP + FP), NA)
    k_selected <- length(res$OTUSRelev)
    AUC_matrix <- res$AUCs
    if (!is.null(AUC_matrix)) {
      auc_relev <- AUC_matrix[OTUs_relev, OTUs_relev, drop = FALSE]
      auc_nonrelev <- AUC_matrix[-OTUs_relev, -OTUs_relev, drop = FALSE]
      mean_auc_relev <- mean(auc_relev[upper.tri(auc_relev)])
      mean_auc_nonrelev <- if(length(auc_nonrelev) > 1) mean(auc_nonrelev[upper.tri(auc_nonrelev)]) else NA
    } else {
      mean_auc_relev <- NA
      mean_auc_nonrelev <- NA
    }

    metricsLR <- list(
      TP = TP, FP = FP, FN = FN, TPR = TPR, FDR = FDR,
      n_detected = length(detected_otus),
      relevant_otus = relevant_otus,
      detected_otus = detected_otus,
      k_selected = k_selected,
      sim_k = k, sim_m = m, sim_n = n, sim_G = G,
      mean_auc_relev = mean_auc_relev,
      mean_auc_nonrelev = mean_auc_nonrelev,
      Time = t1Mio
    )
  }

  # ---- Save combined metrics ----
  saveRDS(list(metrics = metricsLR, Aldex = metricsALDEx, MALU = metricsMALU),
          file = metrics_file)

  message("Finished sim ", sim_id)
  return(metrics_file)  # return filename as small confirmation
}

# ---------- Parallel execution ----------
# Choose number of workers (leave 1-2 cores free)
workers <- max(1, parallel::detectCores() - 1)
plan(multisession, workers = workers)

# Run in parallel, reproducible seeds
res_files <- future_lapply(
  seq_len(nrow(scenarios)),
  function(i) {
    run_one_sim(scenarios[i, ])
  },
  future.seed = TRUE
)

# Optional: inspect which metrics files were produced
res_files <- compact(res_files)
length(res_files)
