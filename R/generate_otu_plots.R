#' Generate comprehensive OTU plots from selection output
#'
#' This function creates multiple types of plots for relevant OTUs based on the results of an OTU selection procedure.
#' It can generate abundance plots, CLR-transformed boxplots, relative proportion distributions, incremental AUC curves,
#' and PCA plots colored by original groups, k-means clusters, or hierarchical clustering.
#'
#' @details
#' The function accepts an `output` list containing the selected OTUs (`OTUSRelev`), the associated data (`dataImp`),
#' and the original OTU information (`OTUS`). Plots are generated without loading the packages explicitly, using
#' fully qualified calls (`ggplot2::`, `dplyr::`, `tidyr::`, `compositions::`, `factoextra::`).
#'
#' Available plot types:
#' \describe{
#'   \item{"abundance"}{Boxplot of imputed abundances of relevant OTUs (log10 scale).}
#'   \item{"clr"}{Boxplot of CLR-transformed abundances for relevant OTUs.}
#'   \item{"proportion_violin"}{Violin plots of relative proportions by group, with mean (red) and median (blue) points.}
#'   \item{"proportion_density"}{Density plots of relative proportions by group, faceted by OTU.}
#'   \item{"proportion_bar"}{Stacked bar plots of mean OTU proportions by group.}
#'   \item{"auc_incremental"}{Line plot of incremental AUC when adding OTUs, with confidence intervals and labels.}
#'   \item{"pca_group"}{PCA plot colored by original groups with ellipses.}
#'   \item{"pca_kmeans"}{PCA plot colored by k-means clusters with ellipses.}
#'   \item{"pca_hclust"}{PCA plot colored by hierarchical clustering (Ward.D2) with ellipses.}
#' }
#'
#' @param output List containing OTU selection results. Must include:
#'   \itemize{
#'     \item `dataImp`: Imputed OTU abundance table (samples x OTUs)
#'     \item `OTUS`: Table with OTU information and associations
#'     \item `OTUSRelev`: Vector of relevant OTU names
#'   }
#' @param group Vector of group labels for each sample (same length as `nrow(dataImp)`).
#' @param plot_types Character vector specifying which plots to generate. If `NULL`, all plots are returned.
#' @param otus_order Optional vector to order OTUs on axes; defaults to `OTUSRelev` order.
#'
#' @return A named list of ggplot objects corresponding to the requested plot types.
#'
#' @examples
#' \dontrun{
#' # Example with HIV dataset:
#' data(HIV)
#' Xnp <- model.matrix(y_HIV ~ MSM_HIV)
#' output1 <- LRRelev(data = x_HIV, sample = rownames(x_HIV), group = y_HIV,
#'   taxa = colnames(x_HIV), otus = colnames(x_HIV), cores = 2,
#'    X = Xnp, method = "hanley")
#' plots <- generate_otu_plots(output1, group=y_HIV)
#' plots$abundance      # Imputed abundance boxplot
#' plots$clr            # CLR boxplot
#' plots$proportion_violin # Violin plots with mean/median
#' plots$auc_incremental   # Incremental AUC plot
#' plots$pca_group         # PCA colored by group
#' }
#'
#' @references
#' - Wickham, H. (2016). *ggplot2: Elegant Graphics for Data Analysis*. Springer.
#' - Aitchison, J. (1986). *The Statistical Analysis of Compositional Data*. Chapman & Hall.
#' - Kassambara, A., & Mundt, F. (2020). *factoextra: Extract and Visualize the Results of Multivariate Data Analyses*. R package.
#'
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_point geom_line geom_ribbon geom_text
#'   geom_violin geom_bar geom_density facet_wrap scale_y_log10 scale_y_continuous
#'   scale_x_continuous theme_bw theme labs element_text arrow unit
#' @importFrom dplyr select mutate filter pull group_by summarise across
#' @importFrom tidyr pivot_longer
#' @importFrom compositions clr
#' @importFrom factoextra fviz_pca_ind
#' @import scales
#' @importFrom stats prcomp dist hclust cutree kmeans
#' @export
generate_otu_plots <- function(output, group, plot_types = NULL, otus_order = NULL) {

  dataImp <- output$dataImp
  dataImp$Grupo <- group
  relev <- output$OTUSRelev
  otus_info <- output$OTUS

  if(is.null(otus_order)) {
    otus_order <- dplyr::filter(otus_info, otus %in% relev) %>% dplyr::pull(otus)
  }

  plots <- list()

  # ---------------- Abundance ----------------
  if(is.null(plot_types) || "abundance" %in% plot_types) {
    df_plot <- dplyr::select(dataImp, dplyr::all_of(relev), Grupo) %>%
      tidyr::pivot_longer(-Grupo, names_to = "OTU", values_to = "Abundance") %>%
      dplyr::mutate(OTU = factor(OTU, levels = otus_order))

    plots$abundance <- ggplot2::ggplot(df_plot, ggplot2::aes(x = OTU, y = Abundance, fill = Grupo)) +
      ggplot2::geom_boxplot(ggplot2::aes(), outlier.alpha = 0.3) +
      ggplot2::scale_y_log10() +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
      ggplot2::labs(title = "Imputed abundances of relevant OTUs",
                    y = "Abundance (log10)",
                    x = "Relevant OTU")
  }

  # ---------------- CLR ----------------
  if(is.null(plot_types) || "clr" %in% plot_types) {
    X <- dplyr::select(dataImp, dplyr::all_of(relev))
    X_clr <- compositions::clr(X)
    X_clr <- as.data.frame(X_clr)
    colnames(X_clr) <- relev
    X_clr$Group <- dataImp$Grupo

    df_plot <- tidyr::pivot_longer(X_clr, -Group, names_to = "OTU", values_to = "CLR") %>%
      dplyr::mutate(OTU = factor(OTU, levels = otus_order))

    plots$clr <- ggplot2::ggplot(df_plot, ggplot2::aes(x = OTU, y = CLR, fill = Group)) +
      ggplot2::geom_boxplot(ggplot2::aes(), outlier.alpha = 0.3) +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
      ggplot2::labs(title = "CLR-transformed relevant OTUs",
                    y = "CLR abundance",
                    x = "Relevant OTU")
  }

  # ---------------- Proportion violin ----------------
  if(is.null(plot_types) || "proportion_violin" %in% plot_types) {
    X_prop <- dplyr::select(dataImp, dplyr::all_of(relev))
    X_prop <- X_prop / rowSums(X_prop)
    X_prop$Group <- dataImp$Grupo

    df_plot <- tidyr::pivot_longer(X_prop, -Group, names_to = "OTU", values_to = "Prop") %>%
      dplyr::mutate(OTU = factor(OTU, levels = otus_order))

    df_stats <- df_plot %>%
      dplyr::group_by(OTU, Group) %>%
      dplyr::summarise(
        MeanProp = mean(Prop),
        MedianProp = stats::median(Prop),
        .groups = "drop"
      )

    plots$proportion_violin <- ggplot2::ggplot(df_plot, ggplot2::aes(x = Group, y = Prop, fill = Group)) +
      ggplot2::geom_violin(alpha = 0.3, color = NA, scale = "width") +
      ggplot2::geom_point(data = df_stats, ggplot2::aes(x = Group, y = MeanProp), color = "red", size = 2, shape = 18) +
      ggplot2::geom_point(data = df_stats, ggplot2::aes(x = Group, y = MedianProp), color = "blue", size = 2, shape = 17) +
      ggplot2::facet_wrap(~OTU, scales = "free_y") +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                     strip.text = ggplot2::element_text(size = 10)) +
      ggplot2::labs(title = "Distribution of relative OTU proportions by group",
                    y = "Relative proportion",
                    x = "Group",
                    subtitle = "Violin: distribution, red = mean, blue = median")
  }

  # ---------------- Proportion density ----------------
  if(is.null(plot_types) || "proportion_density" %in% plot_types) {
    plots$proportion_density <- ggplot2::ggplot(df_plot, ggplot2::aes(x = Prop, fill = Group)) +
      ggplot2::geom_density(alpha = 0.3) +
      ggplot2::facet_wrap(. ~ OTU, scales = "free") +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
      ggplot2::labs(title = "Relative proportion of relevant OTUs",
                    y = "Density",
                    x = "Relative proportion",
                    subtitle = "Distribution by OTU and group")
  }

  # ---------------- Proportion bar ----------------
  if(is.null(plot_types) || "proportion_bar" %in% plot_types) {
    df_group <- X_prop %>%
      dplyr::group_by(Group) %>%
      dplyr::summarise(dplyr::across(dplyr::all_of(relev), mean)) %>%
      tidyr::pivot_longer(cols = dplyr::all_of(relev), names_to = "OTU", values_to = "MeanProp") %>%
      dplyr::mutate(OTU = factor(OTU, levels = otus_order))

    plots$proportion_bar <- ggplot2::ggplot(df_group, ggplot2::aes(x = Group, y = MeanProp, fill = OTU)) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
      ggplot2::theme_bw() +
      ggplot2::labs(title = "Mean proportion of relevant OTUs by group",
                    x = "Group",
                    y = "Proportion (%)",
                    fill = "OTU")
  }

  # ---------------- Incremental AUC ----------------
  if(is.null(plot_types) || "auc_incremental" %in% plot_types) {
    df_raw <- output$OTUS

    fila_par <- dplyr::filter(df_raw, otus == relev[2]) %>%
      dplyr::mutate(
        otus = paste0(relev[1], " +\n ", relev[2]),
        n_OTUs = 2,
        Added_OTU = paste0(relev[1], " +\n ", relev[2])
      )

    df_rest <- dplyr::filter(df_raw, otus %in% relev[c(-1,-2)]) %>%
      dplyr::mutate(
        n_OTUs = 3:(length(relev)),
        Added_OTU = otus
      )

    df_final <- rbind(fila_par, df_rest)

    plots$auc_incremental <- ggplot2::ggplot(df_final, ggplot2::aes(n_OTUs, assoc)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = ciLower, ymax = ciUpper), alpha = 0.2) +
      ggplot2::geom_line(size = 1.2) +
      ggplot2::geom_point(size = 3) +
      ggplot2::geom_text(ggplot2::aes(label = Added_OTU), vjust = -0.7, size = 3) +
      ggplot2::theme_bw() +
      ggplot2::scale_x_continuous(breaks = df_final$n_OTUs) +
      ggplot2::labs(title = "Incremental AUC when adding relevant OTUs",
                    x = "Number of OTUs included",
                    y = "AUC")
  }

  # ---------------- PCA ----------------
  if(any(is.null(plot_types) || plot_types %in% c("pca_group","pca_kmeans","pca_hclust"))) {
    X_clr <- as.data.frame(compositions::clr(dplyr::select(dataImp, dplyr::all_of(relev))))
    pca <- stats::prcomp(X_clr, center = TRUE, scale. = FALSE)

    df_plot <- data.frame(pca$x[,1:2]) %>%
      dplyr::mutate(
        Group = dataImp$Grupo,
        Kmeans = factor(stats::kmeans(X_clr, centers = length(unique(group)), nstart = 25)$cluster),
        Hclust = factor(stats::cutree(stats::hclust(stats::dist(X_clr), method = "ward.D2"), k = length(unique(group))))
      )

    if(is.null(plot_types) || "pca_group" %in% plot_types) {
      plots$pca_group <- factoextra::fviz_pca_ind(
        pca,
        geom.ind = "point",
        pointshape = 21,
        fill = df_plot$Group,
        col.ind = "black",
        alpha.ind = 0.7,
        addEllipses = TRUE
      )
    }

    if(is.null(plot_types) || "pca_kmeans" %in% plot_types) {
      plots$pca_kmeans <- factoextra::fviz_pca_ind(
        pca,
        geom.ind = "point",
        pointshape = 21,
        fill = df_plot$Kmeans,
        col.ind = "black",
        alpha.ind = 0.7,
        addEllipses = TRUE
      )
    }

    if(is.null(plot_types) || "pca_hclust" %in% plot_types) {
      plots$pca_hclust <- factoextra::fviz_pca_ind(
        pca,
        geom.ind = "point",
        pointshape = 21,
        fill = df_plot$Hclust,
        col.ind = "black",
        alpha.ind = 0.7,
        addEllipses = TRUE
      )
    }
  }

  return(plots)
}
