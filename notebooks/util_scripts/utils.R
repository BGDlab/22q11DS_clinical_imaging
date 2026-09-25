# Shared setup, constants, and helper functions sourced by the analysis scripts
# (04_cohort_demographics.qmd onward). Defines the model-fitting wrappers used
# for case-control, age/sex-interaction, and permutation-based spatial
# correlation tests described in the manuscript's Statistical Analysis section,
# along with the brain-atlas plotting helpers used to build Figures 1-4.

# NOTE: this file deliberately does NOT call rm(list=ls()). A sourced library
# should never clear its caller's environment -- doing so silently destroys the
# user's objects if utils.R is re-sourced mid-session. Scripts that want a clean
# slate should call rm(list=ls()) themselves before sourcing this file.
set.seed(2514)

library(knitr)
library(tidyverse)
library(glue)
library(broom)
library(broom.mixed)
library(nlme)
library(lmerTest)
library(effectsize)
library(cowplot)
library(ggpubr)
library(reticulate)
library(effsize)
library(ggseg)
library(emmeans)
library(table1)
library(flextable)
library(kableExtra)
library(gridExtra)

# Global brain measures modeled throughout the analysis scripts (order matches
# manuscript Table 1 / Figure 1): ICV, cortical GMV, subcortical GMV, WMV,
# ventricular volume, cerebellum volume, total surface area, mean cortical thickness.
global_measures <- c("eTIV", "GMV", "sGMV", "WMV", "Ventricles", "Cerebellum.all",
                     "totalSA", "meanCT")

global_measure_names <- c("ICV", "GMV", "sGMV", "WMV", "Ventricles", "Cerebellum", 
                           "SA", "CT")

# Functions/Themes
base_font <- 9
# Shared ggplot theme applied to all manuscript figure panels for consistent
# typography, sizing, and panel-tag ("A", "B", ...) placement.
figure_theme <- theme(
  text = element_text(family = "Helvetica"),
  plot.title  = element_text(size = base_font+1),
  axis.title  = element_text(size = base_font),
  legend.title = element_text(size = base_font),
  legend.text = element_text(size = base_font-1),
  plot.tag = element_text(size = base_font+2,
                          face = "bold",
                          hjust = 0,
                          vjust = 1),
  legend.key.size = unit(0.5, "cm"),
  plot.margin = margin(t = 2, r = 2, b = 2, l = 2),
  plot.tag.location = "plot",
  plot.tag.position = c(0, 1.02),
)


# Converts a point size to the units ggplot's geom_text()/annotate() expect, so
# in-panel annotation text matches the pt sizes used elsewhere in figure_theme.
annotation_pt <- function(pt) pt / 2.845276

# Main Model ---------------------------------
# Fits the primary case-control (or other two-group) linear model independently
# for each brain feature and returns one row of results per feature, with the
# model t-statistic converted to Cohen's d. This is the workhorse behind the
# global/regional effect-size analyses (Figure 1) and sensitivity analyses.
fit_group_comparison_models <- function(data, cov_str, features, predictor = "dx", group_control = "CN", group_case = "22q11DS",rand_str = "", model_type = "centile",
                      mpr_only_cols = NA, simple = F) {
  
  # Define independent predictor
  data$predictor <- factor(data[[predictor]], 
                           levels = c(group_control,group_case), 
                           labels = c(FALSE, TRUE))
  # Drop NA predictors
  data <- data[!is.na(data$predictor),]

  for (feature in features) {
    # Deviation-score column depends on the requested model type
    if (model_type == "centile") {
      measure <- glue("{feature}Transformed.q.wre")
    } else if (model_type == "zscore") {
      measure <- glue("{feature}Transformed.q.zscore")
    } else if (model_type == "volume") {
      measure <- glue("{feature}Transformed.normalised")
    }
    if (measure %in% colnames(data)) {
      # Cortical thickness models use the MPRAGE-only Euler number, since CT
      # extraction is unreliable outside high-resolution MPRAGE sequences
      if (feature %in% mpr_only_cols) {
        cov_str_tmp <- str_replace_all(cov_str,"euler","euler_mpr")
      } else {
        cov_str_tmp <- cov_str
      }
      form <- as.formula(glue("{measure} ~ predictor{cov_str_tmp}{rand_str}"))
      data_mod <- data[!is.na(data[[measure]]) & !is.infinite(data[[measure]]),]
      # Skip any model that doesn't have any data in one of the levels
      if (any(table(data_mod$predictor) == 0)) {
        next
      }

      # Run lm or lmer
      if (rand_str == "") {
        res <- lm(form, data = data_mod)
      } else {
        res <- lmer(form, data = data_mod)
      }

      # Collect results
      if (simple) {
        res_tidy <- tidy(res) %>%
          filter(term == "predictorTRUE") %>%
          mutate(label = feature)
      } else {
        res_tidy <- tidy(res) %>%
          mutate(label = feature,
                 Predictor = glue("{predictor}_{group_control}_vs_{group_case}"),
                 Case = sum(data_mod$predictor == TRUE),
                 Control = sum(data_mod$predictor == FALSE))
        
        if (rand_str != "") {
          res_tidy$isSingular <- isSingular(res)
        }
      }
      
      # Convert the model t-statistic to Cohen's d
      d <- as_tibble(t_to_d(res_tidy$statistic,df.residual(res))) %>%
        dplyr::select(Cohens_d = d,
                      Cohens_d_ci_lo = CI_low,
                      Cohens_d_ci_hi = CI_high)

      res_tidy <- bind_cols(res_tidy, d)

      if (feature == features[1]) {
        res_all <- res_tidy
      } else {
        res_all <- rbind(res_all,res_tidy)
      }
    }
  }
  
  
  # Final ordering
  res_all <- res_all %>% 
    relocate(p.value, estimate, std.error, statistic, 
             .after = label) %>%
    relocate(Cohens_d,Cohens_d_ci_lo, Cohens_d_ci_hi, .before = p.value)
  
  return(res_all)
}

# Subsets fit_group_comparison_models() output to one term (e.g. the diagnosis effect) and applies
# BeFDR correction across the tested features
summarize_model_results <- function(res, term) {
  res_final <- res[res$term == term,] %>%
                  mutate(p.fdr = p.adjust(p.value, method = "fdr"),
                         sig = factor(p.fdr <= 0.05, 
                              levels = c(FALSE, TRUE), 
                              labels = c("n.s.","p.fdr < 0.05")),
                         sig_label = case_when(
                              p.fdr <= 0.001 ~ "***",
                              p.fdr <= 0.01 ~ "**",
                              p.fdr <= 0.05 ~ "*",
                              .default = NA)) %>%
                  relocate(p.fdr, sig, .after = p.value)
  return(res_final)
}

# Tallies how many features moved between significance/direction categories
# between a sensitivity-analysis result (res_test) and the primary-analysis
# result (res_ref) it's being compared to
summarize_model_changes <- function(res_test, res_ref) {

  # precompute the three row-wise comparisons (res_test and res_ref must be row-aligned)
  sign_consistent <- sign(res_test$Cohens_d) == sign(res_ref$Cohens_d)
  test_sig        <- res_test$sig == "p.fdr < 0.05"
  sig_matches     <- res_test$sig == res_ref$sig

  tibble::tibble(
    category = c(
      "Significant and directionally consistent",
      "Significant but directionally different",
      "No longer significant, but directionally consistent",
      "No longer significant and directionally different",
      "Newly significant, but directionally consistent",
      "Newly significant and directionally different",
      "Never significant and not directionally consistent",
      "Never significant, but directionally consistent"
    ),
    count = c(
      sum(sign_consistent  & test_sig  & sig_matches,  na.rm = TRUE),
      sum(!sign_consistent & test_sig  & sig_matches,  na.rm = TRUE),
      sum(sign_consistent  & !test_sig & !sig_matches, na.rm = TRUE),
      sum(!sign_consistent & !test_sig & !sig_matches, na.rm = TRUE),
      sum(sign_consistent  & test_sig  & !sig_matches, na.rm = TRUE),
      sum(!sign_consistent & test_sig  & !sig_matches, na.rm = TRUE),
      sum(!sign_consistent & !test_sig & sig_matches,  na.rm = TRUE),
      sum(sign_consistent  & !test_sig & sig_matches,  na.rm = TRUE)
    )
  )
}

print_top_significant_results <- function(res, n = 10) {
  res %>% 
    mutate(d = round(Cohens_d,2),
           ci = glue("{round(Cohens_d_ci_lo,2)} to {round(Cohens_d_ci_hi,2)}")) %>%
    dplyr::select(label, d, ci, p.value, p.fdr) %>% 
    arrange(p.value) %>%
    mutate(p.value = signif(p.value,3),
           p.fdr = signif(p.fdr,3)) %>%
    filter(p.fdr <= 0.05) %>%
    head(n) %>%
    kable(digits = Inf,
          col.names = c("Region", "Cohen's d", "95% CI", "p", "p (FDR)"),
          caption = glue("Case-Control Significant Regions in CHOP (Top {n})"))
}


# Summarizes regional effect sizes (min/max Cohen's d, proportion FDR-significant)
# grouped by cortical lobe, for the by-region result tables.
summarize_results_by_cortex <- function(data, caption = caption) {
  cortex_labels <- read_csv("../data/atlas_cortices.csv", show_col_types = F)

  data <- data %>%
    mutate(base_label = str_replace_all(label,"lh_",""),
           base_label = str_replace_all(base_label,"rh_",""),
           base_label = str_replace_all(base_label,"Left_",""),
           base_label = str_replace_all(base_label,"Right_",""),
           measure = case_when(startsWith(base_label,"SA.") ~ "SA",
                               startsWith(base_label,"GM.") ~ "GM",
                               startsWith(base_label,"CT.") ~ "CT",
                               startsWith(base_label,"SUBC.") ~ "GM"),
           base_label = str_remove_all(base_label,"SA\\.|GM\\.|CT\\.|SUBC\\.")) %>% 
    left_join(cortex_labels, by = c("base_label" = "label"))
  
  results_by_cortex <- expand_grid(measure = unique(na.omit(data$measure)),
                                   cortex = unique(na.omit(data$cortex))) %>%
    add_column(min_d = NA, max_d = NA, prop_sig = NA) %>%
    filter(!(measure %in% c("SA","CT") & cortex == "subcortex"))
  
  for (i in 1:nrow(results_by_cortex)) {
    cortex <- results_by_cortex$cortex[i]
    tmp <- data[data$cortex == cortex & 
                  data$measure == results_by_cortex$measure[i] & 
                  !is.na(data$cortex),]
    results_by_cortex$min_d[i] <- min(tmp$Cohens_d, na.rm = T)
    results_by_cortex$max_d[i] <- max(tmp$Cohens_d, na.rm = T)
    results_by_cortex$prop_sig[i] <- sum(tmp$p.fdr <= 0.05)/nrow(tmp)
  }
  
  kable(results_by_cortex, col.names = c("Measure","Region", "Cohen's d (Min)",
                                   "Cohen's d (Max)", "Prop. Sig"),
        caption = caption)
}
# Outliers ------------------------------------
# Renders brain maps of a single value (e.g. proportion of patients below the
# 2.5th centile) across GM/SUBC/SA/CT categories, used for the extreme-phenotype
# maps in Figure 2.
plot_extreme_deviation_maps <- function(data, fill_col, color_scheme, plot_lims, fill_name,
                                 colorbar_title = "",sig_col = NA, cats = c("GM","SUBC","SA","CT"),
                                 cort_only = F, add_colorbar = T, add_fillbar = T, exclude_y_label = F,
                                 formatted = T, lh_only = F, p_col = "p.value",
                                 percent = T) {
  # Make plots for GM, CT, and SA
  if (add_colorbar) {
    leg <- "bottom"
  } else {
    leg <- "none"
  }
  plots <- list()
  for (m_cat in cats) {
    cat_name <- case_match(m_cat,
                           "GM" ~ "Gray Matter Volume",
                           "SUBC" ~ "Subcortical Volume",
                           "SA" ~ "Surface Area",
                           "CT" ~ "Cortical Thickness")
    
    data_plot <- data %>% 
      filter(grepl(glue("{m_cat}\\."),label)) %>%
      mutate(label = str_remove_all(label, glue("{m_cat}\\."))) 
    
    if (!is.na(sig_col)) {
      data_plot <- data_plot %>%
        arrange(get(sig_col))
    }
    if (m_cat == "SUBC") {
      plot   <- plot_brain_atlas(data_plot, aseg_plot, fill_col, color_scheme, plot_lims,
                           colorbar_title = colorbar_title,sig_col = sig_col,
                           p_col = p_col, percent = percent)
    } else {
      plot   <- plot_brain_atlas(data_plot, dk_plot, fill_col, color_scheme, plot_lims,
                           colorbar_title = colorbar_title,sig_col = sig_col,
                           position = "stacked", lh_only = lh_only,
                           p_col = p_col, percent = percent)
    }
    
    # Cortical Plot
    plots[[glue("{m_cat} {fill_col}")]] <- plot$plot +
      ggtitle(" ") + ylab(" ") +
      theme(legend.position = "none",
            plot.title.position = "plot") + figure_theme +
      theme(legend.justification = c(1, 0),
            legend.key.size = unit(0.5, "cm"))
    if (!add_colorbar) {
      plots[[glue("{m_cat} {fill_col}")]] <- plots[[glue("{m_cat} {fill_col}")]] + guides(color = "none")
    }
    if (!add_fillbar) {
      plots[[glue("{m_cat} {fill_col}")]] <- plots[[glue("{m_cat} {fill_col}")]] + guides(fill = "none")
    }

    if (m_cat == cats[1] & exclude_y_label == FALSE) {
      plots[[glue("{m_cat} {fill_col}")]] <- plots[[glue("{m_cat} {fill_col}")]] +
        ylab(fill_name)
    }
    plots[[glue("{m_cat} {fill_col}")]] <- plots[[glue("{m_cat} {fill_col}")]] +
      ggtitle(cat_name)
  }

  if (formatted) {
    if ("SUBC" %in% cats) {
      if (length(cats) == 4) {
        plot_n_row <- 1
        plot_n_col <- 4
        w = c(1,1,1,0.9)
      } else if (length(cats) == 2) {
        plot_n_row <- 1
        plot_n_col <- 2
        w = c(1,0.75)
      } else if (length(cats) == 3) {
        plot_n_row <- 1
        plot_n_col <- 3
        w = c(1,1,0.75)
      }
    } else {
      plot_n_row <- 1
      plot_n_col <- 3
      w = 1
    }
    p <- ggarrange(plotlist =  plots, nrow = plot_n_row, ncol = plot_n_col, 
                   common.legend = TRUE, legend = "bottom", widths = w,
                   align = "v")
  } else {
    return(plots)
  }
}

# Plotting ---------------------------------
# Core ggseg brain-map renderer
plot_brain_atlas <- function(data, atlas, fill_col = "Cohens_d",  color_scheme = c("royalblue","white","firebrick"),
                       plot_lims = c(-10,1), colorbar_title = "",sig_col = NA, p_col = "p.value",
                       position = "identity", lh_only = F, skip_sagittal = T, percent = F) {
  
  #If seeking to plot results by significance, columns of atlas need to be ordered 
  # such that the significant borders take precedence over the n.s. border
  atlas <- merge(atlas, data, by = "label", all.x = T)
  if (!is.na(sig_col)) {
    atlas <- atlas[order(!is.na(atlas[[p_col]]),atlas[[sig_col]], -atlas[[p_col]]),]
    if (sum(atlas[[sig_col]] == "n.s.", na.rm = T) == 0) {
      fill_colors <- c("black")
    } else {
      fill_colors <- c("grey","black")
    }
  } else {
    atlas$sig_tmp <- FALSE
    sig_col <- "sig_tmp"
    fill_colors <- c("black","grey")
  }
  if (skip_sagittal) {
    atlas <- atlas %>% filter(side != "sagittal")
  }
  
  #Generate ggplot
  if (lh_only) {
    if (position == "stacked") {
      plot_out <- ggplot() + 
        geom_brain(atlas=atlas, size = 0.5, hemi = "left",
                   position = position_brain(hemi + side ~ .),
                   mapping=aes(fill  = get(fill_col),
                               color = get(sig_col)))
    } else {
      plot_out <- ggplot() + 
        geom_brain(atlas=atlas, size = 0.5, hemi = "left",
                   position = position_brain(. ~ hemi + side),
                   mapping=aes(fill  = get(fill_col),
                               color = get(sig_col)))
    }
  } else if (position == "identity") {
    plot_out <- ggplot() + 
      geom_brain(atlas=atlas, size = 0.5, 
                 mapping=aes(fill  = get(fill_col),
                             color = get(sig_col)))
  } else if (position == "stacked") {
    plot_out <- ggplot() + 
      geom_brain(atlas=atlas, size = 0.5, 
                 position = position_brain(hemi ~ side),
                 mapping=aes(fill  = get(fill_col),
                             color = get(sig_col)))
  } else if (position == "stacked2") {
    plot_out <- ggplot() + 
      geom_brain(atlas=atlas, size = 0.5, 
                 position = position_brain(side ~ hemi),
                 mapping=aes(fill  = get(fill_col),
                             color = get(sig_col)))
  }
  plot_out <- plot_out +
    scale_color_manual(, values = fill_colors, na.value="grey", 
                       guide = "none",
                       na.translate=FALSE) +
    theme_brain2(text.colour = "black", text.family = "")  + 
    guides(fill = guide_colourbar(title = colorbar_title,
                                  title.position = "left",
                                  title.vjust = 0.8,
                                  order = 1),
           color = guide_legend(title = "",
                                order = 2,
                                override.aes = list(fill = NA))) +
    theme(legend.title = element_text(margin = margin(r = 10)))
  
  if (percent) {
    plot_out <- plot_out +
      scale_fill_gradientn(colours = color_scheme, na.value="grey",
                           limits = plot_lims, breaks = c(plot_lims[1],0,plot_lims[2]),
                           labels = scales::label_percent())
  } else {
    plot_out <- plot_out +
      scale_fill_gradientn(colours = color_scheme, na.value="grey",
                           limits = plot_lims, breaks = c(plot_lims[1],0,plot_lims[2]))
    
  }
  
  
  # Save legend
  plot_out <- plot_out + figure_theme
  plot_legend <- get_legend(plot_out)
  plot_out <- plot_out  + theme(legend.position="none")
  
  return(list(plot = plot_out, legend = plot_legend))
}

# Assembles the full main-figure layout (global panel + subcortical atlas +
# GM/SA/CT cortical atlases, with a shared legend) used for Figure 1
assemble_main_figure_layout <- function(data, raw_data, fill_col, color_scheme, plot_lims,
                          plot_type = "zscore", comp_col = "dx", comp_col_1 = "CN", comp_col_2 = "22q11DS",
                          colorbar_title = "",sig_col = NA, cats = c("GM","SA","CT"),
                          cort_only = FALSE, formatted = TRUE, bar_global = F) {
  
  # Make plots for GM, CT, and SA
  plots_cortical <- list("Blank" = NULL)
  
  for (cat in cats) {
    cat_name <- case_match(cat,
                           "GM" ~ "Gray Matter Volume",
                           "SA" ~ "Surface Area",
                           "CT" ~ "Cortical Thickness")
    
    data_plot <- data %>% 
      filter(grepl(glue("{cat}\\."),label)) %>%
      mutate(label = str_remove_all(label, glue("{cat}\\.")))
    
    # Cortical Plot
    plot_cortical   <- plot_brain_atlas(data_plot, dk_plot, fill_col, color_scheme, plot_lims,
                                  colorbar_title = colorbar_title,sig_col = sig_col,
                                  position = "stacked")
    
    plots_cortical[[cat]] <- plot_cortical$plot + 
      ggtitle(cat_name) + 
      theme(legend.position = "bottom",
            plot.title.position = "plot",
            legend.text=element_text(size=7))
    
    leg <- as_ggplot(get_legend(plots_cortical[[cat]]))
    
  }
  
  p_cort <- ggarrange(plotlist = plots_cortical, ncol = 4,
                      common.legend = T, widths = c(0.15,1,1,1),
                      legend = "none", labels = c("C","","",""))
  
  # Add subcortical plot
  data_plot <- data %>% 
    filter(grepl(glue("SUBC\\."),label)) %>%
    mutate(label = str_remove_all(label, glue("SUBC\\.")))
  plot_subcortical   <- plot_brain_atlas(data_plot, aseg_plot, fill_col, color_scheme, plot_lims,
                                   colorbar_title = colorbar_title,sig_col = sig_col)
  
  plot_subcortical <- plot_subcortical$plot + 
    ggtitle("Subcortical GMV") + 
    theme(legend.position = "none",
          plot.title.position = "plot",
          legend.text=element_text(size=7),
          plot.margin = margin(l = 25, r = 5.5, t = 15, b = 5.5, unit = "pt"))
  
  # Plot Global
  if (plot_type == "zscore") {
    raw_data_suf <- "Transformed.q.zscore"
    sig_pos <- 4.7
    plot_ylab <- "Deviation Score"
  } else if (plot_type == "centile") {
    raw_data_suf <- "Transformed.q.wre"
    sig_pos <- 1.0
    plot_ylab <- "Centile"
  } else if (plot_type == "volume") {
    raw_data_suf <- "Transformed.normalised"
    sig_pos <- 75
    plot_ylab <- "Normalized Value"
  }
  if (bar_global == F) {
    plot_global <- plot_global_measures(data, raw_data, plot_type, comp_col, comp_col_1, comp_col_2) +
      theme(plot.margin = margin(l = 30, r = 5.5, t = 10, b = 5.5, unit = "pt"),
            legend.margin = margin(0,0,0,0, unit = "pt"),
            legend.title=element_blank())
  } else {    
    tmp <- data %>% 
            mutate(label = str_replace_all(label,"\\.all",""),
                   label = str_replace_all(label,"_Cortex",""),
                   label = str_replace_all(label,"_"," "),
                   label = str_replace_all(label,"eTIV","ICV"),
                   label = str_replace_all(label,"totalSA","SA"),
                   label = str_replace_all(label,"meanCT","CT"),
                   label = factor(label, levels = global_measure_names)) %>%
            filter(!is.na(label))
    
    plot_global <- ggplot(tmp, aes(x = label, y = .data[[fill_col]], 
                                    ymin = .data[[ci_lo_col]], 
                                    ymax = .data[[ci_hi_col]])) +
      geom_col(aes(fill = .data[[sig_col]]), color = "black") +
      geom_errorbar(width = 0.5) + 
      geom_hline(yintercept = 0, linetype = "dashed") +
      ylab(colorbar_title) +
      theme_pubclean() +
      scale_fill_manual(values = c("white","grey"), drop=FALSE) +
      xlab("") + 
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
            legend.position="none")  + 
      guides(fill=guide_legend(title=""))
  }
  
  if (cort_only) {
    return(p_cort)
  } else if (!formatted) {
    return(list("Global" = plot_global,
                "Subcortical" = p_cort,
                "Cortical" = plots_cortical,
                "Legend" = leg))
  } else {
    p_subcort <- ggarrange(plot_global,plot_subcortical,  ncol = 2, labels = c("A","B"),
                           widths = c(0.7, 0.3))
    p <- ggarrange(p_subcort, p_cort, leg, nrow = 3, heights = c(1.2, 1, 0.2))
    return(p)
  }
}

# Boxplots of the global brain measures (ICV, GMV, sGMV, WMV, Ventricles,
# Cerebellum, SA, CT) by group, with FDR-significance brackets from `data`
plot_global_measures <- function(data, raw_data, plot_type, comp_col, comp_col_1, comp_col_2) {

  # Plot Global
  if (plot_type == "zscore") {
    raw_data_suf <- "Transformed.q.zscore"
    sig_pos <- 4.7
    plot_ylab <- "Deviation Score"
  } else if (plot_type == "centile") {
    raw_data_suf <- "Transformed.q.wre"
    sig_pos <- 1.0
    plot_ylab <- "Centile"
  } else if (plot_type == "volume") {
    raw_data_suf <- "Transformed.normalised"
    sig_pos <- 75
    plot_ylab <- "Normalized Value"
  }
  
  tmp <- raw_data %>% 
    select(participant, any_of(comp_col), sex, age_days, (ends_with(raw_data_suf) &
                                                            starts_with(global_measures))) %>%
    pivot_longer(cols = contains(raw_data_suf), names_to = "Region", values_to = "Value") %>%
    mutate(Region = str_remove_all(Region,raw_data_suf),
           Region = factor(Region, 
                           levels = global_measures,
                           labels = global_measure_names)) %>%
    na.omit()
  
  sig_tab <- data %>% 
    select(label, Cohens_d, p.value, p.fdr) %>% 
    filter(label %in% global_measures) %>% 
    mutate(group1 = comp_col_1, group2 = comp_col_2,
           .y. = "Value") %>% 
    select(Region = label,.y., group1, group2, p.value, p.fdr) %>%
    mutate(p.signif = case_when(p.fdr <= 0.001 ~ "***",
                                p.fdr <= 0.01 ~ "**",
                                p.fdr <= 0.05 ~ "*"),
           Region = reduce2(global_measures, global_measure_names,  .init = Region,
                            ~ str_replace(..1, ..2, ..3)),
           y.position = sig_pos)

  plot_global <- ggplot(tmp, aes(x = Region, y = Value)) +
    geom_point(aes(color = .data[[comp_col]]),
               position = position_jitterdodge(), size = 0.5) +
    geom_boxplot(aes(fill = .data[[comp_col]]),
                 outliers = F, alpha = 0.3) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    ylab(plot_ylab) +
    theme_pubclean() +
    scale_fill_manual(values = c("transparent","transparent"), drop=FALSE) +
    scale_color_manual(values = c("darkgrey","firebrick"), drop=FALSE) +
    xlab("") +
    stat_pvalue_manual(sig_tab, x = "Region", label = "p.signif") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          legend.position = "bottom",
          legend.justification = "right",
          legend.box.spacing = unit(-2, "pt"),
          legend.margin = margin(0, 0, 0, 0, unit = "pt"),
          legend.spacing.x = unit(2, "pt"))

  return(plot_global)
}

# Interaction Models ----
# Fits predictor*interactor models per feature
fit_interaction_model <- function(data, cov_str, features, rand_str = "",
                                       predictor = "dx", group_control = "CN", 
                                       group_case = "22q11DS", interactor = "age_days",  
                                       model_type = "centile", gam = F,
                                       mpr_only_cols = NA) {
  
  # Define independent predictor
  data$predictor <- factor(data[[predictor]], 
                           levels = c(group_control,group_case), 
                           labels = c(group_control, group_case))
  
  data$interactor <- data[[interactor]]
  
  # Drop NA predictors
  data <- data[!is.na(data$predictor) & !is.na(data$interactor),]
  
  prop_thresh_low <- 0.05
  prop_thresh_high <- 1 - prop_thresh_low
  
  if (model_type == "zscore") {
    prop_thresh_low <- qnorm(prop_thresh_low)
    prop_thresh_high <- qnorm(prop_thresh_high)
  }
  
  for (feature in features) {
    # select centile feature
    if (model_type == "centile") {
      measure <- glue("{feature}Transformed.q.wre")
    } else if (model_type == "zscore") {
      measure <- glue("{feature}Transformed.q.zscore")
    } else if (model_type == "volume") {
      measure <- glue("{feature}Transformed.normalised")
    }
    
    if (feature %in% mpr_only_cols) {
      cov_str_tmp <- str_replace_all(cov_str,"euler","euler_mpr")
    } else {
      cov_str_tmp <- cov_str
    }
    
    # Define model predictors + covs + rand effects
    model_str <- glue("~ predictor*interactor{cov_str_tmp}{rand_str}")

    if (measure %in% colnames(data)) {
      form <- as.formula(glue("{measure}{model_str}"))
      print(form)
      data_mod <- data[!is.na(data[[measure]]),]
      #Skip any model that doesn't have any data in one of the levels
      if (any(table(data_mod$predictor) == 0)) {
        next
      }
      # Run lm or lmer
      if (rand_str == "") {
        res <- lm(form, data = data_mod)
      } else if (gam) {
        res <- gam(form, data = data_mod)
      }else {
        res <- lmer(form, data = data_mod)
      }
      
      # Collect results
      res_tidy <- tidy(res)
      
      #Cohen's D
      d <- as_tibble(t_to_d(res_tidy$statistic,res$df.residual)) %>%
        dplyr::select(Cohens_d = d,
                      Cohens_d_ci_lo = CI_low,
                      Cohens_d_ci_hi = CI_high)
      
      res_tidy <- bind_cols(res_tidy, d)
      res_tidy <- res_tidy %>% relocate(estimate, std.error, statistic, p.value, .after = Cohens_d_ci_hi)
      
      res_tidy$BIC <- BIC(res)
      res_tidy$label <- feature
      res_tidy <- res_tidy %>% relocate(label, .before = term)
      
      if (feature == features[1]) {
        res_all <- res_tidy
      } else {
        res_all <- rbind(res_all,res_tidy)
      }
    }
  }
  return(res_all)
}
# Correlation Functions ---------
# Freedman-Lane permutation of case-control Cohen's d, used to build the null
# distribution for the CHOP vs ENIGMA-22q spatial correlation test
freedman_lane_perm <- function(df_perm, perm_features, n_perm = 1000,
                               n_core = 32, suf = "Transformed.q.zscore",
                               covs = "sex + age_days + euler", mixed = F) {

  df_perm <- df_perm %>% select(participant, sex, dx, age_days, age_days_sq, euler,
                                any_of(glue("{perm_features}{suf}")),
                                contains("eTIV"), site) %>%
                                na.omit()

  df_perm[sapply(df_perm, is.infinite)] <- NA
  perms <- replicate(n_perm, sample(1:nrow(df_perm)))
  out_perms <- as_tibble(matrix(NA, nrow = length(perm_features),
                                ncol = n_perm + 1), .name_repair = "minimal")
  colnames(out_perms) <- c("label",glue("perm_{1:n_perm}"))
  out_perms$label <- perm_features
  cluster <- makeCluster(n_cores)
  registerDoParallel(cluster)

  for (n in 1:nrow(out_perms)) {
    results <- list()
    feature <- perm_features[n]
    # Fit the reduced model
    if (mixed) {
      form <- as.formula(glue("{feature}{suf} ~ {covs} + (1 | site)"))
      reduced_model <- lmer(form, df_perm)
    } else {
      form <- as.formula(glue("{feature}{suf} ~ {covs}"))
      reduced_model <- lm(form, df_perm)
    }

    # Permute the residuals and combine with fitted effects
    results <- foreach(i = 1:ncol(perms),
                       .packages = c("tidyverse","glue","broom","effectsize"),
                       .export = c("feature", "suf","covs")) %dorng% {
                         # Calculate permuted outcome
                         permuted_resid <- residuals(reduced_model)[perms[,i]] + fitted.values(reduced_model)
                         df_perm$perm_outcome <- permuted_resid

                         # Refit model to permuted outcome
                         if (mixed) {
                           form <- as.formula(glue("perm_outcome ~ dx + {covs} + (1 | site)"))
                           permuted_model <- lmer(form, df_perm)
                         } else {
                           form <- as.formula(glue("perm_outcome ~ dx + {covs}"))
                           permuted_model <- lm(form, df_perm)
                         }

                         # Calculate permuted Cohen's d
                         results[i] <- t_to_d(coef(summary(permuted_model))["dx22q11DS","Pr(>|t|)"],
                                              df.residual(permuted_model))$d
                       }
    out_perms[n,glue("perm_{1:n_perm}")] <- as.list(unlist(results))
  }

  stopCluster(cluster)
  return(out_perms)
}

# Compares an observed spatial correlation (true_brain_map vs comp_brain_map) to
# the null distribution generated by freedman_lane_perm(), per tissue category.
freedman_lane_corr_test <- function(true_brain_map, perms, comp_brain_map,
                        term_1 = "Cohens_d", term_2 = "Cohens_d",
                        cats = c("GM","SA","CT"), method = "pearson",
                        direction = "either") {
  out_res <- tibble(measure = cats, corr = NA, pval = NA)
  for (cat in cats) {
    # Extract sub-measures for the true map
    true_map_cat <- true_brain_map %>% filter(grepl(glue("{cat}\\."),label))

    # Extract sub-measures for the comparison map
    comp_map_cat <- comp_brain_map %>%
      filter(grepl(glue("{cat}\\."),label)) %>%
      arrange(match(label, true_map_cat$label))

    # Extract sub-measures for the null maps
    null_map_cat <- perms %>%
      filter(grepl(glue("{cat}\\."),label)) %>%
      arrange(match(label, true_map_cat$label)) %>%
      select(-label)

    # Calculate the true correlation
    true_corr <- cor(true_map_cat %>% pull(term_1),
                     comp_map_cat %>% pull(term_2),
                     method = method)

    # Calculate the null correlations
    null_corrs <- cor(null_map_cat,
                      comp_map_cat %>% pull(term_2),
                      method = method)

    if (direction == "either") {
      pval <- sum(abs(null_corrs) >= abs(true_corr))/nrow(null_corrs)
    } else if (direction == "greater") {
      pval <- sum(null_corrs >= true_corr)/nrow(null_corrs)
    }  else if (direction == "lesser") {
      pval <- sum(null_corrs <= true_corr)/nrow(null_corrs)
    } 
    
    # Save results
    out_res[out_res$measure == cat,"corr"] <- true_corr 
    out_res[out_res$measure == cat,"pval"] <- pval 
  }
  return(out_res)
}

