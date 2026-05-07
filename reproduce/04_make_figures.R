### 04_make_figures.R --- 
#----------------------------------------------------------------------
## 
### Code:

library(ggplot2)
library(gridExtra)

source("./simulations/summarize.output.R")

###############################################################################
##### 4A
###############################################################################

results.4A <- 
    rbind(
        rbind(
            summarize.output(sim.setting = "4A", cens.percentage = "high", intervention.A = 1, n = 500, M = 1000, standard.np = TRUE),
            summarize.output(sim.setting = "4A", cens.percentage = "high", intervention.A = 1, n = 500, M = 1000, print.oracle.coverage = TRUE),
            summarize.output(sim.setting = "4A", cens.percentage = "high", intervention.A = 1, n = 500, M = 1000, baseline.tmle = TRUE),
            summarize.output(sim.setting = "4A", cens.percentage = "high", intervention.A = 1, n = 500, M = 1000, misspecify.T = TRUE, misspecify.C = FALSE),
            summarize.output(sim.setting = "4A", cens.percentage = "high", intervention.A = 1, n = 500, M = 1000, misspecify.T = 2, misspecify.C = 2),
            summarize.output(sim.setting = "4A", cens.percentage = "high", intervention.A = 1, n = 500, M = 1000, use.hal = "no-interaction-3")## ,
        )[, cens := "~40% cens."],
        rbind(
            summarize.output(sim.setting = "4A", cens.percentage = "30%", intervention.A = 1, n = 500, M = 1000, standard.np = TRUE),
            summarize.output(sim.setting = "4A", cens.percentage = "30%", intervention.A = 1, n = 500, M = 1000, print.oracle.coverage = TRUE),
            summarize.output(sim.setting = "4A", cens.percentage = "30%", intervention.A = 1, n = 500, M = 1000, baseline.tmle = TRUE),
            summarize.output(sim.setting = "4A", cens.percentage = "30%", intervention.A = 1, n = 500, M = 1000, misspecify.T = TRUE, misspecify.C = FALSE),
            summarize.output(sim.setting = "4A", cens.percentage = "30%", intervention.A = 1, n = 500, M = 1000, misspecify.T = 2, misspecify.C = 2),
            summarize.output(sim.setting = "4A", cens.percentage = "30%", intervention.A = 1, n = 500, M = 1000, use.hal = "no-interaction-3")## ,
        )[, cens := "~30% cens."],
        rbind(
            summarize.output(sim.setting = "4A", cens.percentage = "low", intervention.A = 1, n = 500, M = 1000, standard.np = TRUE),
            summarize.output(sim.setting = "4A", cens.percentage = "low", intervention.A = 1, n = 500, M = 1000, print.oracle.coverage = TRUE),
            summarize.output(sim.setting = "4A", cens.percentage = "low", intervention.A = 1, n = 500, M = 1000, baseline.tmle = TRUE),
            summarize.output(sim.setting = "4A", cens.percentage = "low", intervention.A = 1, n = 500, M = 1000, misspecify.T = TRUE, misspecify.C = FALSE),
            summarize.output(sim.setting = "4A", cens.percentage = "low", intervention.A = 1, n = 500, M = 1000, misspecify.T = 2, misspecify.C = 2),
            summarize.output(sim.setting = "4A", cens.percentage = "low", intervention.A = 1, n = 500, M = 1000, use.hal = "no-interaction-3")## ,
        )[, cens := "~20% cens."],
        rbind(
            summarize.output(sim.setting = "4A", cens.percentage = "10%", intervention.A = 1, n = 500, M = 1000, standard.np = TRUE),
            summarize.output(sim.setting = "4A", cens.percentage = "10%", intervention.A = 1, n = 500, M = 1000, print.oracle.coverage = TRUE),
            summarize.output(sim.setting = "4A", cens.percentage = "10%", intervention.A = 1, n = 500, M = 1000, baseline.tmle = TRUE),
            summarize.output(sim.setting = "4A", cens.percentage = "10%", intervention.A = 1, n = 500, M = 1000, misspecify.T = TRUE, misspecify.C = FALSE),
            summarize.output(sim.setting = "4A", cens.percentage = "10%", intervention.A = 1, n = 500, M = 1000, misspecify.T = 2, misspecify.C = 2),
            summarize.output(sim.setting = "4A", cens.percentage = "10%", intervention.A = 1, n = 500, M = 1000, use.hal = "no-interaction-3")## ,
        )[, cens := "~10% cens."]
    )

results.4A <- melt(results.4A, id.vars = c("sim.setting", "cens", "estimator", "which"))

results.4A[estimator == "standard np", estimator := "unadjusted NP"]
results.4A[estimator == "baseline-tmle", estimator := "baseline TMLE"]
results.4A[estimator == "tmle", estimator := "TMLE"]
results.4A[estimator == "hal-tmle", estimator := "HAL-TMLE"]

results.4A[which == "correct" & estimator == "TMLE", estimator := paste0("cox-", estimator, " (a)")]
results.4A[which == "missp" & estimator == "TMLE", estimator := paste0("cox-", estimator, " (b)")]
results.4A[which == "both missp" & estimator == "TMLE", estimator := paste0("cox-", estimator, " (c)")]

results.4A[, value := as.numeric(value)]

results.4A[, line := 0]
results.4A[variable == "bias", line := 0]
results.4A[variable == "cov", line := 0.95][variable == "oracle.cov", line := 0.95]
results.4A[variable == "mse" & !(estimator %in% grep("g-comp", estimator, value = TRUE)),
           line := min(value), by = c("cens", "variable")]
results.4A[, line2 := line][variable == "cov", line2 := 1][variable == "oracle.cov", line2 := 1]
results.4A[variable == "mse", line2 := NA]
results.4A[variable == "mse", line := NA]

results.4A[, cens := factor(cens, levels = c("~10% cens.", "~20% cens.", "~30% cens.", "~40% cens."))]
results.4A[variable == "cov", variable := "cov. (95%)"]
results.4A[variable == "oracle.cov", variable := "oracle cov. (95%)"]

results.4A[, estimator := factor(estimator, levels = c("unadjusted NP",
                                                       "baseline TMLE",
                                                       "cox-TMLE (c)",
                                                       "cox-TMLE (b)",
                                                       "cox-TMLE (a)",
                                                       "HAL-TMLE"))]

results.4A[, variable := factor(variable, levels = c("bias", "mse", "oracle cov. (95%)", "cov. (95%)"))]


(p.4A <- ggplot(results.4A[#cens == "low" &
     !is.na(estimator) &
     variable %in% c("bias", "mse", "oracle cov. (95%)", "cov. (95%)")]) +
     geom_vline(aes(xintercept = line2), alpha = 0, linetype = "dashed", color = "red") + 
     geom_vline(aes(xintercept = line), alpha = 0.3, linetype = "dashed", color = "red") + 
     geom_point(aes(y = estimator, x = value, shape = cens, color = cens)) +
     theme_bw() + xlab("") + ylab("") +
     ggtitle("Setting 1: Nonlinear dependence on past number of recurrent events") + 
     facet_grid(.~variable, scales = "free") +
     theme(strip.background = element_rect(fill = 'white', colour = 'black'),
           plot.caption = element_text(size = 12, face = "bold"),
           legend.title = element_blank(),
           legend.position = "top",
           legend.justification="right",
           legend.margin = margin(0,0,0,0),
           legend.box.margin = margin(-20,30,-5,-5),
           plot.title = element_text(hjust = -0.315)) + 
     labs(caption = expression("cox-TMLE (a) = all correctly specified, cox-TMLE (b) = "*lambda^{'y'}*" misspecified, cox-TMLE (c) = "*lambda^{'c'}*","*lambda^{'y'}*" misspecified")))

ggsave(file = paste0("./figures/",
                     "fig-simulation-results-4A-2025",
                     ".pdf"),
       p.4A, width = 26, height = 8, units = "cm")



###############################################################################
##### 7A
###############################################################################

results.7A <- 
    rbind(
        rbind(
            summarize.output(sim.setting = "7A", cens.percentage = "high", intervention.A = 1, n = 500, M = 1000, standard.np = TRUE),
            summarize.output(sim.setting = "7A", cens.percentage = "high", intervention.A = 1, n = 500, M = 1000, print.oracle.coverage = TRUE),
            summarize.output(sim.setting = "7A", cens.percentage = "high", intervention.A = 1, n = 500, M = 1000, baseline.tmle = TRUE),
            summarize.output(sim.setting = "7A", cens.percentage = "high", intervention.A = 1, n = 500, M = 1000, misspecify.T = 2, misspecify.C = TRUE),
            summarize.output(sim.setting = "7A", cens.percentage = "high", intervention.A = 1, n = 500, M = 1000, misspecify.T = FALSE, misspecify.C = TRUE),
            summarize.output(sim.setting = "7A", cens.percentage = "high", intervention.A = 1, n = 500, M = 1000, use.hal = "no-interaction-8")## ,
        )[, cens := "~40% cens."],
        rbind(
            summarize.output(sim.setting = "7A", cens.percentage = "30%", intervention.A = 1, n = 500, M = 1000, standard.np = TRUE),
            summarize.output(sim.setting = "7A", cens.percentage = "30%", intervention.A = 1, n = 500, M = 1000, print.oracle.coverage = TRUE),
            summarize.output(sim.setting = "7A", cens.percentage = "30%", intervention.A = 1, n = 500, M = 1000, baseline.tmle = TRUE),
            summarize.output(sim.setting = "7A", cens.percentage = "30%", intervention.A = 1, n = 500, M = 1000, misspecify.T = 2, misspecify.C = TRUE),
            summarize.output(sim.setting = "7A", cens.percentage = "30%", intervention.A = 1, n = 500, M = 1000, misspecify.T = FALSE, misspecify.C = TRUE),
            summarize.output(sim.setting = "7A", cens.percentage = "30%", intervention.A = 1, n = 500, M = 1000, use.hal = "no-interaction-8")## ,
        )[, cens := "~30% cens."],
        rbind(
            summarize.output(sim.setting = "7A", cens.percentage = "low", intervention.A = 1, n = 500, M = 1000, standard.np = TRUE),
            summarize.output(sim.setting = "7A", cens.percentage = "low", intervention.A = 1, n = 500, M = 1000, print.oracle.coverage = TRUE),
            summarize.output(sim.setting = "7A", cens.percentage = "low", intervention.A = 1, n = 500, M = 1000, baseline.tmle = TRUE),
            summarize.output(sim.setting = "7A", cens.percentage = "low", intervention.A = 1, n = 500, M = 1000, misspecify.T = 2, misspecify.C = TRUE),
            summarize.output(sim.setting = "7A", cens.percentage = "low", intervention.A = 1, n = 500, M = 1000, misspecify.T = FALSE, misspecify.C = TRUE),
            summarize.output(sim.setting = "7A", cens.percentage = "low", intervention.A = 1, n = 500, M = 1000, use.hal = "no-interaction-8")## ,
        )[, cens := "~20% cens."],
        rbind(
            summarize.output(sim.setting = "7A", cens.percentage = "10%", intervention.A = 1, n = 500, M = 1000, standard.np = TRUE),
            summarize.output(sim.setting = "7A", cens.percentage = "10%", intervention.A = 1, n = 500, M = 1000, print.oracle.coverage = TRUE),
            summarize.output(sim.setting = "7A", cens.percentage = "10%", intervention.A = 1, n = 500, M = 1000, baseline.tmle = TRUE),
            summarize.output(sim.setting = "7A", cens.percentage = "10%", intervention.A = 1, n = 500, M = 1000, misspecify.T = 2, misspecify.C = TRUE),
            summarize.output(sim.setting = "7A", cens.percentage = "10%", intervention.A = 1, n = 500, M = 1000, misspecify.T = FALSE, misspecify.C = TRUE),
            summarize.output(sim.setting = "7A", cens.percentage = "10%", intervention.A = 1, n = 500, M = 1000, use.hal = "no-interaction-8")## ,
        )[, cens := "~10% cens."]
    )

results.7A <- melt(results.7A, id.vars = c("sim.setting", "cens", "estimator", "which"))

results.7A[estimator == "standard np", estimator := "unadjusted NP"]
results.7A[estimator == "baseline-tmle", estimator := "baseline TMLE"]
results.7A[estimator == "tmle", estimator := "TMLE"]
results.7A[estimator == "hal-tmle", estimator := "HAL-TMLE"]

results.7A[which == "correct" & estimator == "TMLE", estimator := paste0("cox-", estimator, " (a)")]
results.7A[which == "cens missp" & estimator == "TMLE", estimator := paste0("cox-", estimator, " (b)")]
results.7A[which == "both missp" & estimator == "TMLE", estimator := paste0("cox-", estimator, " (c)")]

results.7A[, value := as.numeric(value)]

results.7A[, line := 0]
results.7A[variable == "bias", line := 0]
results.7A[variable == "cov", line := 0.95][variable == "oracle.cov", line := 0.95]
results.7A[variable == "mse" & !(estimator %in% grep("g-comp", estimator, value = TRUE)),
           line := min(value), by = c("cens", "variable")]
results.7A[, line2 := line][variable == "cov", line2 := 1][variable == "oracle.cov", line2 := 1]
results.7A[variable == "mse", line2 := NA]
results.7A[variable == "mse", line := NA]

results.7A[, cens := factor(cens, levels = c("~10% cens.", "~20% cens.", "~30% cens.", "~40% cens."))]
results.7A[variable == "cov", variable := "cov. (95%)"]
results.7A[variable == "oracle.cov", variable := "oracle cov. (95%)"]

results.7A[, estimator := factor(estimator, levels = c("unadjusted NP",
                                                       "baseline TMLE",
                                                       "cox-TMLE (c)",
                                                       "cox-TMLE (b)",
                                                       "cox-TMLE (a)",
                                                       "HAL-TMLE"))]

results.7A[, variable := factor(variable, levels = c("bias", "mse", "oracle cov. (95%)", "cov. (95%)"))]


(p.7A <- ggplot(results.7A[#cens == "low" &
     !is.na(estimator) &
     variable %in% c("bias", "mse", "oracle cov. (95%)", "cov. (95%)")]) +
     geom_vline(aes(xintercept = line2), alpha = 0, linetype = "dashed", color = "red") + 
     geom_vline(aes(xintercept = line), alpha = 0.3, linetype = "dashed", color = "red") + 
     geom_point(aes(y = estimator, x = value, shape = cens, color = cens)) +
     theme_bw() + xlab("") + ylab("") +
     ggtitle("Setting 2: Unobserved heterogeneity") + 
     facet_grid(.~variable, scales = "free") +
     theme(strip.background = element_rect(fill = 'white', colour = 'black'),
           plot.caption = element_text(size = 12, face = "bold"),
           legend.title = element_blank(),
           legend.position = "top",
           legend.justification="right",
           legend.margin = margin(0,0,0,0),
           legend.box.margin = margin(-20,30,-5,-5),
           plot.title = element_text(hjust = -0.185)) + 
     labs(caption = expression("cox-TMLE (a) = "*lambda^{'y'}*" partly correct and "*lambda^{'c'}*" correct, cox-TMLE (b) = "*lambda^{'c'}*" misspecified, cox-TMLE (c) = "*lambda^{'c'}*","*lambda^{'y'}*" misspecified")))


ggsave(file = paste0("./figures/",
                     "fig-simulation-results-7A-2025",
                     ".pdf"),
       p.7A, width = 26, height = 8, units = "cm")


ggsave(file = paste0("./figures/",
                     "fig-simulation-results-4A-7A-2025",
                     ".pdf"),
       grid.arrange(p.4A,
                    p.7A,
                    nrow = 2), width = 26, height = 8*2, units = "cm")


######################################################################
### 04_make_figures.R ends here
