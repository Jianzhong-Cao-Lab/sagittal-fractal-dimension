eval4sa <- function(
    predprob, preddata, etime, estatus, model, dataset, timepoints,
    plotcalimethod, bw4nne, q4quantile
) {
  survobj <- as.formula(paste0("Surv(", etime, ", ", estatus, ") ~ 1"))
  aucdf <- list()
  rocdf <- list()
  bsdf <- list()
  bscdf <- list()
  caldf <- list()

  for (i in seq_along(timepoints)) {
    datai <- data.frame(vi = 1 - predprob[, i])
    colnames(datai) <- model
    set.seed(511543)
    score_obj <- riskRegression::Score(
      as.list(datai),
      formula = survobj,
      data = preddata,
      metrics = c("auc", "brier"),
      summary = "IPA",
      plots = c("roc", "calibrate"),
      times = timepoints[i],
      conf.int = TRUE
    )
    aucdf[[i]] <- score_obj$AUC$score %>% mutate(dataset = dataset)
    rocdf[[i]] <- score_obj$ROC$plotframe %>% mutate(dataset = dataset)
    bsdf[[i]] <- score_obj$Brier$score %>% mutate(dataset = dataset)
    bscdf[[i]] <- score_obj$Brier$contrasts %>% mutate(dataset = dataset)
    caldf[[i]] <- riskRegression::plotCalibration(
      score_obj,
      cens.method = "local",
      method = plotcalimethod,
      bandwidth = bw4nne,
      q = q4quantile,
      plot = FALSE
    )$plotFrames[[1]] %>%
      mutate(dataset = dataset, times = timepoints[i], model = model)
  }

  dataauc <- bind_rows(aucdf)
  dataroc <- bind_rows(
    rocdf,
    data.frame(
      model = model,
      times = timepoints,
      risk = -Inf,
      TPR = 0,
      FPR = 0,
      dataset = dataset
    )
  ) %>% arrange(TPR)
  databs <- bind_rows(bsdf)
  databsc <- bind_rows(bscdf)
  datacal <- bind_rows(caldf)
  rownames(datacal) <- NULL

  plotroc <- dataroc %>%
    left_join(dataauc, by = c("model", "times", "dataset")) %>%
    mutate(
      tauc = paste0(
        "T=", times, ", AUC=", round(AUC, 3),
        "(", round(lower, 3), "~", round(upper, 3), ")"
      ),
      tauc = forcats::as_factor(tauc)
    ) %>%
    ggplot(aes(x = FPR, y = TPR, group = tauc, color = tauc)) +
    geom_line(linewidth = 1) +
    geom_abline(color = "grey") +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = "1-Specificity", y = "Sensitivity", color = "") +
    theme_bw() +
    theme(
      legend.position = c(1, 0),
      legend.justification = c(1, 0),
      legend.background = element_blank(),
      legend.key = element_blank()
    )

  plotcal <- datacal %>%
    mutate(t = forcats::as_factor(paste0("T=", times))) %>%
    ggplot(aes(x = Pred, y = Obs, group = t, color = t)) +
    geom_line(linewidth = 1) +
    geom_point(size = 3, alpha = ifelse(plotcalimethod == "nne", 0, 1)) +
    geom_abline() +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    labs(x = "Predicted risk", y = "Estimated actual risk", color = "") +
    theme_bw() +
    theme(
      legend.position = c(1, 0),
      legend.justification = c(1, 0),
      legend.background = element_blank(),
      legend.key = element_blank()
    )

  list(
    auc = dataauc,
    roc = dataroc,
    rocplot = plotroc,
    brierscore = databs,
    brierscoretest = databsc,
    calibration = datacal,
    calibrationplot = plotcal
  )
}

risk_pred <- function(model, newdata) {
  model$predict_newdata(newdata)$crank
}

surv_pred <- function(model, newdata, times) {
  t(model$predict_newdata(newdata)$distr$survival(times))
}

chf_pred <- function(model, newdata, times) {
  t(model$predict_newdata(newdata)$distr$cumHazard(times))
}
