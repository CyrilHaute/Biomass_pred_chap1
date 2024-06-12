# formula = fmla
# dframe = training
# bw = 100
# kernel = "fixed"
# coords = coords
# ntree = 100
# mtry = NULL
# importance = "impurity"
# nthreads = NULL
# forests = TRUE
# geo.weighted = FALSE
# print.results = TRUE
# path_forest = path_forest

grf <- function (formula, dframe, bw, kernel, coords, ntree = 500, 
                 mtry = NULL, importance = "impurity", nthreads = NULL, forests = TRUE, 
                 geo.weighted = TRUE, print.results = TRUE, path_forest) 
{
  start.time <- Sys.time()
  f <- formula(formula)
  RNames <- attr(terms(f), "term.labels")
  DepVarName <- row.names(attr(terms(f), "factors"))[1]
  Y.DF <- dframe[DepVarName]
  Y <- Y.DF[[1]]
  ModelVarNo <- length(RNames)
  K = ModelVarNo + 1
  ntrees <- ntree
  Obs <- nrow(dframe)
  if (is.null(mtry)) {
    mtry = max(floor(ModelVarNo/3), 1)
  }
  if (print.results) {
    message("\nNumber of Observations: ", Obs)
    message("Number of Independent Variables: ", ModelVarNo)
  }
  if (kernel == "adaptive") {
    Ne <- bw
    if (print.results) {
      message("Kernel: Adaptive\nNeightbours: ", Ne)
    }
  }else {
    if (kernel == "fixed") {
      if (print.results) {
        message("Kernel: Fixed\nBandwidth: ", bw)
      }
    }
  }
  Gl.Model <- eval(substitute(ranger::ranger(formula, data = dframe, 
                                             num.trees = ntree, mtry = mtry, importance = importance, 
                                             num.threads = nthreads)))
  Predict <- predict(Gl.Model, dframe, num.threads = nthreads)
  yhat <- Predict$predictions
  if (print.results) {
    message("\n--------------- Global ML Model Summary ---------------\n")
    print(Gl.Model)
    message("\nImportance:\n")
    print(Gl.Model$variable.importance)
    g.RSS <- sum((Y - yhat)^2)
    g.mean.y <- mean(Y)
    g.TSS <- sum((Y - g.mean.y)^2)
    g.r <- 1 - (g.RSS/g.TSS)
    g.AIC <- 2 * K + Obs * log(g.RSS/Obs)
    g.AICc <- g.AIC + ((2 * K * (K + 1))/(Obs - K - 1))
    message("\nMean Square Error (Not OOB): ", round(g.RSS/Obs, 
                                                     3))
    message("R-squared (Not OOB) %: ", round(100 * g.r, 
                                             3))
    message("AIC (Not OOB): ", round(g.AIC, 3))
    message("AICc (Not OOB): ", round(g.AICc, 3))
  }
  DistanceT <- dist(coords)
  Dij <- as.matrix(DistanceT)

  local_pred <- pbmcapply::pbmclapply(1:Obs, function(m) {
    # for (m in 1:Obs) {
    DNeighbour <- Dij[, m]
    DataSet <- data.frame(dframe, DNeighbour = DNeighbour)
    DataSetSorted <- DataSet[order(DataSet$DNeighbour),]
    if (kernel == "adaptive") {
      SubSet <- DataSetSorted[1:Ne, ]
      Kernel_H <- max(SubSet$DNeighbour)
    }else {
      if (kernel == "fixed") {
        SubSet <- subset(DataSetSorted, DNeighbour <= 
                           bw)
        Kernel_H <- bw
      }
    }
    Wts <- (1 - (SubSet$DNeighbour/Kernel_H)^2)^2
    if (geo.weighted == TRUE) {
      Lcl.Model <- eval(substitute(ranger(formula, data = SubSet, 
                                          num.trees = ntree, mtry = mtry, importance = importance, 
                                          case.weights = Wts, num.threads = nthreads, 
                                          ...)))
      local.predicted.y <- Lcl.Model$predictions[[1]]
      counter <- 1
      while (is.nan(local.predicted.y)) {
        Lcl.Model <- eval(substitute(ranger(formula, 
                                            data = SubSet, num.trees = ntree, mtry = mtry, 
                                            importance = importance, case.weights = Wts, 
                                            num.threads = nthreads, ...)))
        local.predicted.y <- Lcl.Model$predictions[[1]]
        counter <- counter + 1
      }
    }else {
      Lcl.Model <- eval(substitute(ranger::ranger(formula, data = SubSet, 
                                                  num.trees = ntree, mtry = mtry, importance = importance, 
                                                  num.threads = nthreads)))
      counter <- 1
    }
    if (forests == TRUE) {
      save(Lcl.Model, file = paste0(path_forest, "/", m, ".Rdata"))
      # LM_Forests <- Lcl.Model
    }
    LM_LEst <- Lcl.Model$variable.importance
    
    l.predict <- predict(Lcl.Model, dframe[m, ], num.threads = nthreads)
    
    LM_GofFit <- data.frame(y=Y[m], 
                            LM_yfitOOB=Lcl.Model$predictions[[1]], 
                            LM_ResOOB=Y[m] - Lcl.Model$predictions[[1]], 
                            LM_yfitPred=l.predict$predictions,
                            LM_ResPred=Y[m] - l.predict$predictions, 
                            LM_MSE=Lcl.Model$prediction.error, 
                            LM_R2 = Lcl.Model$r.squared,
                            LPerm=counter)
    
    return(list(LM_LEst, LM_GofFit))
    
  }, mc.cores = parallel::detectCores() - 1)
  
  if (forests == TRUE) {
    grf.out <- list(Global.Model = Gl.Model, Locations = coords, 
                    Local.Variable.Importance = do.call(rbind, lapply(local_pred, "[[", 1)), LGofFit = do.call(rbind, lapply(local_pred, "[[", 2)), 
                    Forests = path_forest)
  } else {
    grf.out <- list(Global.Model = Gl.Model, Locations = coords, 
                    Local.Variable.Importance = LM_LEst, LGofFit = LM_GofFit)
  }
  if (print.results) {
    message("\n--------------- Local Model Summary ---------------\n")
    message("\nResiduals OOB:\n")
    print(summary(grf.out$LGofFit$LM_ResOOB))
    message("\nResiduals Predicted (Not OOB):\n")
    print(summary(grf.out$LGofFit$LM_ResPred))
  }
  lvi <- data.frame(Min = apply(grf.out$Local.Variable.Importance, 
                                2, min), Max = apply(grf.out$Local.Variable.Importance, 
                                                     2, max), Mean = apply(grf.out$Local.Variable.Importance, 
                                                                           2, mean), StD = apply(grf.out$Local.Variable.Importance, 
                                                                                                 2, sd))
  l.RSS.OOB <- sum(grf.out$LGofFit$LM_ResOOB^2, na.rm = TRUE)
  l.RSS.Pred <- sum(grf.out$LGofFit$LM_ResPred^2, na.rm = TRUE)
  mean.y <- mean(grf.out$LGofFit$y, na.rm = TRUE)
  TSS <- sum((grf.out$LGofFit$y - mean.y)^2, na.rm = TRUE)
  l.r.OOB <- 1 - (l.RSS.OOB/TSS)
  g.AIC.OOB <- 2 * K + Obs * log(l.RSS.OOB/Obs)
  g.AICc.OOB <- g.AIC.OOB + ((2 * K * (K + 1))/(Obs - K - 
                                                  1))
  l.r.Pred <- 1 - (l.RSS.Pred/TSS)
  g.AIC.Pred <- 2 * K + Obs * log(l.RSS.Pred/Obs)
  g.AICc.Pred <- g.AIC.Pred + ((2 * K * (K + 1))/(Obs - K - 
                                                    1))
  if (print.results) {
    message("\nLocal Variable Importance:\n")
    print(lvi)
    message("\nMean squared error (OOB): ", round(l.RSS.OOB/Obs, 
                                                  3))
    message("R-squared (OOB) %: ", round(100 * l.r.OOB, 
                                         3))
    message("AIC (OOB): ", round(g.AIC.OOB, 3))
    message("AICc (OOB): ", round(g.AICc.OOB, 3))
    message("Mean squared error Predicted (Not OOB): ", 
            round(l.RSS.Pred/Obs, 3))
    message("R-squared Predicted (Not OOB) %: ", round(100 * 
                                                         l.r.Pred, 3))
    message("AIC Predicted (Not OOB): ", round(g.AIC.Pred, 
                                               3))
    message("AICc Predicted (Not OOB): ", round(g.AICc.Pred, 
                                                3))
  }
  lModelSummary = list()
  lModelSummary$l.VariableImportance <- lvi
  lModelSummary$l.MSE.OOB <- l.RSS.OOB/Obs
  lModelSummary$l.r.OOB <- l.r.OOB
  lModelSummary$l.MSE.Pred <- l.RSS.Pred/Obs
  lModelSummary$l.r.Pred <- l.r.Pred
  grf.out$LocalModelSummary <- lModelSummary
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  if (print.results) {
    message("\nCalculation time (in seconds): ", round(time.taken, 
                                                       4))
  }
  return(grf.out)
}

# object = model_fit
# new.data = as.data.frame(testing)
# x.var.name = "X"
# y.var.name = "Y"
# local.w = 0.5
# global.w = 0.5


predict.grf <- function(object, 
                        new.data,
                        x.var.name,
                        y.var.name,
                        local.w = 1,
                        global.w = 0){
  
  Obs <- nrow(new.data)
  
  list_forest <- list.files(object$Forests, full.names = TRUE)
  
  pred <- unlist(pbmcapply::pbmclapply(1:Obs, function(i) {
    
    x <- new.data[i, which(names(new.data) == x.var.name)]
    
    y <- new.data[i, which(names(new.data) == y.var.name)]
    
    locations <- as.data.frame(object$Locations)
    
    D <- sqrt((x - locations[, 1])^2 + (y - locations[,2])^2)
    
    local.model.ID <- which.min(D)
    
    g.predict <- predict(object[[1]], new.data[i, ])
    
    g.prediction <- g.predict$predictions
    
    load(paste0("outputs/biomass_prediction/sprf_local/", local.model.ID, ".Rdata"))
    
    l.predict <- predict(Lcl.Model, new.data[i, ])
    
    l.prediction <- l.predict$predictions
    
    predictions <- global.w * g.prediction[1] + local.w * l.prediction[1]
    
  }, mc.cores = parallel::detectCores() - 1))
  
  return(pred)
}
