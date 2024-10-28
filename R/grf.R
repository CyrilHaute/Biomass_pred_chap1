formula = fmla
dframe = training
bw = 100
kernel = "fixed"
coords = coords
ntree = nrow(training) * 10
geo.weighted = FALSE
mtry = NULL
importance = "impurity"
nthreads = NULL
forests = TRUE
geo.weighted = FALSE
print.results = TRUE

grf <- function(formula, 
                dframe, 
                bw, 
                kernel, 
                coords, 
                ntree = 500, 
                mtry = NULL, 
                importance = "impurity", 
                nthreads = NULL, 
                forests = TRUE, 
                geo.weighted = TRUE,  
                print.results = TRUE){
  
  # Start timing the function execution
  start.time <- Sys.time()
  
  # Convert formula text to a formula object
  f <- formula(formula)
  
  # Extract variable names from the formula
  RNames <- attr(terms(f), "term.labels")
  
  # Get the name of the dependent variable
  DepVarName <- row.names(attr(terms(f), "factors"))[1]
  
  # Create a data frame for the dependent variable
  Y.DF <- dframe[DepVarName]
  
  # Convert the dependent variable data frame to a vector
  Y <- as.numeric(as.character(Y.DF[[1]]))
  
  # Determine the number of independent variables and add 1 for degrees of freedom
  ModelVarNo <- length(RNames)
  K = ModelVarNo + 1
  
  # Set the number of trees in the model
  ntrees <- ntree
  
  # Count the number of observations in the data
  Obs <- nrow(dframe)
  
  # Define mtry if it is not provided [max(floor(Number of Variables/3), 1)]
  if (is.null(mtry)) {mtry= max(floor(ModelVarNo/3), 1)}
  
  # Print initial information if required
  if(print.results) {
    message("\nNumber of Observations: ", Obs)
    message("Number of Independent Variables: ", ModelVarNo)
  }
  
  # Configure the kernel type and its parameters
  if(kernel == 'adaptive'){
    
    Ne <- bw
    if(print.results) {message("Kernel: Adaptive\nNeightbours: ", Ne)}
    
  }else{
    
    if(kernel == 'fixed'){
      
      if(print.results) {
        
        message("Kernel: Fixed\nBandwidth: ", bw)
        
        }
    }
    
  }
  
  print("Fit global random forest")
  # Fit the global random forest model using the ranger package
  Gl.Model <- eval(substitute(ranger::ranger(formula, data = dframe, num.trees=ntree, mtry= mtry, importance=importance, num.threads = nthreads)))
  
  # Get predictions from the global model
  Predict <- predict(Gl.Model, dframe, num.threads = nthreads)
  
  yhat <- as.numeric(as.character(Predict$predictions))
  
  # Print global model summary if required
  if(print.results) {
    message("\n--------------- Global ML Model Summary ---------------\n")
    print(Gl.Model)
    
    message("\nImportance:\n")
    print(Gl.Model$variable.importance)
    
    #calculate pseudoR2
    g.RSS <- sum((Y-yhat)^2)
    g.mean.y <- mean(Y)
    g.TSS<-sum((Y-g.mean.y)^2)
    
    g.r<-1-(g.RSS/g.TSS)
    
    g.AIC <- 2*K + Obs*log(g.RSS/Obs)
    
    g.AICc <- g.AIC + ((2*K*(K +1)) / (Obs - K - 1))
    
    message("\nMean Square Error (Not OOB): ", round(g.RSS/Obs,3))
    message("R-squared (Not OOB) %: ", round(100 * g.r,3))
    message("AIC (Not OOB): ", round(g.AIC,3))
    message("AICc (Not OOB): ", round(g.AICc,3))
  }
  
  # Calculate distances between observations based on coordinates

  DistanceT <- dist(coords)
  Dij <- as.matrix(DistanceT)

  print("Fit and save local random forest")
  fit <- pbmcapply::pbmclapply(1:Obs, function(m) {
    
    DNeighbour <- Dij[,m]
  
    #Get the data
    DataSet <- data.frame(dframe, DNeighbour = DNeighbour)
    
    #Sort by distance
    DataSetSorted <- DataSet[order(DataSet$DNeighbour),]
    
    if(kernel == 'adaptive'){
      
      #Keep Nearest Neighbours
      SubSet <- DataSetSorted[1:Ne,]
      Kernel_H <- max(SubSet$DNeighbour)
      
    }else{
      
      if(kernel == 'fixed'){
        
        SubSet <- subset(DataSetSorted, DNeighbour <= bw)
        Kernel_H <- bw
        
      }
      
    }

    #Bi-square weights
    Wts <- (1-(SubSet$DNeighbour/Kernel_H)^2)^2
    
    #Calculate WLM
    if (geo.weighted == TRUE) {
      Lcl.Model <- eval(substitute(ranger::ranger(formula, data = SubSet, num.trees=ntree, mtry= mtry, importance=importance, case.weights=Wts, num.threads = nthreads)))
      
      local.predicted.y <- Lcl.Model$predictions[[1]]
      counter <- 1
      while (is.nan(local.predicted.y)) {
        Lcl.Model<-eval(substitute(ranger::ranger(formula, data = SubSet, num.trees=ntree, mtry= mtry, importance=importance, case.weights=Wts, num.threads = nthreads)))
        local.predicted.y <- Lcl.Model$predictions[[1]]
        counter <- counter + 1
      }
    } else{
      Lcl.Model<-eval(substitute(ranger::ranger(formula, data = SubSet, num.trees=ntree, mtry= mtry, importance=importance, num.threads = nthreads)))
      counter <- 1
    }
    
    
    
    if (forests == TRUE) {LM_Forests <- Lcl.Model}
    
    LM_LEst <- Lcl.Model$variable.importance
    
    #Observed y
    l.predict <- predict(Lcl.Model, dframe[m,], num.threads = nthreads)
    
    LM_GofFit <- data.frame(y=Y[m], 
                            LM_yfitOOB=Lcl.Model$predictions[[1]], 
                            LM_ResOOB=as.numeric(as.character(Y[m])) - as.numeric(as.character(Lcl.Model$predictions[[1]])), 
                            LM_yfitPred=l.predict$predictions,
                            LM_ResPred=as.numeric(as.character(Y[m])) - as.numeric(as.character(l.predict$predictions)), 
                            LM_MSE=Lcl.Model$prediction.error, 
                            LPerm=counter)
    
    # rm(l.predict, Lcl.Model, Wts, SubSet, Kernel_H, DataSetSorted, DataSet, DNeighbour)
    # gc()

    if (forests == TRUE) {
      
      # saveRDS(LM_Forests, file = paste0("results/local_rf/local_random_forest", "_", m, ".rds"))
      return(list(LM_GofFit, LM_LEst, LM_Forests))

    }else{
      
      return(list(LM_GofFit, LM_LEst))
      
    }

  }, mc.cores = parallel::detectCores() - 1)
  
  print("Get local variable importance and local random forest performance")

  Global.Model <- Gl.Model
  Locations <- coordinates
  Local.Variable.Importance <- lapply(fit, '[[', 2)
  Local.Variable.Importance <- do.call(rbind, Local.Variable.Importance)
  LGofFit <- lapply(fit, '[[', 1)
  LGofFit <- do.call(rbind, LGofFit)
  # Forests <- lapply(local_rf, '[[', 3)
  
  # grf.out <- list(Global.Model = Global.Model, Locations = Locations, Local.Variable.Importance = Local.Variable.Importance, LGofFit = LGofFit, Forests = Forests)
  grf.out <- list(Global.Model = Global.Model, Locations = Locations, Local.Variable.Importance = Local.Variable.Importance, LGofFit = LGofFit)

  if(print.results) {

    message("\n--------------- Local Model Summary ---------------\n")

    message("\nResiduals OOB:\n")
    print(summary(grf.out$LGofFit$LM_ResOOB))

    message("\nResiduals Predicted (Not OOB):\n")

    print(summary(grf.out$LGofFit$LM_ResPred))

  }
  lvi <- data.frame(Min = apply(grf.out$Local.Variable.Importance, 2, min), Max = apply(grf.out$Local.Variable.Importance, 2, max),
                    Mean = apply(grf.out$Local.Variable.Importance, 2, mean), StD = apply(grf.out$Local.Variable.Importance, 2, sd))


  l.RSS.OOB <- sum(grf.out$LGofFit$LM_ResOOB^2)
  l.RSS.Pred<-sum(grf.out$LGofFit$LM_ResPred^2)

  mean.y<-mean(as.numeric(as.character(grf.out$LGofFit$y)))
  TSS<-sum((as.numeric(as.character(grf.out$LGofFit$y))-mean.y)^2)

  l.r.OOB<-1-(l.RSS.OOB/TSS)
  g.AIC.OOB <- 2*K + Obs*log(l.RSS.OOB/Obs)
  g.AICc.OOB <- g.AIC.OOB + ((2*K*(K + 1)) / (Obs - K - 1))



  l.r.Pred<-1-(l.RSS.Pred/TSS)
  g.AIC.Pred <- 2*K + Obs*log(l.RSS.Pred/Obs)
  g.AICc.Pred <- g.AIC.Pred + ((2*K*(K +1)) / (Obs - K - 1))

  if(print.results) {

    message("\nLocal Variable Importance:\n")
    print(lvi)
    message("\nMean squared error (OOB): ", round(l.RSS.OOB/Obs,3))
    message("R-squared (OOB) %: ", round(100* l.r.OOB,3))
    message("AIC (OOB): ", round(g.AIC.OOB,3))
    message("AICc (OOB): ", round(g.AICc.OOB,3))
    message("Mean squared error Predicted (Not OOB): ", round(l.RSS.Pred/Obs,3))
    message("R-squared Predicted (Not OOB) %: ", round(100* l.r.Pred,3))
    message("AIC Predicted (Not OOB): ", round(g.AIC.Pred,3))
    message("AICc Predicted (Not OOB): ", round(g.AICc.Pred,3))
  }

  lModelSummary = list()
  lModelSummary$l.VariableImportance <- lvi
  lModelSummary$l.MSE.OOB <- l.RSS.OOB/Obs
  lModelSummary$l.r.OOB <- l.r.OOB
  lModelSummary$l.MSE.Pred <- l.RSS.Pred/Obs
  lModelSummary$l.r.Pred <- l.r.Pred

  grf.out$LocalModelSummary <- lModelSummary

  # Calculate and print the time taken to run the function
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  
  if(print.results) {message("\nCalculation time (in seconds): ", round(time.taken,4))}
  
  # Return the output list
  return(grf.out)
}

# object = model_fit
# new.data = test
# x.var.name = "X"
# y.var.name = "Y"
# local.w = 0.5
# global.w = 0.5


predict.grf <- function(object, new.data, x.var.name, y.var.name, local.w=1, global.w=0, ...) {
  
  Obs <- nrow(new.data)
  
  file_local_rf <- list.files("results/local_rf", full.names = TRUE)
  locations <- object$Locations
  
  # predictions <- vector(mode="numeric", length=Obs)
  
  predictions <- pbmcapply::pbmclapply(1:Obs, function(i) {
    
    x <- new.data[i, which(names(new.data)==x.var.name)]
    y <- new.data[i, which(names(new.data)==y.var.name)]
    
    D <- sqrt((x-locations[,1])^2 + (y-locations[,2])^2)
    
    local.model.ID <- which.min(D)

    local.model <- readRDS(file_local_rf[[local.model.ID]])
    
    g.predict <- predict(object[[1]], new.data[i,])
    g.prediction <- g.predict$predictions
    l.predict <- predict(local.model, new.data[i,])
    l.prediction <- l.predict$predictions
    
    predictions <- global.w * as.numeric(as.character(g.prediction[1])) + local.w * as.numeric(as.character(l.prediction[1]))
    
  }, mc.cores = parallel::detectCores() - 1)
  
  predictions <- unlist(predictions)
  
  return(predictions)
}
