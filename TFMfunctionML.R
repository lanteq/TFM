require(KnowSeq)
require(caret)
require(e1071)
require(stringr)
require(arulesCBA) #discretizar
require(FSelector) #Ganancia de información
require(lmtest)

StratifiedCV_AnalysisFS <- function(X, Y, ngenes, fs, sub_fs="backward", clasif, L, lfc = 2, cov = 2, balanced=F, after=F,down=F){
  
  vectorOfF1score <- double(length = ngenes-1)
  prediction <- matrix(nrow = ngenes-1, ncol = dim(X)[1])
  first10genes_perfold <- list()
  confusionMatrixFold <- list()
  FoldValidacion<-createFolds(Y, k = L)
  j <- 1
  for(l in FoldValidacion){ # Para cada segmento
    XValCV <- X[l,]
    XTrnCV <- X[-l,]
    YTrnCV <- Y[-l]
    
    # Balanceado de clases pre seleccion
    if(balanced==T && after == F && down == T){
      down_train <- downSample(x = XTrnCV,y = YTrnCV) 
      XTrnCV <- down_train[,-ncol(down_train)]
      YTrnCV <- down_train$Class
    }
    else if(balanced==T && after == F && down == F){
      up_train <- upSample(x = XTrnCV,y = YTrnCV)
      XTrnCV <- up_train[,-ncol(up_train)]
      YTrnCV <- up_train$Class
    }
    
    DEGsInfo <- DEGsExtraction(t(XTrnCV), YTrnCV, lfc = lfc, cov = cov, pvalue = 0.001)
    DEGsMatrix <- DEGsInfo$DEG_Results$DEGs_Matrix
    MLMatrix <- t(DEGsMatrix)
    
    
    # Selección de características  
    if(fs=="mrmr"){
      FSRanking <- featureSelection(MLMatrix, YTrnCV, mode = "mrmr", vars_selected = colnames(MLMatrix))
      first10genes_perfold[[j]] <- names(FSRanking[1:ngenes])
      datosTrain <- as.data.frame(MLMatrix[,names(FSRanking[1:ngenes])])
      XTest <- XValCV[,names(FSRanking[1:ngenes])]
    }
    else if(fs=="rf"){
      FSRanking <- featureSelection(MLMatrix, YTrnCV, mode = "rf", vars_selected = colnames(MLMatrix))
      first10genes_perfold[[j]] <- FSRanking[1:ngenes]
      datosTrain <- as.data.frame(MLMatrix[,FSRanking[1:ngenes]])
      XTest <- XValCV[,FSRanking[1:ngenes]]
    }
    else if(fs=="semanticRL"){
      GOsList <- geneOntologyEnrichment(colnames(MLMatrix), geneType = "GENE_SYMBOL", ontologies = "BP", pvalCutOff=0.6)
      GENESperGO<-GOsList$`BP Ontology GOs`$`Gene Symbols`
      GENESperGO_split<-str_split(GENESperGO, ", ")
      
      datos <- cbind(MLMatrix,YTrnCV)
      datos <- data.frame(datos)
      datos$YTrnCV <- as.factor(datos$YTrnCV)
  
      
      GENESperGO_split <- unique(GENESperGO_split)
      Num.GOs <- length(GENESperGO_split)
      if(sub_fs=="backward") {
        for(i in 1:Num.GOs){
          pvalores <- c()
          while (length(GENESperGO_split[[i]])>2 | max(pvalores)>0.0001) {
            modelFULL<-list()
            list_compared<-list()
            model_reduced<-list()
            myformula <- as.formula(paste("YTrnCV ~ ", paste(GENESperGO_split[[i]], collapse = " + ")))
            modelFULL <- nnet::multinom(myformula, data = datos, trace = FALSE)
            model_reduced<-lapply(GENESperGO_split[[i]], function(j) update(modelFULL,as.formula(paste(".~.-", j))))
            list_compared<-lapply(model_reduced,function(j) lrtest(modelFULL,j))
            pvalores <- sapply(list_compared, "[[",2,5) #pvalores de los distintos test enfrentados
            pos <- which.max(pvalores)
            GENESperGO_split[[i]]<- GENESperGO_split[[i]][-pos]
          }
        }
        genes_selected0 <- unique(unlist(GENESperGO_split))
      } else {
        selectedGeneperGO <-vector("list", length = Num.GOs)
        for(i in 1:Num.GOs){
          model_one <- nnet::multinom(YTrnCV ~ 1, data = datos, trace = FALSE)
          pvalores <- 0
          s<-0
          while (s<2 & min(pvalores)<0.0001 & length(GENESperGO_split[[i]])>0){
            list_compared<-list()
            model_reduced<-list()
            model_reduced<-lapply(GENESperGO_split[[i]], function(j) update(model_one,as.formula(paste(".~.+", j))))
            list_compared<-lapply(model_reduced,function(j) lrtest(j,model_one))
            pvalores <- sapply(list_compared, "[[",2,5)
            pos <- which.min(pvalores)
            selectedGeneperGO[[i]] <- GENESperGO_split[[i]][pos]
            GENESperGO_split[[i]]<- GENESperGO_split[[i]][-pos]
            model_one <- model_reduced[[pos]]
            s<-s+1
          }
        }
        genes_selected0 <- unique(unlist(selectedGeneperGO))
      }
      genes_selected <- c()
      model_one <- nnet::multinom(YTrnCV ~1, data = datos, trace = FALSE)
      for (i in 1:ngenes) {
        model_reduced<-lapply(genes_selected0, function(j) update(model_one,as.formula(paste(".~.+", j))))
        test_compared <- lapply(model_reduced,function(j) lrtest(j,model_one))
        pvalores <- sapply(test_compared, "[[",2,5)
        pos <- which.min(pvalores)
        genes_selected[i] <- genes_selected0[pos]
        genes_selected0 <- genes_selected0[-pos]
        model_one <- model_reduced[[pos]]
      }
      first10genes_perfold[[j]] <- genes_selected
      datosTrain <- as.data.frame(MLMatrix[,genes_selected])
      XTest <- XValCV[,genes_selected]
    }
    else if(fs=="semantic"){
      datos <- cbind(MLMatrix,YTrnCV)
      datos <- as.data.frame(datos)
      datos$YTrnCV <- as.factor(datos$YTrnCV)
      data_disc <- discretizeDF.supervised(YTrnCV ~ ., datos)
      
      ig <- information.gain(YTrnCV~., data_disc)
      listgenes <- gsub("\\.","", rownames(ig))
      listgenes <- gsub("-","", listgenes)
      
      GOsList <- geneOntologyEnrichment(listgenes, geneType = "GENE_SYMBOL", ontologies = "BP", pvalCutOff=0.6)
      GENESperGO<-GOsList$`BP Ontology GOs`$`Gene Symbols`
      GENESperGO_split<-str_split(GENESperGO, ", ")
      Num.GOs <- length(GENESperGO_split)
      GO_list <- list()
      GOs <- list()
      GOs_scores <- double(length = Num.GOs)
      for(i in seq(1:Num.GOs)){
        pos <- match(GENESperGO_split[[i]],listgenes)
        GO_list <- NULL
        GO_list[GENESperGO_split[[i]]]<-ig[pos,]
        GOs[[i]] <- GO_list
        GOs_scores[i] <- mean(unlist(GOs[[i]]))
      }
      
      genes_selected <- c()
      while (length(genes_selected)<ngenes){
        
        K <- which.max(GOs_scores)
        genes_selected <- append(genes_selected,names(which.max(unlist(GOs[[K]]))))
        
        for (i in seq(1:Num.GOs)){
          GOs[[i]] <- GOs[[i]][setdiff(names(GOs[[i]]),genes_selected)]
          GOs_scores[i] <- mean(unlist(GOs[[i]]))
        }
      }
      first10genes_perfold[[j]] <- genes_selected
      datosTrain <- as.data.frame(MLMatrix[,genes_selected])
      XTest <- XValCV[,genes_selected]
    }
    
    # Balanceado de clases post seleccion
    if(balanced==T && after == T && down == T){
      down_train <- downSample(x = datosTrain,y = YTrnCV)
      datosTrain <- down_train[,-ncol(down_train)]
      YTrnCV <- down_train$Class
    }
    else if(balanced==T && after == T && down == F){
      up_train <- upSample(x = datosTrain,y = YTrnCV)
      datosTrain <- up_train[,-ncol(up_train)]
      YTrnCV <- up_train$Class
    }
    
    #Normalizacion
    maximum <-apply(datosTrain, 2, max)
    minimum <- apply(datosTrain, 2, min)
    for(i in seq(1:ngenes)){
      datosTrain[,i] <- ((datosTrain[,i] - minimum[i]) / (maximum[i] - minimum[i])) * 2 - 1
    }
    for(i in seq(1:ngenes)){
      XTest[,i] <- ((XTest[,i] - minimum[i]) / (maximum[i] - minimum[i])) * 2 - 1
    }
    
    
    if(clasif=="knn"){
      K<-15
      vectorOfAccsCV <- double(length = K)
      for (i in seq(1:(ngenes-1))) {
        for(k in seq(1:K)){ # Para cada valor de k
          knn_mod_CV = knn3(datosTrain[,1:(i+1)], as.factor(YTrnCV), k = k)
          values <- predict(knn_mod_CV, XTest[,1:(i+1)], type = "class")
          predicts <- as.double(as.character(values))
          confMatrix <- confusionMatrix(as.factor(predicts), Y[l])
          vectorOfAccsCV[k] <- confMatrix$overall[1]
        }
        BestK <- which.max(vectorOfAccsCV)
        BestAcc <- max(vectorOfAccsCV)
        knn_mod = knn3(datosTrain[,1:(i+1)], as.factor(YTrnCV), k = BestK)
        prediction[i,l] <- predict(knn_mod, XTest[,1:(i+1)], type = "class")
      }
      
    }
    
    else if(clasif=="svm"){
      fitControl <- trainControl(method = "cv", number = 10)
      grid_radial <- expand.grid(sigma = c(0,0.01, 0.02, 0.025, 0.03, 0.04,0.05, 0.06, 0.07,0.08, 0.09, 0.1, 0.25, 0.5, 0.75,0.9), C = c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 1, 1.5, 2,5))
      dataForTunning <- cbind(datosTrain,YTrnCV)
      for (i in seq(1:(ngenes-1))) {
        Rsvm_sb <- train(YTrnCV ~ ., data = dataForTunning[,c(1:(i+1),dim(dataForTunning)[2])],type = "C-svc", method = "svmRadial", trControl = fitControl,tuneGrid = grid_radial)
        bestC <- Rsvm_sb$bestTune$C
        bestG <- Rsvm_sb$bestTune$sigma
        svm_model <- svm(datosTrain[,1:(i+1)],as.factor(YTrnCV),kernel='radial', cost=bestC,gamma=bestG,probability=TRUE)
        prediction[i,l] <- as.double(as.character(predict(svm_model,XTest[,1:(i+1)],probability=TRUE)))
      }
      
    }
    confusionMatrixFold[[j]] <- confusionMatrix(as.factor(prediction[ngenes-1,l]),Y[l])
    j <- j+1
  }
  
  #Analizamos  el rendimiento de nuestra clasificación:
  for (i in seq(1:(ngenes-1))){
    confMatrix <- confusionMatrix(as.factor(prediction[i,]), Y)
    suma = apply(confMatrix[[4]][,c(1,3)], 1, sum)
    producto = apply(confMatrix[[4]][,c(1,3)], 1, prod)
    vectorofF1 = (2*producto)/suma
    vectorOfF1score[i] = sum(vectorofF1)/length(vectorofF1)
  }
  return(list(F1score=vectorOfF1score,genes_Selected=first10genes_perfold, datos=datosTrain, predicciones=prediction, reY=YTrnCV, confMatrixs=confusionMatrixFold, mlmatrix=MLMatrix))
}



