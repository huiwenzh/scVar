#'Identifying biological variability in scRNA-seq data
#'@author Huiwen Zheng
#'@param dat  Raw count matrix
#'@param norm Logic parameter to indicate normalisation (default to F)
#'@param metric Name of the metric to run for analysis
#'@param org Organisms name, only required for DM
#'@param group covariates for the count matrix.
#'@return A named vector with corresponding gene expression variability 
#'
#'@export
scVar <- function(dat, norm=F, metric=c('SD',"MAD","IQR","CV","FF",'scran','LCV','DESeq2','edgeR','Seurat_mvp','Seurat_vst',  'DM',"glmGamPoi","BASiCS"), org='mm9', group=NULL){
  if (ncol(dat)<50){
    message('Small sample size may lead to incorrect estimation')
  }
  dat <- dat[Matrix::rowSums(dat)!=0,]
  if (norm){
    library(scran)
    library(SingleCellExperiment)
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = dat),  colData=DataFrame(label=group))
    sce <- scran::computeSumFactors(sce)
    sce <- scater::logNormCounts(sce)
    dat <- SingleCellExperiment::logcounts(sce)
    if (metric =='SD'){
      vars <- apply(dat,1,function(x) sd(x))
    }
    else if (metric == 'MAD'){
      vars <- apply(dat,1,function(x) mad(x))
    }
    else if (metric == 'IQR'){
      vars <- apply(dat,1,function(x) IQR(x))
    }
    else if (metric == 'CV'){
      means <- rowMeans(as.matrix(dat))
      vars <- apply(dat,1,function(x) sd(x))/means
    }
    else if (metric == 'FF'){
      means <- rowMeans(as.matrix(dat))
      vars <- apply(dat,1,function(x) var(x))/means
    }
    else if (metric == 'LCV'){
      CVs <- apply(dat,1,function(x) sd(x))/rowMeans(as.matrix(dat))
      vars.df <- cbind(rowMeans(as.matrix(dat)),CVs)
      colnames(vars.df) <- c("Mean","CV")
      # ordered cv by mean expression
      vars.df <- vars.df[order(vars.df[,1]),]
      # an example based on paper suggested parameter - follow python package # with window as 500
      vars <- c()
      window_size = 500
      for (k in 1:nrow(vars.df)){
        if (k <= window_size/2){
          pos <-  match(vars.df[k,2],vars.df[(1:window_size),2])   
        }else if (k >window_size/2 & k< dim(vars.df)[1]- window_size/2){
          window_order <- vars.df[c((k-window_size/2):(k+window_size/2-1)),2]
          window_order <- window_order[order(window_order)]
          
          pos <- match(vars.df[k,2], window_order)
        }
        else {pos <- match(vars.df[k,2],vars.df[c((dim(vars.df)[1]-window_size+1):(dim(vars.df)[1])),2])}
        lcv_value <- pos/(window_size)*100
        names(lcv_value) <- rownames(vars.df)[k]  
        vars <- c(vars, lcv_value)
        
      }
      vars <- vars[match(rownames(dat),names(vars))]
    }
  }
  if (!norm){
    if (metric == 'edgeR'){
      y <- edgeR::DGEList(dat,group = group)
      y <- edgeR::estimateDisp(y,robust=TRUE)
      vars <- sqrt(y$tagwise.dispersion)
      names(vars) <- rownames(dat)
    }
    else if (metric == 'DESeq2'){
      coldata <- data.frame(type = group)
      dds <- DESeq2::DESeqDataSetFromMatrix(dat, colData=coldata, design = ~type)
      dds <- DESeq2::estimateSizeFactors(dds)
      dds <- DESeq2::estimateDispersions(dds)
      vars <- log1p(DESeq2::dispersions(dds))
      names(vars) <- rownames(dat)
    }
    else if (metric == 'glmGamPoi'){
      coldata <- data.frame(type = group)
      dds <- DESeq2::DESeqDataSetFromMatrix(dat, colData=coldata, design = ~type)
      dds <- DESeq2::estimateSizeFactors(dds)
      dds <- DESeq2::estimateDispersions(dds, fitType = "glmGamPoi")
      vars <- log1p(mcols(dds)$dispGeneEst)
      names(vars) <- rownames(dat)
    }
    else if (metric == 'scran'){
      sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = dat))
      sce <- scran::computeSumFactors(sce)
      sce <- scater::logNormCounts(sce)
      dec <- scran::modelGeneVar(sce)
      vars <- dec$bio
      names(vars) <- dec@rownames
    }
    else if (metric == 'DM'){
      dat <- as.matrix(dat)
      Norm_factor <- DESeq2::estimateSizeFactorsForMatrix(dat)
      dat <- t( t(dat) / Norm_factor)
      # remove mean normalised value less than 10
      dat <- dat[rowMeans(dat)<10,]
      dat_mean <- rowMeans(dat)
      dat_sd <- matrixStats::rowSds(dat)
      dat_cv <- dat_sd/dat_mean
      dat_cv2 <- dat_cv^2
      dat_mcv2 <- cbind(log10(dat_mean),log10(dat_cv2))
      dat_mcv2 <- dat_mcv2[order(dat_mcv2[,1]),]
      # the rolling median is calculated by size 50 and the number of overlapping genes between adjacent windows is 25
      z <- zoo::zoo(dat_mcv2)
      roll_median <- zoo::rollapply(z,width=50,by=25, FUN=median) 
      # Mean corrected residual of the squared CV of gene i to its corresponding rolling median 
      DM <- c() # r(i)
      for (i in 1:nrow(dat_mcv2)){
        if(floor(i/25)<nrow(roll_median)){
          corrected.1 <- dat_mcv2[i,2]- roll_median[floor(i/25)+1,2]}else{corrected.1 <- dat_mcv2[i,2]- roll_median[nrow(roll_median),2]}
        DM <- c(DM,corrected.1)
      }
      names(DM) <- rownames(dat_mcv2)
      # To correct for the effect of gene length
      # Find gene length 
      gene_length <- goseq::getlength(names(DM),org,"geneSymbol")
      names(gene_length) <- names(DM)
      # DM only based on the genes with a known gene length 
      gene_length_na <- gene_length[!is.na(gene_length)]  
      DM.1 <-  DM[names(DM)%in%names(gene_length_na)]
      
      # Calculate second rolling median 
      DM_filtered1 <- cbind(DM.1, log10(gene_length_na))
      DM_filtered1 <- DM_filtered1[order(DM_filtered1[,1]),]
      
      z1 <- zoo::zoo(DM_filtered1)  
      roll_median1 <- zoo::rollapply(z1,width=50,by=25, FUN=median) #g(i)
      
      vars <- c() # DM = r(i)-g(i) 
      for (i in 1:length(DM.1)){
        if(floor(i/25)<dim(roll_median1))
          corrected.2 <- DM_filtered1[i,2]- roll_median1[floor(i/25)+1,2]
        else 
          corrected.2 <- DM_filtered1[i,2]- roll_median1[nrow(roll_median1),2]
        vars <- c(vars,corrected.2)
      }
      names(vars) <- names(DM.1)
    }
    else if (metric == 'Seurat_mvp'){
      obj <- Seurat::CreateSeuratObject(counts=dat)    
      obj <- Seurat::NormalizeData(obj)
      mvp <- Seurat::FindVariableFeatures(obj,selection.method = "mvp")
      vars <- Seurat::HVFInfo(mvp,selection.method = 'mvp')$dispersion.scaled
      names(vars) <- rownames(Seurat::HVFInfo(mvp,selection.method =  'mvp'))
      # remove completely no expression genes - same as others
      vars <- vars[Seurat::HVFInfo(mvp,selection.method =  'mvp')$mean > 0 ]
      vars <- vars[which(!duplicated(names(vars)))]
    }
    else if (metric == 'Seurat_vst'){
      obj <- Seurat::CreateSeuratObject(counts=dat)    
      obj <- Seurat::NormalizeData(obj)
      vst <- Seurat::FindVariableFeatures(obj,selection.method = "vst")
      vars <- Seurat::HVFInfo(vst,selection.method = 'vst')$variance.standardized
      names(vars) <- rownames(Seurat::HVFInfo(vst,selection.method =  'vst'))
      # remove completely no expression genes - same as others
      vars <- vars[Seurat::HVFInfo(vst,selection.method =  'vst')$mean > 0 ]
      vars <- vars[which(!duplicated(names(vars)))]
    }
    else if (metric == 'BASiCS'){
      sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = dat), colData = data.frame(MouseID = group))
      set.seed(1)
      Chain <- BASiCS::BASiCS_MCMC(
        Data = sce,
        N = 10000, Thin = 20, Burn = 1000,
        PrintProgress = F, Regression = TRUE,WithSpikes = FALSE 
      )
      ChainSummary <- BASiCS::Summary(Chain)
      vars <- ChainSummary@parameters[["delta"]][,1]
    }
  }
  return(vars)
}
