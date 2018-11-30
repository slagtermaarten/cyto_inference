RunFastMNN <- function (object.list, 
                        assay = NULL, features = 2000, reduction.name = "mnn",
    reduction.key = "mnn_", verbose = TRUE, ...) {
    # SeuratWrappers:::CheckPackage(package = "batchelor", repository = "bioconductor")
    if (!all(sapply(X = object.list, FUN = inherits, what = "Seurat"))) {
        stop("'object.list' must be a list of Seurat objects",
            call. = FALSE)
    }
    if (length(x = object.list) < 2) {
        stop("'object.list' must contain multiple Seurat objects for integration",
            call. = FALSE)
    }
    assay <- assay %||% DefaultAssay(object = object.list[[1]])
    for (i in 1:length(x = object.list)) {
        DefaultAssay(object = object.list[[i]]) <- assay
    }
    if (is.numeric(x = features)) {
        if (verbose) {
            message(paste("Computing", features, "integration features"))
        }
        features <- SelectIntegrationFeatures(object.list = object.list,
            nfeatures = features, assay = rep(assay, length(object.list)))
    }
    objects.sce <- lapply(X = object.list, FUN = function(x, f) {
        return(as.SingleCellExperiment(x = subset(x = x, features = f)))
    }, f = features)
    ## 2020-05-04 18:53 added this to make merging work
    for (i in 2:length(objects.sce)) {
      objects.sce[[i]] <- objects.sce[[i]][rownames(objects.sce[[1]]), ]
    }
    integrated <- merge(x = object.list[[1]], 
                        y = object.list[2:length(x = object.list)])
    out <- do.call(what = batchelor::fastMNN, args = c(objects.sce,
        list(...)))
    rownames(x = SingleCellExperiment::reducedDim(x = out)) <- colnames(x = integrated)
    colnames(x = SingleCellExperiment::reducedDim(x = out)) <- paste0(reduction.key,
        1:ncol(x = SingleCellExperiment::reducedDim(x = out)))
    integrated[[reduction.name]] <- CreateDimReducObject(embeddings = SingleCellExperiment::reducedDim(x = out),
        assay = DefaultAssay(object = integrated), key = reduction.key)
    Tool(object = integrated) <- out
    integrated <- LogSeuratCommand(object = integrated)
    return(integrated)
}

