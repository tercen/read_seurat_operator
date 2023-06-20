library(tercen)
library(dplyr)
library(Seurat)

ctx <- tercenCtx()

if(!any(ctx$cnames == "documentId")) stop("Column factor documentId is required") 

filename = tempfile()
writeBin(ctx$client$fileService$download(ctx$cselect("documentId")[[1]]), filename)
on.exit(unlink(filename))

doc <- ctx$client$fileService$get(ctx$cselect("documentId")[[1]])
ext <- tolower(tools::file_ext(doc$name)) 

if(ext %in% c("rdata", "rda")) {
  dat <- load(filename)
  if(length(dat) > 1) stop("A single object must be part of the RData file.")
  merged_seurat <- get(dat)  
} else if(filetype == "rds") {
  merged_seurat <- readRDS(filename)
} else {
  stop("Invalid document format.")
}

if(class(merged_seurat) != "Seurat") stop("Imported object is not a Seurat object.")

## sparse matrix to data frame
spm <- GetAssayData(merged_seurat)
tr <- try(slot(spm, "j"))
if(inherits(tr, "try-error")) {
  dff <- data.frame(
    i = spm@i + 1L,  # m@i is 0-based, not 1-based like everything else in R
    j = as.integer(rep(1:spm@Dim[2], diff(spm@p))),  # m@j is 0-based, not 1-based like everything else in R
    x = spm@x
  )
} else {
  dff <- as.data.frame(summary(spm)) %>% 
    as_tibble() %>%
    mutate(i = as.integer(i), j = as.integer(j))
}

df_out <- dff %>%
  rename(gene_id = i, cell_id = j, value = x) %>%
  mutate(.ci = 0L) %>%
  ctx$addNamespace()

gene_names <- dimnames(spm)[[1]]
df_gene <- tibble(gene_id = seq_along(gene_names), gene_names = gene_names) %>% 
  ctx$addNamespace() %>%
  rename_with(~ gsub("__remove_ns__.", "", .x, fixed = TRUE))

cell_names <- dimnames(spm)[[2]]
df_cell <- tibble(cell_id = seq_along(cell_names), cell_names = cell_names)

meta <- as_tibble(merged_seurat@meta.data, rownames = "cell_names")
df_cell <- df_cell %>%
  left_join(meta, by = "cell_names") %>% 
  ctx$addNamespace() %>%
  rename_with(~ gsub("__remove_ns__.", "", .x, fixed = TRUE))

data_relation <- df_out %>% 
  rename_with(~ gsub("__remove_ns__.", "", .x, fixed = TRUE)) %>%
  as_relation()
gene_relation <- df_gene %>% as_relation()
cell_relation <- df_cell %>% as_relation()

gid <- colnames(df_gene)[grep("gene_id", colnames(df_gene))]
cid <- colnames(df_cell)[grep("cell_id", colnames(df_cell))]

rel_out <- data_relation %>%
  left_join_relation(gene_relation, gid, gid) %>%
  left_join_relation(cell_relation, cid, cid) %>%
  as_join_operator(list(), list())

save_relation(rel_out, ctx)
