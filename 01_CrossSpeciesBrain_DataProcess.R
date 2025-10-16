library(Seurat, help, pos = 2, lib.loc = NULL)
RenameGenesSeurat <- function(obj = ls.Seurat[[1]], newnames = HGNC.updated[[1]]$Suggested.Symbol) { # Replace gene
  print("Run this before integration. It only changes obj@assays$RNA@counts, @data and @scale.data.")
  RNA <- obj@assays$RNA
  
  if (nrow(RNA) == length(newnames)) {
    if (length(RNA@counts)) RNA@counts@Dimnames[[1]]      <- newnames
    if (length(RNA@data)) RNA@data@Dimnames[[1]]          <- newnames
    if (length(RNA@scale.data)) RNA@scale.data@Dimnames[[1]] <- newnames
  } else {"Unequal gene sets: nrow(RNA) != nrow(newnames)"}
  obj@assays$RNA <- RNA
  return(obj)
}
projectPath = '/mnt2/wanggd_group/zjj/BGCscRNA/ZhaoKai'
setwd(projectPath)
dataPath = paste0(projectPath,'/Data/GSE127774/')

data1 = readRDS(paste0(dataPath,'GSE127774_CER_seurat.rds'))
data1=UpdateSeuratObject(data1)

data2 = readRDS(paste0(dataPath,'GSE127774_CN_seurat.rds'))
data2=UpdateSeuratObject(data2)

data3 = readRDS(paste0(dataPath,'GSE127774_ACC_seurat.rds'))
data3=UpdateSeuratObject(data3)

dataTemp = data3
dataTemp <- RunUMAP(dataTemp, reduction = "pca", dims = 1:30)
### 转换基因名
library(biomaRt)
# 查询ensembl
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# 获取转换表
DefaultAssay(data3)  # 应该显示 "integrated"
length(rownames(dataTemp[["RNA"]]))        # 应该是 13010
length(rownames(dataTemp[["integrated"]])) # 应该是 2000
genes <- rownames(dataTemp[["integrated"]])
gene_map <- getBM(
    attributes = c('ensembl_gene_id', 'hgnc_symbol'),
    filters    = 'ensembl_gene_id',
    values     = genes,
    mart       = mart
)
# 查看前几行
head(gene_map)
dim(gene_map)
# 删除没有symbol的行，只保留有symbol的
gene_map <- gene_map[gene_map$hgnc_symbol != "", ]
# 为了防止symbol重复，可以这样处理（取第一个出现的）
gene_map <- gene_map[!duplicated(gene_map$ensembl_gene_id), ]
gene_map <- gene_map[!duplicated(gene_map$hgnc_symbol), ]
#查看目的基因是否存在
gene_map[gene_map$hgnc_symbol %in% c("IGF1", "IGF1R"), ]

# 匹配新symbol并赋值给Seurat对象
subDataTemp = subset(dataTemp,features = gene_map[,1])
