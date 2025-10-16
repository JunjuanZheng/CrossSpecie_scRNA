#将下面的代码复制到文件run_monocle3.r 执行下面的命令即可进行拟时序分析
#注意输入的Seurat object的metadata中需要有一列celltype表示细胞注释类型
Rscript run_monocle3.r --rds pbmc.rds --outdir resultdir
####################################
get_earliest_principal_node <- function(cds, time_bin="2"){
    cell_ids <- which(colData(cds)[, "seurat_clusters"] == time_bin)
    closest_vertex <-cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
    closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
    root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))]
    root_pr_nodes
}

run_3d=function(cds,outdir){
    cds_3d <- reduce_dimension(cds, max_components = 3,preprocess_method = "PCA",reduction_method="UMAP")
    cds_3d <- cluster_cells(cds_3d)
    cds_3d <- learn_graph(cds_3d)
    cds_3d <- order_cells(cds_3d, root_pr_nodes=igraph::V(principal_graph(cds)[["UMAP"]])$name[1])
    cds_3d_plot_obj <- plot_cells_3d(cds_3d, color_cells_by="orig.ident")
    ggsave(paste0(outdir,"/Trajectory3D.pdf"), plot = cds_3d_plot_obj, width = 8, height = 6)
}

run_monole=function(obj,outdir,assay="RNA",batch="orig.ident",threads=8){
    if(!dir.exists(outdir)){
      dir.create(outdir,recursive=TRUE)
    }
    expr_matrix = GetAssayData(obj,assay=assay,slot="counts")
    expr_matrix = as(expr_matrix,"sparseMatrix")
    p_data <- obj@meta.data
    f_data <- data.frame(gene_short_name = rownames(obj))
    rownames(f_data)=rownames(obj)
    cds <- new_cell_data_set(expr_matrix,
                          cell_metadata = p_data,
                          gene_metadata = f_data)

    cds <- preprocess_cds(cds, num_dim = 50)
    pdf(paste0(outdir,"/pca_components.pdf"))
    plot_pc_variance_explained(cds)
    dev.off()

    cds <- align_cds(cds, alignment_group = batch)
    cds <- reduce_dimension(cds,preprocess_method = "PCA",reduction_method="UMAP")
    pdf(paste0(outdir,"/umap.pdf"))
    plot_cells(cds)
    dev.off()

    cds <- cluster_cells(cds)
    #轨迹学习Learn the trajectory graph
    cds <- learn_graph(cds)
    p = plot_cells(cds, color_cells_by = "celltype", label_groups_by_cluster=FALSE,
               label_leaves=FALSE, label_branch_points=FALSE)
    ggsave(paste0(outdir,"/Trajectory.pdf"), plot = p, width = 8, height = 6)

    #选择root
    # rownames(cds@principal_graph_aux[["UMAP"]]$dp_mst) <- NULL
    # colnames(cds@int_colData@listData$reducedDims@listData$UMAP) <- NULL

    #如果有时间信息,可以自动设置root
    #cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))
    #cds <- order_cells(cds,root_pr_nodes=igraph::V(principal_graph(cds)[["UMAP"]])$name[1])
    if(is.null){
        cds <- order_cells(cds,root_pr_nodes=igraph::V(principal_graph(cds)[["UMAP"]])$name[1])
    }else{
        cds <- order_cells(cds,root_cells=root_cell_id)
    }
    p = plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, 
               label_leaves = FALSE,  label_branch_points = FALSE)
    ggsave(paste0(outdir,"/Trajectory_Pseudotime.pdf"), plot = p, width = 8, height = 6)
    qsave(cds, file = paste0(outdir,"/cds.qs"))

    #3D轨迹
    #run_3d(cds,outdir)

    #寻找拟时轨迹差异基因
    #graph_test分析最重要的结果是莫兰指数（morans_I），其值在-1至1之间，0代表此基因没有
    #空间共表达效应，1代表此基因在空间距离相近的细胞中表达值高度相似。
    Track_genes <- graph_test(cds, neighbor_graph="principal_graph", cores=threads)
    Track_genes <- Track_genes[,c(5,2,3,4,1,6)] %>% filter(q_value < 0.05)
    write.csv(Track_genes, paste0(outdir,"/Trajectory_genes.csv"), row.names = F)

    Track_genes_sig <- Track_genes %>% top_n(n=10, morans_I) %>%
    pull(gene_short_name) %>% as.character()

    p <- plot_genes_in_pseudotime(cds[Track_genes_sig,], color_cells_by="seurat_clusters", 
                              min_expr=0.5, ncol = 2)
    ggsave(paste0(outdir,"/Genes_Jitterplot.pdf"), plot = p, width = 8, height = 6)


    p <- plot_cells(cds, genes=Track_genes_sig, show_trajectory_graph=FALSE,
                label_cell_groups=FALSE,  label_leaves=FALSE)
    p$facet$params$ncol <- 5
    ggsave(paste0(outdir,"/Genes_Featureplot.pdf"), plot = p, width = 20, height = 8)
    
    #寻找共表达基因模块
    Track_genes <- read.csv(paste0(outdir,"/Trajectory_genes.csv"))
    genelist <- pull(Track_genes, gene_short_name) %>% as.character()
    gene_module <- find_gene_modules(cds[genelist,], resolution=1e-1, cores = threads)
    write.csv(gene_module, paste0(outdir,"/Genes_Module.csv"), row.names = F)
    cell_group <- tibble::tibble(cell=row.names(colData(cds)), 
                                 cell_group=colData(cds)$seurat_clusters)
    agg_mat <- aggregate_gene_expression(cds, gene_module, cell_group)
    row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
    p <- pheatmap::pheatmap(agg_mat, scale="column", clustering_method="ward.D2")
    ggsave(paste0(outdir,"/Genes_Module.pdf"), plot = p, width = 8, height = 8)

    pseudotime <- pseudotime(cds, reduction_method = 'UMAP')
    pseudotime <- pseudotime[rownames(obj@meta.data)]
    #obj$pseudotime <- pseudotime
    #p = FeaturePlot(obj, reduction = "umap", features = "pseudotime")
    #ggsave(paste0(outdir,"/Pseudotime_Seurat.pdf"), plot = p, width = 8, height = 6)
    df = data.frame(barcodes=colnames(obj),pseudotime=pseudotime)
    write.csv(df,paste0(outdir,"/pseudotime.csv"),row.names=F)

}

library(argparse)
parser <- ArgumentParser()
parser$add_argument("--rds", type = "character", required = TRUE, help = "RDS文件路径")
parser$add_argument("--outdir", type = "character", required = TRUE, help = "输出目录路径")
parser$add_argument("--assay", type = "character", required = TRUE, help = "分析的assay名称")
parser$add_argument('--threads',help='threads',default=8,type='integer')
parser$add_argument('--batch',help='批次,default[orig.ident]',required=FALSE,default='orig.ident')
parser$add_argument("--cytotrace2",type = "character",required=FALSE,help="cytotrace2结果文件")

args <- parser$parse_args()
rds_file <- args$rds
outdir <- args$outdir
assay <- args$assay
threads = args$threads
batch = args$batch
cytotrace2 = args$cytotrace2

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(qs))
suppressPackageStartupMessages(library(monocle3))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(patchwork))

if(!is.null(cytotrace2)){
  df=read.csv(cytotrace2,row.names=1)
  df=arrange(df,desc(CytoTRACE2_Score))
  root_cell_id=rownames(df)[1]
}else{
  root_cell_id=NULL
}

obj = readRDS(rds_file)
Idents(obj)=obj$celltype
run_monole(obj,outdir,assay=assay,batch=batch,threads=threads)
