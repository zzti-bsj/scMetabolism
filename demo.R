library(phangorn)
library(scMetabolism)
library(ggplot2)
library(rsvd)
library(AUCell)
library(GSEABase)
library(GSVA)
library(Seurat)


#########definition
# 使用diffData 和 diffDotPlot函数来构建代谢上调下调的点图
diffData <- function(obj, pathway, phenotype, norm = "y"){
  input.norm = norm
  input.pathway <- pathway
  input.parameter<-phenotype

  metadata<-obj@meta.data
  metabolism.matrix <- obj@assays$METABOLISM$score

  cat("\nPlease Cite: \nYingcheng Wu, Qiang Gao, et al. Cancer Discovery. 2021. \nhttps://pubmed.ncbi.nlm.nih.gov/34417225/   \n\n")


  metadata[,input.parameter]<-as.character(metadata[,input.parameter])
  metabolism.matrix_sub<-t(metabolism.matrix[input.pathway,])

  #arrange large table
  gg_table<-c()
  for (i in 1:length(input.pathway)){
    gg_table<-rbind(gg_table, cbind(metadata[,input.parameter], input.pathway[i], metabolism.matrix_sub[,i]))
  }
  gg_table<-data.frame(gg_table)

  #get median value
  gg_table_median<-c()
  input.group.x<-unique(as.character(gg_table[,1]))
  input.group.y<-unique(as.character(gg_table[,2]))


  for (x in 1:length(input.group.x)){
    for (y in 1:length(input.group.y)){
      gg_table_sub<-subset(gg_table, gg_table[,1] == input.group.x[x] & gg_table[,2] == input.group.y[y])
      gg_table_median<-rbind(gg_table_median, cbind(input.group.x[x], input.group.y[y], median(as.numeric(as.character(gg_table_sub[,3])))))

    }
  }
  gg_table_median<-data.frame(gg_table_median)
  gg_table_median[,3]<-as.numeric(as.character(gg_table_median[,3]))


  #normalize
  gg_table_median_norm<-c()
  input.group.x<-unique(as.character(gg_table[,1]))
  input.group.y<-unique(as.character(gg_table[,2]))


  range01 <- function(x){(x-min(x))/(max(x)-min(x))}

  if (input.norm == "y")
    for (y in 1:length(input.group.y)){
      gg_table_median_sub<-subset(gg_table_median, gg_table_median[,2] == input.group.y[y])
      norm_value<- range01(as.numeric(as.character(gg_table_median_sub[,3])))
      gg_table_median_sub[,3]<-norm_value
      gg_table_median_norm<-rbind(gg_table_median_norm, gg_table_median_sub)
    }

  if (input.norm == "x")
    for (x in 1:length(input.group.x)){
      gg_table_median_sub<-subset(gg_table_median, gg_table_median[,1] == input.group.x[x])
      norm_value<- range01(as.numeric(as.character(gg_table_median_sub[,3])))
      gg_table_median_sub[,3]<-norm_value
      gg_table_median_norm<-rbind(gg_table_median_norm, gg_table_median_sub)
    }

  if (input.norm == "na") gg_table_median_norm<-gg_table_median


  gg_table_median_norm<-data.frame(gg_table_median_norm)
  gg_table_median_norm[,3]<-as.numeric(as.character(gg_table_median_norm[,3]))

  gg_table_median_norm
}

diffDotPlot <- function(gg_table_median_norm, title) {
  library(wesanderson)
  pal <- wes_palette("Zissou1", 100, type = "continuous")
  ggplot(data=gg_table_median_norm, aes(x=gg_table_median_norm[,1], y=gg_table_median_norm[,2], color = gg_table_median_norm[,3])) +
    geom_point(data=gg_table_median_norm, aes(size = gg_table_median_norm[,3])) + #geom_line() +
    #theme_bw()+theme(aspect.ratio=0.5, axis.text.x = element_text(angle = 45, hjust = 1)) +
    ylab("Metabolic Pathway")+ xlab('celltypes') + labs(title = title) +
    theme_bw()+theme(axis.text.x = element_text(angle = 45, hjust = 1), #aspect.ratio=1,
                     panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
    scale_color_gradientn(colors = c("blue", "white", "red"), limits = c(-1, 1)) +
    labs(color = "Value", size = "Value") +
    #facet_wrap(~tissueunique, ncol = 1) +
    #theme_bw()+theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    NULL
}


# 联合箱线图 jointBoxPlot
jointBoxPlot <- function(obj_list, pathway, phenotype, ncol = 1, title = 'Metabolic DotPlot'){
  # gg_table_list
  gg_table_df <- data.frame()
  for (obj in obj_list){
    input.pathway<-pathway
    input.parameter<-phenotype

    cat("\nPlease Cite: \nYingcheng Wu, Qiang Gao, et al. Cancer Discovery. 2021. \nhttps://pubmed.ncbi.nlm.nih.gov/34417225/   \n\n")

    metadata <-obj@meta.data
    metabolism.matrix <- obj@assays$METABOLISM$score

    metadata[,input.parameter]<-as.character(metadata[,input.parameter])
    metabolism.matrix_sub<-t(metabolism.matrix[input.pathway,])

    #arrange large table
    gg_table<-c()
    for (i in 1:length(input.pathway)){
      gg_table<-rbind(gg_table, cbind(metadata[,input.parameter], input.pathway[i], metabolism.matrix_sub[,i]))
    }
    # 修改gg_table, 增加sample_origin列
    gg_table<-data.frame(gg_table)
    gg_table$sample_origin <- unique(obj$sample_origin)
    gg_table[,3]<-as.numeric(as.character(gg_table[,3]))
    
    gg_table_df <- rbind(gg_table_df, gg_table)
  }

  # 合并两个dataframe
  print(gg_table_df[1:4,])
  print(unique(gg_table_df[,4]))

  # 长度
  print(length(gg_table_df))
  print(nrow(gg_table_df))
  print(ncol(gg_table_df))

  library(wesanderson)
  pal <- wes_palette("Zissou1", 100, type = "continuous")
  ggplot(data=gg_table_df, aes(x=gg_table_df[,1], y=gg_table_df[,3], fill = gg_table_df[,4])) +
    geom_boxplot(outlier.shape=NA)+
    ylab("Metabolic Pathway")+
    xlab(input.parameter)+
    labs(title = title, hjust = 0.5) +
    theme_bw()+theme(axis.text.x = element_text(angle = 45, hjust = 1), #aspect.ratio=1,
                     panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
                     plot.title = element_text(hjust = 0.5)) +
    # scale_color_gradientn(colours = pal) +
    facet_wrap(~gg_table_df[,2], ncol = ncol, scales = "free") +
    labs(fill = 'Group') +
    #theme_bw()+theme(axis.text.x = element_text(angle = 45, hjust = 1)) +

    NULL
}


############ callable
# read sce data
sce <- readRDS(sprintf("%s/sce_Annotation.RDS", codeDataDirR))

# 根据年龄设置sample_origin
# unique(sce$orig.ident)
# old <- c("HM1", "HM2", "HM3", "HM4", "HM6", "HM7", "HM9", "HM11", "HM12", "HM15", "HM23")
# young <- c("HM5", "HM8", "HM10", "HM13", "HM21", "HM22")
# sce$sample_origin <- ifelse(sce[["orig.ident"]][[1]] %in% old, 'old', 'young')
# saveRDS(sce, sprintf("%s/sce_Annotation.RDS", codeDataDirR))

young_sce <- subset(sce, sample_origin == 'young')
old_sce <- subset(sce, sample_origin == 'old')

# ======== 设置sub_set
# sub_sce <- young_sce
# sample <- 'young'

# 1. 使用sce的counts来计算代谢通路的评分
young_metabolism <- sc.metabolism.Seurat(young_sce, method = "AUCell", imputation = F, ncores = 2, metabolism.type = "KEGG")
old_metabolism <- sc.metabolism.Seurat(old_sce, method = "AUCell", imputation = F, ncores = 2, metabolism.type = "KEGG")

#要研究的代谢通路列表
metabolism_types <- list(
  "Carbohydrate Metabolism" = c( # 糖代谢
    "Glycolysis / Gluconeogenesis", # 糖酵解/糖异生
    "Pentose Phosphate Pathway",    # 磷酸戊糖途径
    "Glycogen Synthesis",           # 糖原合成
    "Glycogen Breakdown"            # 糖原分解
  ),
  "Mitochondrial Energy Metabolism" = c( # 线粒体代谢
    "Citrate cycle (TCA cycle)",    # 柠檬酸循环（三羧酸循环）
    "Pyruvate metabolism",          # 丙酮酸代谢
    "Oxidative phosphorylation",    # 氧化磷酸化
    "Beta-oxidation"                # β-氧化
  ),
  "Lipid Metabolism" = c( # 脂代谢
    "Fatty acid biosynthesis",       # 脂肪酸生物合成
    "Fatty acid elongation",         # 脂肪酸延长
    "Fatty acid degradation",        # 脂肪酸降解
    "Synthesis and degradation of ketone bodies",  # 酮体的合成和降解
    "Steroid biosynthesis",          # 固醇生物合成
    "Primary bile acid biosynthesis", # 初级胆汁酸生物合成
    "Steroid hormone biosynthesis",  # 类固醇激素生物合成
    "Glycerolipid metabolism",       # 甘油酯代谢
    "Glycerophospholipid metabolism", # 甘油磷脂代谢
    "Ether lipid metabolism",        # 醚脂代谢
    "Sphingolipid metabolism",       # 鞘脂代谢
    "Arachidonic acid metabolism",   # 花生四烯酸代谢
    "Linoleic acid metabolism",      # 亚油酸代谢
    "alpha-Linolenic acid metabolism", # α-亚麻酸代谢
    "Biosynthesis of unsaturated fatty acids"  # 不饱和脂肪酸的生物合成
  ),
  "Amino Acid Metabolism" = c( # 氨基酸代谢
    "Alanine, aspartate and glutamate metabolism", # 丙氨酸、天冬氨酸和谷氨酸代谢
    "Glycine, serine and threonine metabolism",    # 甘氨酸、丝氨酸和苏氨酸代谢
    "Cysteine and methionine metabolism",          # 半胱氨酸和蛋氨酸代谢
    "Valine, leucine and isoleucine degradation",  # 缬氨酸、亮氨酸和异亮氨酸降解
    "Valine, leucine and isoleucine biosynthesis", # 缬氨酸、亮氨酸和异亮氨酸生物合成
    "Lysine degradation",                          # 赖氨酸降解
    "Arginine biosynthesis",                       # 精氨酸生物合成
    "Arginine and proline metabolism",             # 精氨酸和脯氨酸代谢
    "Histidine metabolism",                        # 组氨酸代谢
    "Tyrosine metabolism",                         # 酪氨酸代谢
    "Phenylalanine metabolism",                    # 苯丙氨酸代谢
    "Tryptophan metabolism",                       # 色氨酸代谢
    "Phenylalanine, tyrosine and tryptophan biosynthesis", # 苯丙氨酸、酪氨酸和色氨酸生物合成
    "beta-Alanine metabolism",                     # β-丙氨酸代谢
    "Taurine and hypotaurine metabolism",          # 牛磺酸和高牛磺酸代谢
    "Glutathione metabolism"                       # 谷胱甘肽代谢
  )
)
metabolizes <- names(metabolism_types)

# 2. visualize
##2.1 DimPlot 绘制umap图
# for(m in metabolizes){
#   # 创建对应代谢图片目录
#   metabolize_path <- file.path('metabolizes', sample, m)
#   if(!dir.exists(metabolize_path)){
#     dir.create(metabolize_path, recursive = T)
#   }
#   for(pathway in metabolism_types[[m]]){
#     # 图片的路径
#     if(pathway == "Glycolysis / Gluconeogenesis"){
#       p_path <- file.path(metabolize_path, sprintf("%s%s", "Glycolysis - Gluconeogenesis", '.png'))
#     }else{
#       p_path <- file.path(metabolize_path, sprintf("%s%s", pathway, '.png'))
#     }
#     # 做代谢的DimPlot
#     print(p_path)
#     p <- DimPlot.metabolism(young_metabolism, pathway = pathway, dimention.reduction.type = "umap", dimention.reduction.run = F, size = 1)
#     ggsave(p_path, p)
#   }
# }

##2.3 DotPlot&BoxPlot（基础BoxPlot和DotPlot的Demo，可作为参考）
###young
# for(m in metabolizes){
#     # 创建对应代谢图片目录
#     pathways <- unlist(metabolism_types[m])
#     # # DotPlot
#     p_path <- file.path('metabolizes', 'young', sprintf("%s%s%s", m, '-DotPlot', '.png'))
#     p <- DotPlot.metabolism(obj = young_metabolism, pathway = pathways, phenotype = "celltypes", norm = "y", title = m)
#     ggsave(p_path, p)
#     # BoxPlot
#     p_path <- file.path('metabolizes', 'young', sprintf("%s%s%s", m, '-BoxPlot', '.png'))
#     p <- BoxPlot.metabolism(obj = young_metabolism, pathway = pathways, phenotype = "celltypes", ncol = 4, title = m)
#     ggsave(p_path, p, height = length(pathways) / 4 * 4, width = 16)
# }

# ###old
# for(m in metabolizes){
#     # 创建对应代谢图片目录
#     pathways <- unlist(metabolism_types[m])
#     # DotPlot
#     p_path <- file.path('metabolizes', 'old', sprintf("%s%s%s", m, '-DotPlot', '.png'))
#     p <- DotPlot.metabolism(obj = old_metabolism, pathway = pathways, phenotype = "celltypes", norm = "y", title = m)
#     ggsave(p_path, p)
#     # BoxPlot
#     p_path <- file.path('metabolizes', 'old', sprintf("%s%s%s", m, '-BoxPlot', '.png'))
#     p <- BoxPlot.metabolism(obj = old_metabolism, pathway = pathways, phenotype = "celltypes", ncol = 4, title = m)
#     ggsave(p_path, p, height = length(pathways) / 4 * 4, width = 16)
# }

## 2.4 StackedBoxPlot（迭代产物，不用于生产，可Delete）
### young/old
# for(m in metabolizes){
#     # 创建对应代谢图片目录
#     pathways <- unlist(metabolism_types[m])
#     # BoxPlot
#     p_path <- file.path('metabolizes', sprintf("%s%s%s", m, '-BoxPlot', '.png'))
#     p1 <- BoxPlot.metabolism(obj = old_metabolism, pathway = pathways, phenotype = "celltypes", ncol = 4, title = sprintf('%s - %s', 'Young', m))
#     p2 <- BoxPlot.metabolism(obj = young_metabolism, pathway = pathways, phenotype = "celltypes", ncol = 4, title = sprintf('%s - %s', 'Old', m))
    
#     # 判断一下，组织图片的方式
#     # 处理宽高占比
#     if(length(pathways) > 8){
#       p <- p1 + p2
#       height <- length(pathways) / 4 * 4
#       width <- 32
#     }else{
#       p <- p1 / p2
#       height <- length(pathways) / 4 * 4 * 2
#       width <- 16
#     }
#     ggsave(p_path, p, height = height, width = width)
# }

##2.5 diffDotPlot & Combined DotPlot （用于生产环境的DotPlot）
for(m in metabolizes){
    # 创建对应代谢图片目录
    pathways <- unlist(metabolism_types[m])
    # DotPlot
    diffData_young <- diffData(obj = young_metabolism, pathway = pathways, phenotype = "celltypes", norm = "y")
    diffData_old <- diffData(obj = old_metabolism, pathway = pathways, phenotype = "celltypes", norm = "y")

    # 假设 diffData_old 和 diffData_young 已经被正确加载为数据框
    # 首先将 NA 替换为 0
    diffData_old$X3[is.na(diffData_old$X3)] <- 0
    diffData_young$X3[is.na(diffData_young$X3)] <- 0

    # 确保行顺序一致，这里假设它们已经是一致的
    # 进行相减操作
    # young - old
    diffData_diff <- data.frame(
      X1 = diffData_young$X1,
      X2 = diffData_young$X2,
      X3 = diffData_young$X3 - diffData_old$X3
    )

    # young,old sample dotplot
    p_young <- DotPlot.metabolism(obj = young_metabolism, pathway = pathways, phenotype = "celltypes", norm = "y", title = m)
    p_old <- DotPlot.metabolism(obj = old_metabolism, pathway = pathways, phenotype = "celltypes", norm = "y", title = m)
    # diff dotplot
    p_diff <- diffDotPlot(diffData_diff, sprintf("%s%s", 'Young vs Old / ', m))
    # save diff dotplot
    # p_path <- file.path('metabolizes', sprintf("%s%s%s", m, '-Diff-DotPlot', '.png'))
    # ggsave(p_path, p_diff)
    # save combined dotplot
    p <- p_young + p_old + p_diff
    combined_path <- file.path('metabolizes', sprintf("%s%s%s", m, 'Combined-DotPlot', '.png'))
    ggsave(combined_path, p, width = 24, height = 8)

}


obj_list <- list(young_metabolism, old_metabolism)
# 2.6 Combined Box Plot（联合箱线图）

##2.6.1 SplitedBoxPlot（单张的BoxPlot）
# for(m in metabolizes){
#   pathways <- unlist(metabolism_types[m])
#   for(pathway in pathways){
#     p <- jointBoxPlot(obj_list = obj_list, pathway = pathway, phenotype = "celltypes", title = 'BoxPlot')
#     # 图片的路径
#     if(pathway == "Glycolysis / Gluconeogenesis"){
#       p_path <- file.path('metabolizes', m, sprintf("%s%s", "Glycolysis - Gluconeogenesis", '.png'))
#     }else{
#       p_path <- file.path('metabolizes', m, sprintf("%s%s", pathway, '.png'))
#     }
#     # 保存图片
#     ggsave(p_path, p, height = 8, width = 8)
#   }
# }

##2.6.2 StackedBoxPlot（堆叠的BoxPlot）
for(m in metabolizes){
  pathways <- unlist(metabolism_types[m])
  p <- jointBoxPlot(obj_list = obj_list, pathway = pathways, phenotype = "celltypes", title = 'BoxPlot', ncol = 4)
  # 图片的路径
  if(pathway == "Glycolysis / Gluconeogenesis"){
    p_path <- file.path('metabolizes', sprintf("%s%s", m, '-CombinedBox.png'))
  } else {
    p_path <- file.path('metabolizes', sprintf("%s%s", m, '-CombinedBox.png'))
  }
  # 处理宽高占比
  # 处理图片
  if(length(pathways) > 8){
    height <- length(pathways) / 4 * 8
    width <- 32
  }else{
    height <- length(pathways) / 4 * 4 * 3
    width <- 24
  }
  ggsave(p_path, p, height = height, width = width)
}
