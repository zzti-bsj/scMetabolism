library(phangorn)
library(scMetabolism)
library(ggplot2)
library(rsvd)
library(AUCell)
library(GSEABase)
library(GSVA)
library(Seurat)

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
sub_sce <- young_sce
sample <- 'young'

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
##DimPlot
for(m in metabolizes){
  # 创建对应代谢图片目录
  metabolize_path <- file.path('metabolizes', sample, m)
  if(!dir.exists(metabolize_path)){
    dir.create(metabolize_path, recursive = T)
  }
  for(pathway in metabolism_types[[m]]){
    # 图片的路径
    if(pathway == "Glycolysis / Gluconeogenesis"){
      p_path <- file.path(metabolize_path, sprintf("%s%s", "Glycolysis - Gluconeogenesis", '.png'))
    }else{
      p_path <- file.path(metabolize_path, sprintf("%s%s", pathway, '.png'))
    }
    # 做代谢的DimPlot
    print(p_path)
    p <- DimPlot.metabolism(sub_sce, pathway = pathway, dimention.reduction.type = "umap", dimention.reduction.run = F, size = 1)
    ggsave(p_path, p)
  }
}

##DotPlot&BoxPlot
###young
for(m in metabolizes){
    # 创建对应代谢图片目录
    pathways <- unlist(metabolism_types[m])
    # # DotPlot
    p_path <- file.path('metabolizes', 'young', sprintf("%s%s%s", m, '-DotPlot', '.png'))
    p <- DotPlot.metabolism(obj = young_metabolism, pathway = pathways, phenotype = "celltypes", norm = "y", title = m)
    ggsave(p_path, p)
    # BoxPlot
    p_path <- file.path('metabolizes', 'young', sprintf("%s%s%s", m, '-BoxPlot', '.png'))
    p <- BoxPlot.metabolism(obj = young_metabolism, pathway = pathways, phenotype = "celltypes", ncol = 4, title = m)
    ggsave(p_path, p, height = length(pathways) / 4 * 4, width = 16)
}
###old
for(m in metabolizes){
    # 创建对应代谢图片目录
    pathways <- unlist(metabolism_types[m])
    # DotPlot
    p_path <- file.path('metabolizes', 'old', sprintf("%s%s%s", m, '-DotPlot', '.png'))
    p <- DotPlot.metabolism(obj = old_metabolism, pathway = pathways, phenotype = "celltypes", norm = "y", title = m)
    ggsave(p_path, p)
    # BoxPlot
    p_path <- file.path('metabolizes', 'old', sprintf("%s%s%s", m, '-BoxPlot', '.png'))
    p <- BoxPlot.metabolism(obj = old_metabolism, pathway = pathways, phenotype = "celltypes", ncol = 4, title = m)
    ggsave(p_path, p, height = length(pathways) / 4 * 4, width = 16)
}
### young/old
for(m in metabolizes){
    # 创建对应代谢图片目录
    pathways <- unlist(metabolism_types[m])
    # BoxPlot
    p_path <- file.path('metabolizes', sprintf("%s%s%s", m, '-BoxPlot', '.png'))
    p1 <- BoxPlot.metabolism(obj = old_metabolism, pathway = pathways, phenotype = "celltypes", ncol = 4, title = sprintf('%s - %s', 'Young', m))
    p2 <- BoxPlot.metabolism(obj = young_metabolism, pathway = pathways, phenotype = "celltypes", ncol = 4, title = sprintf('%s - %s', 'Old', m))
    
    # 判断一下，组织图片的方式
    # 处理宽高占比
    if(length(pathways) > 8){
      p <- p1 + p2
      height <- length(pathways) / 4 * 4
      width <- 32
    }else{
      p <- p1 / p2
      height <- length(pathways) / 4 * 4 * 2
      width <- 16
    }
    ggsave(p_path, p, height = height, width = width)
}

# diffDotPlot & Combined DotPlot
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
