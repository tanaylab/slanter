load('data.Rda')

source('slant.r')

png(filename='large_unclustered_pheatmap.png')
pheatmap::pheatmap(data, show_rownames=F, show_colnames=F, cluster_rows=F, cluster_cols=F)
dev.off()

png(filename='large_clustered_pheatmap.png')
pheatmap::pheatmap(data, show_rownames=F, show_colnames=F)
dev.off()

png(filename='large_unclustered_sheatmap.png')
sheatmap(data, show_rownames=F, show_colnames=F, cluster_rows=F, cluster_cols=F)
dev.off()

png(filename='large_reordered_sheatmap.png')
sheatmap(data, show_rownames=F, show_colnames=F)
dev.off()

png(filename='large_modified_sheatmap.png')
sheatmap(data, show_rownames=F, show_colnames=F, clusters='modify')
dev.off()

png(filename='large_replaced_sheatmap.png')
sheatmap(data, show_rownames=F, show_colnames=F, clusters='replace')
dev.off()

png(filename='large_cut_replaced_sheatmap.png')
sheatmap(data, show_rownames=F, show_colnames=F, clusters='replace', cutree_rows=4, cutree_cols=4)
dev.off()

data = data[1:7,1:7]

png(filename='small_unclustered_pheatmap.png')
pheatmap::pheatmap(data, show_rownames=T, show_colnames=T, cluster_rows=F, cluster_cols=F)
dev.off()

png(filename='small_clustered_pheatmap.png')
pheatmap::pheatmap(data, show_rownames=T, show_colnames=T)
dev.off()

png(filename='small_unclustered_sheatmap.png')
sheatmap(data, show_rownames=T, show_colnames=T, cluster_rows=F, cluster_cols=F)
dev.off()

png(filename='small_reordered_sheatmap.png')
sheatmap(data, show_rownames=T, show_colnames=T)
dev.off()

png(filename='small_modified_sheatmap.png')
sheatmap(data, show_rownames=T, show_colnames=T, clusters='modify')
dev.off()

png(filename='small_replaced_sheatmap.png')
sheatmap(data, show_rownames=T, show_colnames=T, clusters='replace')
dev.off()
