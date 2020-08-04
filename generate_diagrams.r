load('data.Rda')

png(filename='clustered_pheatmap.png')
pheatmap::pheatmap(data, show_rownames=F, show_colnames=F)
dev.off()

source('slant.r')

png(filename='unclustered_sheatmap.png')
sheatmap(data, show_rownames=F, show_colnames=F)
dev.off()

png(filename='clustered_sheatmap.png')
sheatmap(data, show_rownames=F, show_colnames=F, cluster_rows=T, cluster_cols=T)
dev.off()
