# scratch job script
cd4_all_rna_umap <- umap_plot(scrna_cd4_all, assay='SCT_INT', cluster_name='seurat_clusters', reduction='wnn.umap')
cd4_all_adt_umap <- umap_plot(scrna_cd4_all, assay='SCT_ADT_INT', cluster_name='seurat_clusters', reduction='wnn.umap')