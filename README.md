# scRNA-seq from FASTQ

这个仓库记录我从 **10x FASTQ 开始学习并实际完成 scRNA-seq 分析流程** 的过程。

## 学习参考

下游分析主要参考：

- [Single-cell Best Practices](https://www.sc-best-practices.org/)

实际流程根据数据特点和本地计算环境进行调整。

## 数据

使用 **10x Genomics PBMC 1k v3** 数据。

选择该数据主要是因为：

- PBMC 是常用的 scRNA-seq 示例数据；
- 数据规模较小，适合在个人电脑上运行；
- 可以完整练习从 FASTQ 到 cell type annotation 的分析流程。

## 分析环境与工具

项目主要在 **Apple Silicon MacBook** 本地运行。

上游主要使用：

```text
piscem
+
alevin-fry
+
pyroe
```

下游主要使用：

```text
SoupX / scDblFinder
↓
AnnData / Scanpy
```

其中：

- `piscem`：建立索引并进行 reads mapping；
- `alevin-fry`：处理 barcode / UMI 并生成表达矩阵；
- `pyroe`：读取 alevin-fry 输出并构建 AnnData；
- `SoupX`：ambient RNA correction；
- `scDblFinder`：doublet detection；
- `Scanpy`：QC、normalization、feature selection、降维、聚类及后续分析。

## 整体流程

```text
10x FASTQ
↓
Genome FASTA + GENCODE GTF
↓
splici reference
↓
piscem index
↓
piscem map-sc
↓
RAD
↓
alevin-fry
↓
Gene × Cell count matrix
↓
pyroe.load_fry()
↓
AnnData
↓
QC
↓
Ambient RNA correction
↓
Doublet detection
↓
QC filtering
↓
Normalization
↓
Feature selection
↓
PCA
↓
Neighbor graph
↓
UMAP
↓
Leiden clustering
↓
Cell type annotation
```

## 仓库结构（更新中）

```text
scrna-from-fastq/
├── README.md
├── need2know.md
│
├── upstream/
│   └── fastq-to-anndata_workflow.md
│
└── downstream/
    ├── analysis/
    │   ├── QC_github.ipynb
    │   ├── QC_workflow.md
    │   ├── normalization_github.ipynb
    │   ├── feature_selection_github.ipynb
    │   ├── dimensionality_reduction_github.ipynb
    │   └── clustering_github.ipynb
    │
    └── notes/
        ├── qc.md
        ├── normalization.md
        ├── feature_selection.md
        ├── dimensionality_reduction.md
        ├── clustering.md
        └── annotation.md
```

### Upstream

[`upstream/fastq-to-anndata_workflow.md`](upstream/fastq-to-anndata_workflow.md)

记录从 FASTQ、参考序列构建、piscem mapping、alevin-fry quantification 到 AnnData 的流程。

### Downstream analysis

[`downstream/analysis/`](downstream/analysis/)

保存使用真实 PBMC 数据进行各分析步骤的 notebook 和 workflow。

### Notes

[`downstream/notes/`](downstream/notes/)

记录各个分析步骤的学习笔记，包括：

```text
QC
↓
Normalization
↓
Feature selection
↓
Dimensionality reduction
↓
Clustering
↓
Cell type annotation
```

### Concepts

[`need2know.md`](need2know.md)

记录分析过程中涉及的基础概念、术语和工具原理。

## 当前进度

- [x] FASTQ → AnnData
- [x] QC
- [x] Ambient RNA correction
- [x] Doublet detection
- [x] QC filtering
- [x] Normalization
- [x] Feature selection
- [x] Dimensionality reduction
- [x] Clustering
- [ ] Cell type annotation

持续更新中。