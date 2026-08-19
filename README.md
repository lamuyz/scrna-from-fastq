# scRNA-seq from FASTQ

这个仓库记录我从 **10x FASTQ 开始学习 scRNA-seq 分析** 的过程。

项目从原始测序数据出发，依次完成 reads mapping、barcode/UMI 处理、定量、质量控制及下游分析，用于梳理从 FASTQ 到单细胞表达矩阵及后续分析的完整流程。

## 学习参考

下游分析主要参考：

**Single-cell Best Practices**  
https://www.sc-best-practices.org/

实际流程会根据所使用的数据和本地计算环境进行调整。

## 数据集

使用 **10x Genomics PBMC 1k v3** 数据。

选择该数据主要考虑：

- PBMC 是常用的 scRNA-seq 示例数据；
- 数据规模较小，适合在个人电脑上运行；
- 可以用于练习从 FASTQ 到 cell type annotation 的标准分析流程。

## 软件选择

项目主要在 **Apple Silicon MacBook** 本地运行。

上游使用：

```text
piscem
+
alevin-fry
+
pyroe
```

下游主要使用：

```text
AnnData
+
Scanpy
```

整体流程：

```text
10x FASTQ
↓
splici reference
↓
piscem mapping
↓
alevin-fry quantification
↓
Gene × Cell count matrix
↓
AnnData
↓
Quality Control
↓
Normalization
↓
Dimensionality reduction
↓
Clustering
↓
Cell type annotation
```

## 现有仓库结构（更新中）

```text
scrna-from-fastq/
├── need2know.md
├── upstream/
│   └── fastq-to-anndata_workflow.md
└── downstream/
    └── QC_workflow.md
```

- [`upstream/fastq-to-anndata_workflow.md`](upstream/fastq-to-anndata_workflow.md)：FASTQ 到 AnnData 的上游流程
- [`downstream/QC_workflow.md`](downstream/QC_workflow.md)：QC 及下游分析记录
- [`need2know.md`](need2know.md)：学习过程中整理的关键概念

## 当前进度

- [x] FASTQ → AnnData
- [x] 基础 QC metrics
- [x] Ambient RNA correction
- [x] Doublet detection
- [x] QC filtering
- [x] Normalization
- [ ] Dimensionality reduction
- [ ] Clustering
- [ ] Cell type annotation

持续更新中。