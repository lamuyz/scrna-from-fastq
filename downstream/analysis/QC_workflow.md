# QC_Workflow

该文件记录 10x PBMC 1k 数据的实际 QC 流程。方法主要参考 [Single-cell Best Practices](https://www.sc-best-practices.org/preprocessing_visualization/quality_control.html)，具体阈值根据当前数据分布调整。

## 1. Input

初始 AnnData：

```text
1231 cells × 78941 genes
```

主要使用 `.X` 中的 `spliced + ambiguous counts`。

## 2. QC metrics

根据 gene symbol 标记 mitochondrial、ribosomal 和 hemoglobin genes，并计算：

```text
total_counts
n_genes_by_counts
pct_counts_mt
pct_counts_ribo
pct_counts_hb
pct_counts_in_top_20_genes
```

随后检查 violin plot、scatter plot 和 histogram。

## 3. MAD outlier inspection

按照教程示例，对以下指标使用 5 MAD：

```text
log1p_total_counts
log1p_n_genes_by_counts
pct_counts_in_top_20_genes
```

mitochondrial outlier 使用：

```text
3 MAD
OR
pct_counts_mt > 8%
```

但当前数据的 median `pct_counts_mt` 约为 11%，如果直接使用 8% cutoff，会标记绝大多数 cells。

因此这里没有直接按 8% 删除，而是把它作为异常参考，后续结合其他 QC 指标、聚类和 marker 继续判断。

## 4. SoupX

使用 unfiltered matrix 和 selected cell matrix 进行 ambient RNA correction。

由于 4 个 selected barcodes 不在 raw-like matrix 中，SoupX 实际使用：

```text
1227 cells × 78941 genes
```

在 selected cells 的副本上进行简单聚类，为 SoupX 提供 cluster labels。该 clustering 只用于 contamination estimation，不用于直接过滤 cells。

SoupX 校正后，将 corrected counts 写回 AnnData。

## 5. Doublet detection

在 SoupX-corrected counts 上运行 `scDblFinder`：

```text
1189 singlets
38 doublets
```

结合 `total_counts`、`n_genes_by_counts` 和 cluster distribution 检查后，去除 38 个 doublets。

## 6. Independent QC review

对 1189 个 singlets 重新进行 PCA、UMAP 和 Leiden clustering，并结合 marker 和 QC metrics 检查各 cluster。

其中一个 cluster 主要由 mitochondrial genes 构成，并表现出明显的 low-quality 特征，因此去除：

```text
37 cells
```

剩余：

```text
1152 cells
```

## 7. Final cell-level filtering

最后只过滤明显异常的 cells：

```text
pct_counts_mt >= 50%

OR

pct_counts_mt >= 25%
AND
n_genes_by_counts < 500

OR

n_genes_by_counts < 200
```

共去除：

```text
12 cells
```

得到：

```text
1140 cells × 78941 genes
```

## 8. Gene filtering

最后过滤只在少量 cells 中出现的 genes：

```python
sc.pp.filter_genes(adata, min_cells=20)
```

最终得到：

```text
1140 cells × 12451 genes
```

## Workflow summary

```text
1231 cells
↓
QC metric inspection
↓
SoupX
↓
1227 cells
↓
scDblFinder
↓
-38 doublets
↓
1189 singlets
↓
independent clustering + QC review
↓
-37 MT-high / low-quality cells
↓
1152 cells
↓
final cell-level QC
↓
-12 cells
↓
1140 cells
↓
filter_genes(min_cells=20)
↓
1140 cells × 12451 genes
```

最终 QC 后的 count matrix 用于后续 normalization。