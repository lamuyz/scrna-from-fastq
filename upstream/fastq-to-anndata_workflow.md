
# scRNA-seq 上游流程：FASTQ → AnnData

整体流程是：

```text
FASTQ
↓
reference 构建
↓
piscem index
↓
piscem map-sc
↓
map.rad
↓
alevin-fry generate-permit-list
↓
alevin-fry collate
↓
alevin-fry quant
↓
USA-mode count matrix
↓
pyroe.load_fry()
↓
AnnData
↓
gene metadata annotation
```

---

## 1. 原始数据：10x Chromium v3 FASTQ

使用的是 10x Genomics PBMC 1k v3 数据。

这套 **Chromium Single Cell 3′ Gene Expression v3 chemistry** 中：

```text
R1 = 28 bp
├── 16 bp Cell barcode
└── 12 bp UMI

R2 = 91 bp
└── cDNA sequence，用于 mapping
```

另外还有：

```text
I1 = Index Read 1
→ sample index，用于区分样本
```

其中：

- **Cell barcode（细胞条形码）**：告诉我们这条文库分子来自哪个 cell。
    
- **UMI = Unique Molecular Identifier，唯一分子标识符**：用于区分原始 RNA molecules，减少 PCR duplicate 被重复计数。
    
- **R2**：携带 cDNA sequence，用于判断它可能来源于哪个 transcript/reference，最终再归到 gene。
    

严格来说，R1 和 R2 各自都是一条 read，它们组成一个 **read pair（成对读段）**，来源于同一个 library molecule。

---

## 2. 准备 genome 和 gene annotation

下载了：

```text
GRCh38.primary_assembly.genome.fa.gz
gencode.v50.primary_assembly.annotation.gtf.gz
```

其中：

**FASTA**（参考基因组的实际 DNA sequence）：

```text
GRCh38 genome sequence
```

**GTF = Gene Transfer Format，基因转移格式**：

```text
哪个 chromosome
哪个位置
是 gene / transcript / exon
gene_id 是什么
gene_name 是什么
```

所以可以理解为：

```text
FASTA
→ 地图本身

GTF
→ 地图上的注释
```

---

## 3. 构建 splici reference

使用：

```bash
pyroe make-spliced+intronic
```

根据 genome FASTA + GTF 构建了 **splici reference**。

**splici = spliced + intronic**

其中：

- **spliced = 已剪接**
    
- **intronic = 内含子相关的**
    

这样 reference 不仅包含已剪接 transcript 相关序列，也保留 intronic sequence，使后续能够利用未剪接 RNA 的信息。

生成的重要文件包括：

```text
gencode_v50_splici_fl86.fa
gencode_v50_splici_fl86_t2g_3col.tsv
gene_id_to_name.tsv
duplicate_entries.tsv
```

其中：

```text
t2g = transcript/reference → gene mapping
```

告诉后面的 quantification：

> 某个被 mapping 到的 reference/transcript 最终属于哪个 gene。

---

## 4. 构建 piscem index

使用 splici FASTA：

```text
gencode_v50_splici_fl86.fa
```

运行：

```bash
piscem build
```

生成 piscem index。

index 的目的不是产生表达矩阵，而是提前把 reference 组织成适合快速查询的数据结构。

后面 piscem mapping 时，就不需要每条 R2 都从头扫描整个 3.7 GB reference。

这里涉及：

- **k-mer**：长度为 k 的连续序列片段。
    
- `k=31`
    
- `m=19`：用于索引和快速候选搜索的较短序列长度。
    

---

## 5. piscem mapping：FASTQ → RAD

运行：

```bash
piscem map-sc
```

使用：

```text
R1
→ barcode + UMI

R2
→ biological read
```

以及：

```text
-g chromium_v3
```

告诉 piscem 这批 reads 的 **geometry（读段结构）**：

```text
R1:
16 bp barcode
12 bp UMI

R2:
全部用于 mapping
```

mapping 结果：

```text
Mapped 58,983,500 / 66,601,887 reads
≈ 88.6%
```

这里的 `58,983,500` 对应成功进入 RAD mapping records 的 read-pair/library-fragment observations，不应该理解成 58,983,500 个 RNA molecules。

输出最重要的文件：

```text
map.rad
```

**RAD = Reduced Alignment Data，精简比对数据。**

它不像 BAM 那样保存完整 alignment 细节，而是保留后续单细胞定量真正需要的 mapping 信息，例如：

```text
cell barcode
UMI
mapping targets
```

因此：

```text
FASTQ
↓
piscem
↓
RAD
```

---

## 6. generate-permit-list：处理 cell barcode

接下来使用：

```bash
alevin-fry generate-permit-list
```

这一步主要处理：

**哪些 barcode 可以进入后续分析，以及部分错误 barcode 应该如何校正。**

因为实际 sequencing 中可能出现：

```text
真实 barcode
AAACCTG...

↓ sequencing error

AAACCTA...
```

所以需要：

```text
统计 barcode
↓
根据 read support 和设置确定可用 barcode
↓
对部分 barcode 做 correction
↓
得到 permit/correction 信息
```

我使用了：

```text
--expect-cells 1000
```

意思不是强制只保留 1000 个 cells，而是告诉程序：

> 这个样本大约是 1000-cell scale。

最终进入定量矩阵的是：

```text
1231 barcodes
```

这里更准确地说是 **1231 个 quantified barcodes**。在这个 10x 数据中，它们可以作为后续 QC 的细胞对象，但在正式 QC 前不必把它们理解成已经完成质量筛选的最终 cells。

---

## 7. collate：按 corrected barcode 整理 records

运行：

```bash
alevin-fry collate
```

原始 RAD records 可能类似：

```text
Cell A
Cell C
Cell B
Cell A
Cell B
Cell A
```

collate 后逻辑上变成：

```text
Cell A
Cell A
Cell A

Cell B
Cell B

Cell C
```

也就是：

> 根据前一步得到的 barcode correction 信息，把属于同一个 corrected cell barcode 的 RAD records 聚到一起。

这样下一步 quantification 就可以高效地：

```text
逐个 barcode/cell
→ 看 gene
→ 看 UMI
```

---

## 8. quant：UMI 定量并生成 count matrix

运行：

```bash
alevin-fry quant
```

**quant = quantification，定量。**

这一步真正回答：

> 对于每个 barcode 和每个 gene，最终有多少个可计数的 RNA molecules / UMI counts？

例如：

```text
Cell A
Gene X

U1
U1
U2
```

如果直接数 read records：

```text
Gene X = 3
```

但两个 `U1` 可能来自同一个原始 RNA molecule 的 PCR amplification。

所以经过 UMI resolution / deduplication 后：

```text
Gene X = 2
```

最终生成：

```text
quants_mat.mtx
quants_mat_rows.txt
quants_mat_cols.txt
```

矩阵大小：

```text
1231 × 236823
```

---

## 9. USA mode：为什么有 236,823 列

因为 quant 使用的是三列 t2g，所以输出是 **USA mode**：

- **U = Unspliced，未剪接**
    
- **S = Spliced，已剪接**
    
- **A = Ambiguous，剪接状态无法明确归类**
    

实际 gene 数：

```text
236823 ÷ 3 = 78941 genes
```

所以原始矩阵本质上是在同一个二维矩阵里保存：

```text
1231 barcodes
×
78941 genes
×
3 count states
```

由于 `.mtx` 是二维格式，三个 state 被展开成：

```text
1231 × 236823
```

也就是说，236,823 列并不是 236,823 个不同 genes，而是：

```text
78941 genes × S/U/A 三种状态
```

---

## 10. pyroe.load_fry：matrix → AnnData

运行：

```python
adata = pyroe.load_fry(
    "pbmc_1k_quant",
    output_format="scRNA"
)
```

这里真正创建了：

```text
adata
```

`adata` 是一个 **AnnData class 的实例对象**。

pyroe 把 USA matrix 重新整理：

```text
adata.X
= Spliced + Ambiguous

adata.layers["unspliced"]
= Unspliced
```

于是得到：

```text
adata.shape
= (1231, 78941)
```

也就是：

```text
1231 observations
×
78941 variables
```

在 scRNA-seq 中：

- **obs = observations，观测对象 → cells/barcodes**
    
- **var = variables，变量/features → genes**
    

因此 AnnData 的结构是：

```text
adata
├── X
│   └── 1231 × 78941
│       Spliced + Ambiguous counts
│
├── layers["unspliced"]
│   └── 1231 × 78941
│       Unspliced counts
│
├── obs
│   └── cell/barcode metadata
│
└── var
    └── gene metadata
```

`X` 和 `unspliced` 的 shape 一样，是因为它们描述的是同一批：

```text
barcodes × genes
```

只是每个格子统计的 RNA 状态不同。

---

## 11. 保存为 H5AD

把内存中的 AnnData：

```python
adata
```

保存到磁盘：

```python
adata.write_h5ad("pbmc_1k_raw.h5ad")
```

这样退出 Python 后，不需要重新从 `.mtx` 构建 AnnData。

`.h5ad` 可以保存：

```text
X
obs
var
layers
以及其他 AnnData metadata
```

---

## 12. Gene metadata annotation：Ensembl ID → gene symbol

原始矩阵中的 gene feature 使用：

```text
ENSG00000296732.1
ENSG000001...
```

也就是 **Ensembl Gene ID**。

我没有直接把这些 ID 替换掉，而是从 GENCODE GTF 中提取：

```text
gene_id → gene_name
```

生成：

```text
gene_id_to_symbol.tsv
```

然后添加到：

```python
adata.var["gene_symbol"]
```

最终类似：

```text
Ensembl Gene ID          gene_symbol

ENSG...                  NKG7
ENSG...                  CD3D
ENSG...                  LYZ
```

这样做的好处是：

```text
Ensembl Gene ID
→ 保留作为稳定的 gene identifier

gene symbol
→ 用于阅读、生物学解释和后续 marker/QC
```

没有修改：

```text
adata.X
```

中的任何 count。

更准确地说，这一步属于 **gene metadata annotation / gene identifier annotation**，不是功能富集意义上的 gene functional annotation，也不是后面要做的 cell type annotation。

最后保存：

```text
pbmc_1k_raw_annotated.h5ad
```

---

## 最终上游产物

到这里，上游得到的是：

```text
pbmc_1k_raw_annotated.h5ad
```

包含：

```text
1231 quantified barcodes
×
78941 genes

adata.X
→ Spliced + Ambiguous counts

adata.layers["unspliced"]
→ Unspliced counts

adata.obs
→ barcode/cell metadata

adata.var
→ Ensembl Gene ID + gene_symbol
```

因此可以把整个上游浓缩成一句：

> **Raw FASTQ → reference construction → read mapping → barcode processing/correction → UMI quantification → sparse count matrix → AnnData → gene metadata annotation。**

