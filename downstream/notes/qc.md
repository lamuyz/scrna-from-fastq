# qc

## 1. QC 的目标是什么？

从 FASTQ 构建出 count matrix 后，不能直接默认每个 barcode 都对应一个高质量单细胞。

上游通常已经做过一轮：

> **cell calling，细胞判定**

主要目的是从大量 barcodes 中识别：

```text
更像真实细胞的 barcode
vs
empty droplets / background barcodes
```

但这并不意味着留下来的 barcode 都是高质量单细胞。

因此下游还需要进行 QC：

```text
count matrix
↓
low-quality cell filtering
↓
ambient RNA correction
↓
doublet detection
↓
得到更加可信的单细胞表达矩阵
```

可以简单理解为：

```text
Cell calling
→ “这个 barcode 像不像一个 cell？”

Downstream QC
→ “这个 cell 是否高质量、是否被污染、是否实际上包含多个 cells？”
```

---

# 2. Low-quality cell filtering｜低质量细胞过滤

这一部分主要利用每个 cell 的 QC metrics 判断其质量。

最重要的三个指标是：

```text
total_counts
n_genes_by_counts
pct_counts_mt
```

---

## 2.1 total_counts

`total_counts` 表示：

> 一个 cell 中所有 genes 的 count 总和。

例如：

```text
Gene A   3
Gene B   5
Gene C   0
Gene D   2
```

那么：

```text
total_counts = 3 + 5 + 0 + 2 = 10
```

对于 UMI-based scRNA-seq，可以粗略理解为：

> 这个 cell 中检测到了多少 RNA molecules / UMI counts。

---

## 2.2 n_genes_by_counts

`n_genes_by_counts` 表示：

> 一个 cell 中有多少个 gene 的 count > 0。

还是上面的例子：

```text
Gene A   3
Gene B   5
Gene C   0
Gene D   2
```

那么：

```text
total_counts = 10
n_genes_by_counts = 3
```

因此：

```text
total_counts
→ RNA 总量有多少？

n_genes_by_counts
→ 检测到了多少种 gene？
```

---

## 2.3 pct_counts_mt

`pct`：

> percentage，百分比

`mt`：

> mitochondrial，线粒体

所以：

> **pct_counts_mt = 线粒体基因 counts 占总 counts 的百分比**

例如：

```text
total_counts = 10000
mitochondrial counts = 1000
```

那么：

```text
pct_counts_mt = 10%
```

---

# 3. 为什么高 mitochondrial percentage 常提示低质量细胞？

细胞受损或死亡时，细胞质中的 RNA 可能大量流失。

```text
cell damage / cell death
↓
cytoplasmic RNA leakage
↓
普通 mRNA 大量减少
```

而线粒体由于自身具有膜结构，其 RNA 可能相对保留得更多。

因此常见：

```text
total_counts ↓
n_genes_by_counts ↓
pct_counts_mt ↑
```

注意：

> pct_counts_mt 升高并不一定代表线粒体 RNA 的绝对数量增加，更可能是普通 RNA 丢失得更多，因此 mitochondrial RNA 所占比例上升。

---

# 4. QC covariates

`covariate` （**协变量**） 在这里可以直接理解为：

> 用来描述和判断 cell quality 的指标 / 特征。

例如：

```text
total_counts
n_genes_by_counts
pct_counts_mt
```

不能只依赖其中某一个指标。

例如：

```text
total_counts 很低
```

可能是低质量细胞，也可能是一个本身 RNA 含量较低的小细胞类型。

同样：

```text
total_counts 很高
```

可能是 doublet，也可能只是某种本身体积较大的真实细胞。

所以 QC 应该 **联合观察多个指标。**

---

# 5. 如何查看 QC 指标分布？

常用：

```text
Histogram
Violin plot
Scatter plot
```

---

## 5.1 Histogram｜直方图

主要用于观察：

> 一个 QC metric 的整体分布。

例如 `total_counts`：

可以看到：

- 主体分布在哪里
    
- 有没有长尾
    
- 有没有极端值
    
- cutoff 大概应该放在哪里
    

---

## 5.2 Violin plot｜小提琴图

也是用来看数据分布。

哪里宽：

> 这个范围内 cells 多。

哪里窄：

> cells 少。

特别适合比较：

```text
不同 sample
不同 batch
不同 condition
```

之间的 QC 指标。

---

## 5.3 Scatter plot｜散点图

`scatter plot`：

> **散点图**

每一个点代表一个 cell。

例如：

```python
sc.pl.scatter(
    adata,
    x="total_counts",
    y="n_genes_by_counts",
    color="pct_counts_mt"
)
```

表示：

```text
x 轴
→ total_counts

y 轴
→ n_genes_by_counts

颜色
→ pct_counts_mt
```

这样可以同时观察多个 QC metrics。

这里的“距离”是：

> **特征 / 表达空间上的距离**

而不是细胞在组织里的物理距离。

---

# 6. Percentile-based threshold ｜基于百分位数的阈值

例如：

```text
2.5th percentile
97.5th percentile
```

把所有 cells 按某个指标排序：

```text
最低 2.5% |------ 中间 95% ------| 最高 2.5%
```

它不是固定 numerical threshold，因为具体数值取决于当前数据。

例如：

```text
Dataset A:
2.5th percentile = 200

Dataset B:
2.5th percentile = 500
```

而：

```text
pct_counts_mt < 15%
```

才是一个固定数值 cutoff。

---

# 7. MAD｜Median Absolute Deviation

MAD 全称：

> **Median Absolute Deviation** ，**中位数绝对偏差**

它用于判断：

> 哪些 observations 离主体分布特别远。

公式：

```text
MAD = median(|Xi - median(X)|)
```

---

## 7.1 MAD 的计算

假设：

```text
900
1000
1100
1200
10000
```

第一步：

```text
median = 1100
```

第二步，计算每个值离 median 多远：

```text
|900 - 1100|   = 200
|1000 - 1100|  = 100
|1100 - 1100|  = 0
|1200 - 1100|  = 100
|10000 - 1100| = 8900
```

得到：

```text
200
100
0
100
8900
```

排序：

```text
0
100
100
200
8900
```

所以：

```text
MAD = 100
```

如果用：

```text
5 MAD
```

则大致范围是：

```text
median ± 5 × MAD
```

即：

```text
1100 ± 500
=
600 ~ 1600
```

因此 `10000` 会被认为是 outlier。

---

# 8. 实际数据中如何选择 MAD？

没有：

> “scRNA-seq 必须使用 5 MAD”

这样的规则。

可以粗略理解：

```text
3 MAD → 更严格
4 MAD → 中间
5 MAD → 更宽松
```

正确思路是：

```text
先查看原始 QC distribution
↓
选择一个相对合理的 MAD
↓
先标记 outliers
↓
看被标记的 cells 是什么情况
↓
结合多个 QC metrics 判断
↓
确认 threshold
↓
再真正过滤
```

而不是：

```text
先决定必须用 5 MAD
↓
机械删除
```

QC 的一个重要原则是：

> 不要为了“数据看起来特别干净”而过度过滤。

因为一些真实 rare cell populations 本身可能就在分布尾部。

---

# 9. 给 mitochondrial / ribosomal / hemoglobin genes 打标签

例如：

```python
adata.var["mt"] = adata.var_names.str.startswith("MT-")
```

执行逻辑：

```text
adata.var_names
↓
获取所有 gene names

.str
↓
对字符串执行操作

.startswith("MT-")
↓
判断每个 gene 是否以 MT- 开头
```

例如：

```text
CD3D     False
IL7R     False
MT-CO1   True
MT-ND1   True
```

然后：

```python
adata.var["mt"] = ...
```

表示：

> 在 `adata.var` 中新增一列 `mt`，使用 True / False 标记 mitochondrial genes。

Human mitochondrial genes 常见：

```text
MT-CO1
MT-CO2
MT-ND1
MT-ND2
...
```

Mouse 命名通常需要注意大小写，例如常见 `mt-`。

---

## 9.1 Ribosomal genes

例如：

```python
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
```

常见：

```text
RPS
→ Ribosomal Protein Small subunit

RPL
→ Ribosomal Protein Large subunit
```

例如：

```text
RPS18
RPL13
```

会被标记：

```text
ribo = True
```

---

## 9.2 Hemoglobin genes

可以使用 gene naming pattern / regular expression 标记：

```text
HBA1
HBA2
HBB
...
```

然后写入：

```text
adata.var["hb"]
```

---

# 10. sc.pp.calculate_qc_metrics()

例如：

```python
sc.pp.calculate_qc_metrics(
    adata,
    qc_vars=["mt", "ribo", "hb"],
    inplace=True,
    percent_top=[20],
    log1p=True
)
```

这是 Scanpy 内置的：

> **QC metrics calculation function**

---

## 10.1 `percent_top=[20]`

表示计算：

> 一个 cell 中表达最高的前 20 个 genes，占 total counts 的百分比。

例如：

```text
total_counts = 10000

top 20 genes counts = 7000
```

那么：

```text
pct_counts_in_top_20_genes = 70%
```

如果极少数 genes 占了几乎所有 counts：

> 可能提示转录本复杂度较低或表达 profile 异常。

---

## 10.2 `log1p=True`

`log1p`：

```text
log(1 + x)
```

会生成例如：

```text
total_counts
log1p_total_counts

n_genes_by_counts
log1p_n_genes_by_counts
```

log transformation 可以：

> 压缩非常大的数值，使高度偏斜的 count distribution 更容易处理。

---

# 11. Ambient RNA｜环境 RNA

样本处理中，一些 cells 可能发生：

```text
cell lysis
↓
RNA leakage
↓
RNA 漂在细胞悬液中
```

这些游离 RNA 就叫：

> **ambient RNA，环境 RNA / 游离 RNA**

它们可能进入其他 droplets。

---

# 12. Empty droplet｜空液滴

表示：

> droplet 中没有完整 cell。

但：

```text
没有 cell
≠
没有 RNA
```

因为 ambient RNA 仍然可能进入 empty droplet。

所以：

```text
empty droplet
+
ambient RNA
↓
仍然可能产生 UMI counts
```

---

# 13. Raw matrix vs Filtered matrix

## Raw matrix

> 尚未经过 cell calling 的广泛 barcode matrix。

通常包含：

```text
真实 cells
+
empty droplets
+
background barcodes
```

## Filtered matrix

经过 cell calling 后：

> 主要保留被认为更像真实 cells 的 barcodes。

所以：

```text
raw matrix
↓
cell calling
↓
filtered matrix
```

---

# 14. SoupX｜Ambient RNA correction

SoupX 是一个：

> **R package**

用于：

> **ambient RNA contamination correction**

它不是简单删除 cells，也不是修改 gene annotation。

它校正的是：

> **gene expression counts**

---

## 14.1 Soup profile

`profile` 在这里可以理解为： **整体组成模式 / 分布模式**

`gene expression profile`：**基因表达谱**

`soup profile`：

> **环境 RNA 的组成谱 / 背景 RNA 表达谱**

例如：

```text
HBB      很多
HBA1     很多
MS4A1    一部分
CD3D     很少
...
```

这一整套 ambient RNA 的 gene composition 就是 soup profile。

---

## 14.2 Contamination fraction

`contamination fraction`： **污染比例**

表示：

> 一个 cell 的总 counts 中，估计有多少比例来自 ambient RNA。

例如：

```text
total_counts = 5000
ambient counts ≈ 500
```

那么：

```text
contamination fraction ≈ 10%
```

---

## 14.3 SoupX 的核心逻辑

首先：

```text
empty droplets
↓
没有真实 cell
↓
其中 RNA 主要反映 ambient RNA
↓
估计 soup profile
```

然后：

```text
估计每个 cell 的 contamination fraction
```

最后：

```text
soup profile
+
contamination fraction
↓
估计每个 cell × gene 中
有多少 counts 来自 ambient RNA
↓
校正 counts
```

因此 SoupX 可以理解为：

> **通过减少估计来自 ambient RNA 的 gene counts，降低环境 RNA 对真实表达矩阵的影响。**

**注意：

不是所有 genes 一起机械减小。而是根据污染模型进行针对性校正。

---

# 15. 为什么 SoupX 可以使用 clustering？

SoupX 可以利用粗略的 cluster 信息帮助判断：

> 某些 gene expression 更像真实表达还是背景污染。

例如：

```text
一个 T-cell cluster
普遍表达：
CD3D
CD3E

但基本不表达：
HBB
```

与此同时 soup profile 中：

```text
HBB 很高
```

如果少量 T cells 出现：

```text
HBB = 1
HBB = 2
```

那么这些 counts 更可能来自 ambient RNA。

所以 clustering 在 SoupX 中：

> 是辅助 contamination estimation 的信息。

不是：

> 用 clustering 直接删掉某个 cluster。

---

# 16. Doublet｜双细胞

理想情况：

```text
1 droplet
+
1 cell
+
1 bead
↓
1 barcode ≈ 1 cell
```

但有时：

```text
1 droplet
+
2 cells
+
1 bead
```

两个 cells 的 RNA 都会被同一个 Cell Barcode 标记。

最后矩阵里只有：

```text
1 barcode
```

但表达 profile 实际上是：

```text
Cell A expression
+
Cell B expression
```

这就是：

> **doublet，双细胞 / 双重捕获**

更一般的情况叫：

> **multiplet，多细胞液滴**

---

# 17. Singlet、Homotypic doublet、Heterotypic doublet

## Singlet

```text
1 droplet
+
1 cell
```

即理想的单细胞捕获。

## Homotypic doublet

两个相同或非常相似的 cell types：

```text
T cell + T cell
```

通常更难识别。

## Heterotypic doublet

两个不同 cell types：

```text
T cell + B cell
```

通常更容易检测，因为表达模式明显混合。

---

# 18. 为什么不能只根据 total_counts 检测 doublet？

Doublet 经常表现：

```text
total_counts ↑
n_genes_by_counts ↑
```

因为两个 cells 的 RNA 加在一起。

但：

```text
counts 高
≠
一定是 doublet
```

有些真实细胞本身：

- 体积较大
    
- RNA 含量较高
    
- 表达 genes 较多
    

因此需要专门的 doublet detection algorithm。

---

# 19. Artificial doublet｜人工双细胞

scDblFinder 的核心思想之一是：

> 先人为制造一些“已知是 doublet”的参考数据。

例如：

```text
真实 Cell A
+
真实 Cell B
↓
artificial doublet
```

再制造很多：

```text
A + B
A + C
B + D
C + E
...
```

这些 artificial doublets：

> 不是真实实验里的细胞，而是算法生成的参考样本。

---

# 20. 为什么 artificial doublets 可以帮助找到真实 doublets？

真实 doublet 本质上也是：

```text
两个 cell expression profiles
混合在一起
```

因此如果算法人工制造：

```text
T cell + B cell
```

那么一个真实的：

```text
T + B doublet
```

在表达特征上通常会和人工 T+B doublet 比较相似。

因此可以通过：

> “这个真实 barcode 看起来有多像人工 doublet？”

来判断 doublet 风险。

---

# 21. PCA 空间中的“距离”

这里的距离：

> **不是物理空间距离。**

而是：

> **特征 / 基因表达模式上的距离。**

原始数据每个 cell 有成千上万个 gene features：

```text
Cell
↓
Gene1
Gene2
Gene3
...
Gene20000
```

PCA：

> Principal Component Analysis，主成分分析

把高维表达信息压缩成较少的主要特征。

于是：

```text
每个 cell
↓
变成 PCA feature space 中的一个点
```

点越近：

> 整体 expression profile 越相似。

---

# 22. Nearest neighbor｜最近邻

`nearest neighbor`：

> **最近邻**

这里指：

> 在表达特征空间中与当前 cell 最相似、距离最近的一批点。

不是组织中真正挨在一起的细胞。

假设某个真实 Cell X 最近的 10 个 neighbors：

```text
1. artificial doublet
2. artificial doublet
3. real cell
4. artificial doublet
5. real cell
6. artificial doublet
7. artificial doublet
8. real cell
9. artificial doublet
10. real cell
```

其中：

```text
6 / 10
```

是 artificial doublets。

说明：

> Cell X 所处的 expression space 更像人工 doublet 区域。

因此其 doublet 嫌疑更高。

---

# 23. scDblFinder 的核心直觉

可以简化理解为：

```text
真实 cells
↓
随机组合
↓
制造 artificial doublets
↓
真实 cells + artificial doublets
放进同一个特征空间
↓
对每个真实 cell 找 nearest neighbors
↓
观察周围有多少 artificial doublets
↓
越像人工 doublet
↓
doublet probability / score 越高
```

注意：

> 实际 scDblFinder 算法比这个简化过程更复杂。

最近邻中的 artificial doublet 比例只是重要信息之一，并不是简单：

```text
6 / 10 artificial doublets
↓
doublet score = 0.6
```

最终还会综合其他 features 和 classification 过程。

---

# 24. scDblFinder 的主要输出

最后可以为每个 cell 得到：

```text
scDblFinder.score
```

表示：

> 越高越像 doublet。

以及：

```text
scDblFinder.class
```

分类为：

```text
singlet
doublet
```

这些结果可以写入：

```text
adata.obs
```

例如：

```text
barcode_1   singlet
barcode_2   singlet
barcode_3   doublet
barcode_4   singlet
```

---

# 25. Doublet 是否应该立刻删除？

不一定。

更稳妥的思路：

```text
scDblFinder
↓
先标记 doublets
↓
后续 clustering / visualization
↓
观察这些 doublets 的分布
↓
结合其他信息再次检查
↓
再决定最终过滤
```

这和前面的 MAD 思路类似：

> **先标记、检查，再删除。**

---

# 26. 完整 QC workflow

目前学习到的完整 QC 可以整理成：

```text
Count matrix
↓
计算基础 QC metrics
│
├── total_counts
├── n_genes_by_counts
├── pct_counts_mt
└── 其他 metrics
↓
查看 distribution
│
├── histogram
├── violin
└── scatter
↓
MAD / threshold
↓
标记 low-quality cells
↓
检查 outliers
↓
过滤明显低质量 cells
↓
Ambient RNA correction
│
├── raw matrix
├── empty droplets
├── soup profile
├── contamination fraction
└── SoupX
↓
得到 ambient-RNA-corrected counts
↓
Doublet detection
│
├── artificial doublets
├── PCA feature space
├── nearest neighbors
└── scDblFinder
↓
标记 doublets
↓
后续 clustering / visualization 中再次检查
↓
最终高质量 cell × gene expression matrix
```

---

# 27. 三种 QC 操作解决的是三个不同问题

这一点非常重要。

### Low-quality cell filtering

问的是：

> **这个 cell 整体质量是不是太差？**

结果通常是：

```text
整个 cell 保留 / 删除
```

---

### SoupX

问的是：

> **这个 cell 的 gene counts 中有多少来自 ambient RNA？**

结果通常是：

```text
cell 保留
但 gene expression counts 被校正
```

---

### scDblFinder

问的是：

> **这个 barcode 是不是实际上混入了两个 cells？**

结果首先是：

```text
singlet / doublet 标记
```

然后再决定是否过滤。

---
