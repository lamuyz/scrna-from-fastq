# Annotation（细胞类型注释）

## 1. Annotation 的目的

单细胞 RNA 测序中，clustering 只能把表达模式相似的 cells 分到同一个 cluster。

例如：

```text
cluster 0
cluster 1
cluster 2
```

但：

```text
cluster ≠ cell type
```

Leiden clustering 并不知道：

```text
cluster 0 是什么细胞？
```

Annotation 的目的：

> 根据 gene expression pattern、marker genes 和已有生物学知识，为 cell 或 cluster 赋予 cell type identity。

例如：

```text
cluster 0 → T cell

cluster 1 → B cell

cluster 2 → Monocyte
```

因此：

```text
Clustering

↓  

哪些cells表达模式相似？


Annotation

↓

这些cells在生物学上是什么？
```

---

## 2. Annotation 在单细胞分析流程中的位置

基础流程：

```text
Raw counts

↓

Quality control

↓

Normalization

↓

Feature selection

↓

PCA

↓

Neighbor graph

↓

Leiden clustering

↓

UMAP

↓

Cell type annotation
```

做到 annotation 后，每个 cell 可以同时具有：

```text
cell ID
+
cluster label
+
cell type label
```

---

## 3. Cellular phenotype 和 cell identity

### 3.1 Cellular phenotype

> 细胞表型

指细胞表现出来的特征，例如：

- gene expression
    
- protein expression
    
- cellular function
    
- morphology
    
- cell state
    

scRNA-seq annotation 主要利用：

```text
gene expression pattern
```

判断 cellular phenotype。

---

### 3.2 Cell identity

> 细胞身份 / 细胞类型身份

例如：

```text
T cell

B cell

NK cell

Monocyte
```

不同 cell type 执行不同功能，因此会表达不同的 gene expression program。

所以：

```text
cell identity

↓

cellular function

↓

gene expression program

↓

marker gene pattern
```

可以利用 gene expression 推断 cell identity。

---

## 4. Marker gene

> 标志基因

更准确的理解是：

> 对识别某种 cell type 有信息量的 gene。

注意：

```text
marker gene

≠

只在这一种cell type中表达的gene
```

很多 cell type 并不存在唯一 marker。

因此 annotation 通常需要：

```text
marker combinations
```

甚至结合：

```text
expression thresholds
```

进行判断。

---

## 5. Positive marker 和 Negative marker

### 5.1 Positive marker

Positive marker：

> 阳性 marker / 正向 marker

指：

> 某种 cell type 应该表达的 gene。

例如：

```text
Marker A +

↓

支持某种cell identity
```

---

### 5.2 Negative marker

Negative marker：

> 阴性 marker / 负向 marker

指：

> 某种 cell type 应该不表达或低表达的 gene。

例如：

```text
Cell type X

Positive:
A +
B +
C +

Negative:
D -
```

如果：

```text
A +
B +
C +
D -
```

整体支持该 annotation。

但如果：

```text
A +
B +
C +
D也很高
```

则说明：

```text
该annotation可能存在矛盾
```

因此 annotation 不只是看：

```text
应该表达什么？
```

还需要看：

```text
不应该表达什么？
```

---

## 6. 为什么通常先 clustering 再人工 annotation？

单个 cell 容易受到：

- sparsity
    
- transcript sampling
    
- technical noise
    

影响。

例如：

```text
一个真正的T cell

↓

刚好没有检测到CD3D

↓

单独判断容易出错
```

而一个 cluster：

```text
多个表达相似的cells

↓

一起观察marker expression

↓

随机缺失影响减弱
```

因此 cluster-level annotation：

```text
more robust to noise
```

---

## 7. Manual annotation（人工注释）

> 人工注释

核心：

```text
gene expression

+

known biological knowledge

↓

人工判断cell type
```

主要有两条路线：

```text
① Known markers

↓

cluster annotation
```

和：

```text
② Cluster DE genes

↓

cell type annotation
```

两种路线可以互相验证。

---

### 7.1 路线一：Known markers → cluster annotation

先根据文献获得 cell-type marker list。

例如：

```text
T cell

→ CD3D
→ CD3E
→ TRBC1
→ TRBC2
```

```text
B cell

→ MS4A1
→ CD79A
→ CD79B
→ CD74
```

```text
NK cell

→ NKG7
→ GNLY
→ KLRD1
```

然后：

```text
known marker set

↓

检查这些markers在哪些clusters表达

↓

找到最符合该marker pattern的cluster

↓

赋予cell type
```

---

### 7.2 路线二：Cluster DE genes → cell type

第二条人工注释路线：

```text
cluster

↓

找这个cluster的DE genes

↓

查看排名靠前的genes

↓

查这些genes对应的biology

↓

推断cell type
```

例如：

```text
cluster 9

↓

KLRD1
PRF1
KLRF1
...

↓

这些genes与NK cell biology有关

↓

cluster 9可能是NK cell
```

---

### 7.3 Marker gene list 从哪里来？

Marker 最好来自：

```text
同物种

+

同组织

+

RNA-based literature
```

例如 human PBMC 数据：

优先寻找：

```text
human PBMC

human peripheral blood

immune cell scRNA-seq
```

相关文献中的 marker genes。

注意：

> FACS 中好用的 protein marker，不一定在 transcriptomic data 中同样好用。

并且：

```text
一个dataset中的好marker

≠

另一个dataset中一定也是好marker
```

---

## 8. UMAP 在 annotation 中的作用

可以把某个 marker 的 expression 映射到 UMAP。

例如：

```python
sc.pl.umap(
    adata,
    color="CD3D"
)
```

作用：

> 查看表达该 marker 的 cells 主要分布在哪里。

例如：

```text
CD3D high cells

↓

主要集中在cluster 4
```

说明：

```text
cluster 4
可能具有T-cell identity
```

但需要记住：

```text
UMAP上离得远

≠

一定是不同cell type
```

以及：

```text
某个marker在UMAP某区域高表达

≠

该区域自动等于某种cell type
```

UMAP 在 annotation 中主要用于：

```text
visual inspection
```

即：

> 可视化检查 marker expression 的空间分布。

---

## 9. Dotplot

普通 scatter plot：

```text
每个点

=

一个cell / 一个观测
```

而 annotation dotplot：

```text
x轴

=

marker genes


y轴

=

clusters
```

每一个点表示：

```text
一个gene

×

一个cluster
```

的汇总表达信息。

通常：

```text
点大小

↓

该cluster中表达这个gene的cell比例
```

颜色：

```text
↓

该gene在这个cluster中的平均表达水平
```

因此 dotplot 可以：

> 比较一组 marker genes 在不同 clusters 中的表达模式。

---

## 10. UMAP 和 Dotplot 的区别

UMAP：

```text
某个marker

↓

主要在哪里表达？
```

Dotplot：

```text
一组markers

↓

在不同clusters中的表达模式是什么？
```

两者作用不同，但可以互相补充。

---

## 11. DE 和 DE ranking

> Differentially Expressed Genes

DE ranking：

> 差异表达基因排名

例如对于 cluster 0：

```text
cluster 0

vs

其它所有cells
```

对 genes 做差异表达统计，并按照统计结果排序。

目的：

> 找到当前 cluster 最值得关注的候选 marker genes。

---

### 11.1 `rank_genes_groups()`

Scanpy：

```python
sc.tl.rank_genes_groups(
    adata,
    groupby="leiden",
    method="wilcoxon"
)
```

例如：

```text
cluster 0

↓

与其它所有cells比较

↓

Gene A
Gene B
Gene C
...
```

每个 cluster 都会得到一个 ranked gene list。

---

### 11.2 Wilcoxon rank-sum test

代码：

```python
method="wilcoxon"
```

使用：

> Wilcoxon rank-sum test，Wilcoxon 秩和检验

它属于：

> 非参数统计检验

主要比较：

```text
当前cluster

vs

其它cells
```

某个 gene 的 expression distribution 是否存在差异。

---

### 11.3. Score 怎么理解？

`rank_genes_groups()` 可以得到：

```text
score
```

对于 Wilcoxon：

```text
较大的正score

↓

该gene的expression distribution

更倾向于在当前cluster中偏高
```

因此：

```text
score越高

↓

作为当前cluster高表达candidate marker的统计证据越强
```

但是：

```text
score高

≠

一定是好的cell-type marker
```

因为高排名 gene 可能是：

- stress gene
    
- cell-cycle gene
    
- interferon-response gene
    
- mitochondrial gene
    
- ribosomal gene
    

因此：

```text
cluster-discriminating gene

≠

cell-type-specific marker
```

---

### 11.4 p-value 不是判断 marker 的唯一标准

不能简单：

```text
p-value最小

↓

最好的marker
```

p-value 主要回答：

> 当前 cluster 和其它 cells 的表达差异是否具有统计显著性。

而判断 marker 还需要考虑：

```text
expression direction

+

effect size

+

cluster specificity

+

表达该gene的cell比例

+

gene biological function
```

所以 annotation 不能只看 p-value。

---

### 11.5 Double-dipping

Double-dipping：

> 同一份数据被用于两次相关分析。

例如：

```text
expression data

↓

先用于clustering

↓

得到cluster
```

然后：

```text
同一份expression data

↓

比较这些cluster之间的DE
```

问题：

cluster 本身就是根据 expression difference 分出来的。

因此再用同一份数据问：

```text
这些clusters是否存在expression difference？
```

很容易得到非常显著的 p-value。

所以：

```text
很小的p-value

≠

证明cluster一定是真实生物学群体
```

DE analysis 在这里主要用于：

```text
寻找candidate markers

+

解释cluster biology
```

---

### 11.6 DE gene 不等于 cluster-specific marker

某个 gene 可能：

```text
cluster 3中高表达
```

但：

```text
cluster 4
cluster 5
cluster 7

也大量表达
```

它仍可能是：

```text
DE gene
```

但不是很好的：

```text
cluster-specific marker
```

因此可以进一步筛选：

```text
当前cluster中

较多cells表达

+

其它clusters中

较少cells表达
```

获得：

> 更具有 cluster specificity 的 candidate markers。

---

## 12. Manual annotation 的完整流程

人工注释可以整理为：

```text
cluster

↓

DE ranking

↓

candidate markers

↓

检查cluster specificity

↓

查gene biology

↓

提出cell type hypothesis

↓

检查known markers

↓

检查positive / negative markers

↓

UMAP + dotplot验证

↓

赋予cell type label
```

如果证据不足：

```text
NK-like

Myeloid-like

Unknown

NK cells (?)
```

都比强行给一个过细的 annotation 更合理。

---

## 13. Automated annotation（自动注释）

核心不是：

```text
算法自己发现了cell type
```

而是：

> 自动利用已有 annotation knowledge，把已有知识转移到新的数据上。

人工注释：

```text
expression

↓

human biological knowledge

↓

cell type
```

自动注释：

```text
expression

↓

classifier / reference knowledge

↓

predicted cell type
```

---

### 13.1 Marker-based classifier

一种自动注释方法：

```text
预先定义cell-type markers

↓

算法检查marker expression

↓

自动分类
```

例如：

- Garnett
    
- CellAssign
    

本质：

> 把人工 marker-based annotation 的部分判断规则自动化。

---

### 13.2 利用更多 genes 的 classifier

另一类方法不只使用少量 marker genes。

而是：

```text
几千个gene expression values

↓

trained classifier

↓

predicted cell type
```

例如：

> CellTypist

模型提前在已有 label 的细胞上训练：

```text
gene expression pattern

↓

known cell type
```

然后对新数据：

```text
unknown cell

↓

trained model

↓

predicted label
```

---

### 13.3. Reference mapping

另一类自动注释方法：

> Annotation by mapping to a reference

涉及：

```text
Reference
```

和：

```text
Query
```

---

### 13.4 Mapping 和 Label transfer

基本逻辑：

```text
query cells

↓

映射到reference建立的表示空间

↓

寻找附近的reference cells
```

之后进行：

```text
label transfer
```

例如：

```text
Query Cell X

附近reference cells：

T
T
T
T
B

↓

大多数是T cell

↓

Query Cell X → T cell
```

---

### 13.5 Joint embedding

Joint embedding：

> 联合嵌入 / 联合表示空间

意思：

> 把 reference 和 query 放进同一个低维表示空间。

例如：

```text
Reference cells

↓

latent coordinates
```

```text
Query cells

↓

同一套latent coordinates
```

这样才能比较：

```text
query cell

附近有哪些reference cells？
```

---

### 13.6 KNN label transfer

KNN：

> K-Nearest Neighbors

中文：

> K近邻

例如一个 query cell 最近的 15 个 reference cells：

```text
14 T cells

1 B cell
```

那么可以：

```text
Query Cell

↓

T cell
```

核心假设：

> 相同 cell type 在合理的 transcriptomic representation 中应该彼此接近。

---

## 14. Coarse annotation 和 Fine annotation

Coarse annotation：

> 粗粒度注释

例如：

```text
T cell
B cell
NK cell
Monocyte
```

Fine annotation：

> 细粒度注释

例如：

```text
CD4 T cell
CD8 T cell

Naive B
Memory B

Classical monocyte
Non-classical monocyte
```

可以理解：

```text
coarse

↓

大类
```

```text
fine

↓

更具体的subtype
```

---

## 15. Normalization 和 Batch correction

必须区分：

```text
Normalization

≠

Batch correction
```

Normalization 主要处理：

- library size
    
- sequencing depth
    
- capture efficiency
    

例如之前使用：

```text
scran normalization
```

属于：

```text
normalization
```

不是 batch correction。

Batch correction：

> 批次校正

主要处理：

- 不同实验批次
    
- 不同 donor
    
- 不同时间
    
- 不同实验平台
    

造成的系统性技术差异。

Reference mapping 时：

```text
reference

+

query
```

通常来自不同数据源，因此需要考虑 batch effect。

---

### 16. Manual、Classifier 和 Reference mapping

## Manual annotation

```text
cluster

↓

markers

↓

人工解释biology

↓

cell type
```

---

## Classifier

```text
cell expression

↓

trained model

↓

predicted cell type
```

---

## Reference mapping

```text
query

↓

映射到reference

↓

寻找reference neighbors

↓

label transfer
```

三者共同点：

> 都是利用 gene expression 推断 cell identity。

---

## 17. Annotation 的核心原则

需要记住：

```text
1. Annotation不是看到一个marker就命名。
```

```text
2. Cell type identification依赖一组一致的marker expression pattern。
```

```text
3. Marker gene不一定是某种cell type独有的gene。
```

```text
4. DE genes只是candidate markers，
不等于cell-type-specific markers。
```

```text
5. UMAP用于观察marker主要在哪里表达。
```

```text
6. Dotplot用于比较一组markers在不同clusters中的表达模式。
```

```text
7. Manual annotation可以从两个方向进行：

known markers
→ clusters

或者

cluster DE genes
→ cell type
```

```text
8. Positive marker和negative marker都可以作为annotation证据。
```

```text
9. Automated annotation利用已有知识，
不是自动发现新的生物学真相。
```

```text
10. 自动annotation仍然需要结合known markers人工验证。
```

```text
11. Annotation允许存在不确定性。

证据不足时不要强行给出过细的cell-type label。
```

---
