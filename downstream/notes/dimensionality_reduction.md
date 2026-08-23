## dimensionality_reduction

### 1. 为什么需要降维

单细胞数据里，一个 cell 往往有成千上万个 gene 维度。

经过 feature selection 后，例如：

```text
20,000 genes
↓
2,000 informative genes
```

虽然已经减少了维度，但 2,000 维依然太高。

所以继续做：

```text
2,000 genes
↓ PCA
30–50 PCs
```

目的是：

> 用更少的新维度，尽量保留原始数据中最重要的 variation（变化）。

---

### 2. Feature selection 和 PCA 的区别

**Feature selection：挑 gene**

```text
20,000 genes
↓
挑出 2,000 genes
```

保留下来的维度还是 gene。

**PCA：重新组合 gene**

```text
2,000 genes
↓
PC1, PC2, PC3...
```

PC 不是 gene，而是很多 gene 的线性组合。

例如：

```text
PC1 = 0.7 × Gene A + 0.7 × Gene B
```

---

### 3. PCA 是什么

> **Principal Component Analysis**  
> 主成分分析

它会寻找数据中 variation 最大的方向。

```text
PC1 → 最大 variation
PC2 → 第二大 variation
PC3 → 第三大 variation
...
```

PC2 与 PC1 正交，可以直观理解为：

> PC2 尽量不重复 PC1 已经捕获的那套变化。

---

### 4. 为什么 Gene A / B 可以被压成一个 PC

例子：

```text
Gene A: -2 -1 0 1 2
Gene B: -2 -1 0 1 2
```

A 和 B 几乎完全同步变化，所以信息高度重复。

PCA 可以把它们组合成：

```text
PC1 ≈ Gene A + Gene B
```

于是原来两个维度：

```text
Gene A
Gene B
```

可以主要由：

```text
PC1
```

来表示。

这就是 PCA 降维的核心直觉：

> **把重复的变化模式压缩到同一个 PC 里。**

---

## 5. loading（载荷）

loading 属于： 

> **gene

它表示： 某个 gene 在某个 PC 中的权重。

例如：

```text
PC1 =
+0.40 × CD3D
+0.35 × CD3E
-0.38 × LST1
-0.32 × S100A8
```

那么：

```text
CD3D    loading = +0.40
CD3E    loading = +0.35
LST1    loading = -0.38
S100A8  loading = -0.32
```

loading 可以帮助回答：

> **这个 PC 主要由哪些 gene 驱动？**

如果正 loading 主要是 T-cell-related genes，而负 loading 主要是 myeloid-related genes，就可以推测：

> PC1 捕获的主要 variation 可能与 T-cell-like 和 myeloid-like expression difference 有关。

但注意：

> PCA 本身不知道什么是 T cell 或 myeloid，是我们后来根据 marker knowledge 解释。

---

## 6. PC score

PC score 属于：

> **cell**

它表示：

> 某个 cell 在某个 PC 轴上的坐标。

例如：

```text
Cell A PC1 score = +3.4
Cell B PC1 score = -3.2
```

说明两个 cell 位于 PC1 的不同方向。

如果 PC1 的正侧由 CD3D/CD3E 驱动、负侧由 LST1/S100A8 驱动，那么：

```text
PC1 score 很正
→ 更偏 CD3D/CD3E 那套表达模式

PC1 score 很负
→ 更偏 LST1/S100A8 那套表达模式
```

但不能机械地说：

```text
PC1 > 0 = 某类细胞
PC1 < 0 = 另一类细胞
```

PC score 是帮助区分和解释 cell variation，不是直接的 cell-type 分类标签。

---

## 7. loading 和 PC score 的区别

最重要的记忆：

```text
loading
= gene 对 PC 的权重

PC score
= cell 在 PC 上的位置
```

或者：

```text
gene --loading--> PC --score--> cell 的位置
```

---

## 8. Scanpy 里对应什么

### `adata.varm["PCs"]`

保存：

```text
gene × PC
```

也就是：

> loadings

例如：

```text
          PC1    PC2    PC3
CD3D      0.4    ...
CD3E      0.35   ...
LST1     -0.38   ...
```

---

### `adata.obsm["X_pca"]`

保存：

```text
cell × PC
```

也就是：

> PC scores

例如：

```text
          PC1    PC2    PC3
Cell1    -3.2    0.5    1.1
Cell2     2.7   -0.4    0.8
```

所以：

```text
adata.X
= cell × gene

adata.obsm["X_pca"]
= cell × PC
```

---

## 9. PC scatter plot

PC scatter plot：

> 从 `X_pca` 里选两个 PC 来画。

例如：

```text
x = PC1 score
y = PC2 score
```

每个点：

> 一个 cell

所以：

```text
PC scatter plot
= X_pca 的两个维度的二维展示
```

如果有 10 个 PC，理论上可以画很多组合，比如：

```text
PC1 vs PC2
PC1 vs PC3
PC2 vs PC3
...
```

但通常不会全部画。

最常先看：

```text
PC1 vs PC2
```

因为它们解释的 variation 最大。

---

## 10. explained variance

每个 PC 会解释原始数据中的一部分 variation。

例如：

```text
PC1 → 24%
PC2 → 18%
PC3 → 13%
...
```

前面的 PC 通常解释更多 variation。

---

## 11. scree plot / elbow

scree plot（碎石图）用于看：

> 每个 PC 解释了多少 variance。

如果：

```text
PC1–PC6
下降明显

PC7以后
开始变平
```

那么 PC6 附近可能是：

> **elbow（肘部）**

可以作为选择保留多少 PC 的参考。

但不能机械使用 elbow，还要结合：

```text
explained variance
+
loadings
+
downstream clustering / neighbors
```

---

## 12. 为什么不是只保留 PC1 和 PC2

PC1/PC2 最适合画二维 PCA 图，但后面的 PC 仍可能包含重要 biological signal。

例如：

```text
PC1 → T vs myeloid
PC2 → B-cell variation
PC3 → NK variation
PC4 → monocyte subtype
PC5 → activation state
```

所以：

> **画图用 PC1/PC2，不等于下游分析只用 PC1/PC2。**

---

# 13. t-SNE

t-SNE：

> **t-distributed Stochastic Neighbor Embedding**  
> t 分布随机邻域嵌入

更强调：

> **local structure（局部邻域结构）**

也就是：

> 哪些 cell 是彼此的近邻。

它有时能更明显展示局部小团块，但不能因为 t-SNE 分成几个岛就直接认为一定有几个真实 cell type。

---

# 14. UMAP

UMAP：

> **Uniform Manifold Approximation and Projection**  
> 统一流形近似与投影

这里的 `MAP` 不是“地图”。

标准流程通常是：

```text
PCA
↓
neighbors
↓
UMAP
```

UMAP 会尽量：

> 保留局部邻居关系，同时通常比 t-SNE 更能保留一些中尺度/整体连接结构。

例如高维结构：

```text
A —— B —— C
```

UMAP 往往更容易表现出这种连续关系。

但：

> UMAP1 / UMAP2 本身没有生物学含义。

不能说：

```text
UMAP1 越大 = 某种生物学指标越高
```

---

# 15. KNN / neighbors

KNN：

> **k-nearest neighbors**  
> k 近邻

意思：

> 对每个 cell，找距离最近的 k 个 cell。

例如：

```text
n_neighbors = 15
```

可以理解成：

> 每个 cell 重点寻找大约 15 个近邻。

---

## 16. 为什么在 PCA 空间找邻居

不是直接用 2,000 个 gene，而通常用：

```text
PC1–PC30
```

因为 PCA 已经：

- 降低维度
    
- 减少冗余
    
- 提取主要 variation
    

所以：

```text
HVG expression
↓
PCA
↓
PC space
↓
KNN
```

更稳定也更高效。

---

## 17. cell-cell neighbor graph

其中：

```text
node（节点） = cell
edge（边） = 两个 cell 是近邻
```

---

## 18. distance

`distance`：

> 两个 cell 在 PCA 空间中的几何距离。

例如：

```text
distance 小
→ cell 更相似

distance 大
→ cell 更不相似
```

---

## 19. connectivity / weight

真实的 neighbor graph 通常不是简单：

```text
连了 = 1
没连 = 0
```

而是 weighted graph（加权图）。

例如：

```text
A-B connectivity = 0.95
A-C connectivity = 0.40
```

表示：

```text
A-B 是很强的近邻关系
A-C 是较弱的近邻关系
```

而且你后来指出的理解是对的：

> **connectivity 本身就是基于 distance 和局部邻域关系构建出来的。**

执行流：

```text
PCA space
↓
计算 distance
↓
找 KNN
↓
根据局部距离关系
↓
构建 weighted neighbor graph
↓
connectivity / weight
```

---

## 20. Scanpy 里 neighbor graph 存哪里

```python
adata.obsp["distances"]
```

保存：

> 近邻距离

```python
adata.obsp["connectivities"]
```

保存：

> 近邻连接强度

它们一般是 sparse matrix（稀疏矩阵），因为每个 cell 只连接少数近邻，不需要保存所有 cell pair。

---

# 21. neighbors 和 UMAP / clustering 的关系

neighbor graph 是一个非常关键的中间结果：

```text
PCA
↓
neighbors
↓
weighted cell-cell neighbor graph
↓
├── UMAP
│   用于低维可视化
│
└── Leiden clustering
    用于聚类
```

所以：

> PCA 负责找到好的表示空间；  
> neighbors 负责建立 cell-cell 相似网络；  
> UMAP 负责把这个网络展示出来。

---

# 最终完整执行流

可以把整个降维过程压缩成：

```text
Normalized data
↓
Feature selection
↓
selected informative genes
↓
PCA
↓
得到：

1. loadings
   gene × PC
   adata.varm["PCs"]

2. PC scores
   cell × PC
   adata.obsm["X_pca"]

↓
看 explained variance / scree plot
↓
选择合适数量的 PCs
↓
KNN / neighbors
↓
计算 distance
↓
构建 weighted cell-cell neighbor graph
↓
connectivity
↓
├── UMAP
├── t-SNE（另一种可视化路线）
└── Leiden clustering
```

总结：

> **PCA 把 gene expression 压缩成 PC；loading 告诉我们 PC 由哪些 gene 驱动，PC score 告诉我们每个 cell 在 PC 空间的位置；neighbors 再利用这些 PC score 建立细胞近邻图，UMAP 把这种近邻结构投影成二维图。**