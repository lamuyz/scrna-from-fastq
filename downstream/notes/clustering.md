# Clustering（细胞聚类）

## 1. Clustering 的目的

单细胞 RNA 测序数据中，每个 cell 都有一个 gene expression profile。

目标：

> 根据细胞之间的 transcriptomic similarity（转录组相似性），将表达模式相似的细胞分到同一个 cluster。

注意：
```

cluster ≠ cell type

```

Clustering 只负责：
```

哪些 cell 表达模式相似？

```

而不是直接判断：
```

这个 cluster 是什么细胞？

```

例如：
```

Clustering:

cell1  
cell2  
cell3

↓

cluster 0

```

后续需要通过 marker genes 判断：
```

cluster 0 → T cell

```

这个过程属于 cell type annotation。

---

# 2. Clustering 在单细胞分析流程中的位置

完整流程：
```

Raw counts  
↓  
Quality control  
↓  
Normalization  
↓  
Feature selection (HVG)  
↓  
Dimensionality reduction (PCA)  
↓  
Neighbor graph construction  
↓  
Clustering (Leiden)  
↓  
UMAP visualization  
↓  
Cell type annotation

```

---

# 3. 为什么不直接对 gene expression matrix 聚类？

原始单细胞数据：
```

cells × genes

```

例如：
```

10000 cells × 20000 genes

```

每个 cell 是一个：
```

20000维向量

```

直接计算 cell similarity 会受到：

- 高维度
- dropout
- noise
- technical variation

影响。

因此通常：
```

gene expression matrix

↓

HVG selection

↓

PCA

↓

neighbor graph

↓

clustering

````

---

# 4. Neighbor graph（细胞邻接图）

## 4.1 什么是 neighbor graph？

中文：

- 细胞邻接图
- 细胞近邻图


它是一种 graph（网络结构）。

组成：

| 图元素 | 单细胞含义 |
|---|---|
| Node（节点） | cell |
| Edge（边） | 两个cell之间的连接 |
| Weight（权重） | cell之间的相似程度 |

---

## 4.2 `sc.pp.neighbors()`

Scanpy：

```python
sc.pp.neighbors(
    adata,
    n_neighbors=15,
    n_pcs=30
)
````

作用：

> 根据 PCA 空间中的 cell 坐标，计算 cell-cell similarity，并建立 neighbor graph。

注意：

```
sc.pp.neighbors()
不是画图函数
```

它只是创建并保存 graph。

结果保存：

```python
adata.obsp
```

例如：

```python
adata.obsp["connectivities"]
```

里面保存：

cell-cell connectivity matrix。

---

## 4.3 Neighbor graph 和 UMAP 的区别

错误理解：

```
UMAP = neighbor graph
```

正确：

```
PCA

↓

neighbor graph
(sc.pp.neighbors)

↓

Leiden clustering

↓

UMAP visualization
```

UMAP：

只是把高维结构压缩到二维进行展示。

---

# 5. Leiden clustering

## 5.1 Leiden是什么？

Leiden：

> 一种基于 graph 的聚类算法（graph-based clustering algorithm）。

它不是直接对 gene expression 聚类。

输入：

```
weighted neighbor graph
```

输出：

```
cluster labels
```

例如：

|cell|Leiden label|
|---|---|
|cell1|0|
|cell2|0|
|cell3|1|
|cell4|2|

---

## 5.2 Community detection

Leiden 的核心思想：

> 在 cell-cell graph 中寻找连接紧密的社区。

例如：

```
        cell2
       /
cell1 ---- cell3


        cell5
       /
cell4 ---- cell6
```

Leiden认为：

```
cell1 cell2 cell3

属于一个community


cell4 cell5 cell6

属于另一个community
```

这些community就是cluster。

---

# 6. k-means 和 Leiden 的区别

## k-means

中文：

> K均值聚类算法

特点：

需要提前指定：

```
k = cluster数量
```

例如：

```python
k=10
```

问题：

1. 不知道真实cluster数量
    
2. 假设cluster形状比较规则
    
3. 不适合复杂单细胞结构
    

---

## Leiden

特点：

- 基于graph
    
- 不需要提前指定cluster数量
    
- 适合单细胞复杂结构
    

因此scRNA-seq中更常用：

```
Leiden clustering
```

---

# 7. Resolution 参数

## 7.1 含义

resolution：

中文：

> 聚类分辨率

Scanpy：

```python
sc.tl.leiden(
    adata,
    resolution=1
)
```

作用：

控制cluster划分的细致程度。

---

## 7.2 resolution 对结果的影响

低 resolution：

```
resolution=0.2
```

cluster较少：

```
immune cell
```

---

高 resolution：

```
resolution=2
```

cluster更多：

```
T cell

CD4 T cell

CD8 T cell

Memory T cell
```

---

规律：

```
resolution ↑

↓

cluster数量 ↑
```

但是：

cluster越多不代表越准确。

---

# 8. n_neighbors 参数

代码：

```python
sc.pp.neighbors(
    adata,
    n_neighbors=15
)
```

含义：

> 每个cell连接多少个最近邻。

---

小 n_neighbors：

例如：

```
n_neighbors=5
```

特点：

- 更关注局部结构
    
- 更容易发现细小亚群
    
- 可能过度分群

---

大 n_neighbors：

例如：

```
n_neighbors=50
```

特点：

- 更关注整体结构
    
- cluster更稳定
    
- 小亚群可能被合并
    

常用：

```
10-30
```

---

# 9. n_pcs 参数

代码：

```python
sc.pp.neighbors(
    adata,
    n_pcs=30
)
```

含义：

> 使用多少个 PCA 主成分计算 cell similarity。

例如：

PCA得到：

```
PC1
PC2
...
PC50
```

设置：

```
n_pcs=30
```

表示：

使用：

```
PC1-PC30
```

计算距离。

---

为什么不用全部PC？

因为后面的PC可能包含：

- noise
    
- technical variation
    

通常根据：

```
PCA variance explained
```

选择。

---

# 10. 如何比较不同 resolution？

实际分析不会只运行一次。

例如：

```
resolution=0.2

resolution=0.5

resolution=1

resolution=2
```

比较：

## ① cluster数量

例如：

|resolution|cluster数量|
|---|---|
|0.2|5|
|0.5|8|
|1|15|
|2|30|

---

## ② cluster是否稳定

如果稍微改变参数：

cluster完全变化：

说明：

```
cluster不稳定
```

---

## ③ 生物学合理性

需要结合后续：

marker genes

判断：

好的：

```
T cell

↓

CD4 T
CD8 T
```

有生物意义。

不好的：

```
T cell A

T cell B

T cell C
```

但marker完全一样。

可能只是过度聚类。

---

# 11. UMAP 中的 color 参数

例如：

```python
sc.pl.umap(
    adata,
    color="leiden"
)
```

含义：

不是：

```
Leiden是一种颜色
```

而是：

```
使用 adata.obs["leiden"] 中保存的cluster标签进行颜色显示
```

例如：

adata.obs：

|cell|leiden|
|---|---|
|cell1|0|
|cell2|0|
|cell3|1|

UMAP：

```
cluster0 → 一个颜色

cluster1 → 一个颜色
```

颜色只是展示方式。
