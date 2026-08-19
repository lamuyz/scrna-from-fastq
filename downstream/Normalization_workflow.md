# Normalization_workflow

> 当前仅记录经典的 **Shifted logarithm（移位对数归一化）流程。Scran normalization 和 analytic Pearson residuals 暂不展开。

## 1. 为什么需要归一化？

QC 后保留下来的细胞，其总 UMI counts 仍可能存在明显差异。

例如：

```text
Cell A: 4,000 UMI
Cell B: 9,000 UMI
```

这种差异不一定完全来自真实的生物学差异，也可能受到 RNA 捕获效率、逆转录效率和测序深度等 **technical sampling effects（技术采样效应）**影响。

因此，normalization（归一化）的目标是：

- 减少不同细胞之间 count depth 的技术差异；
    
- 改善 gene expression 的数值分布和方差结构；
    
- 为后续 HVG、PCA、邻居图和聚类提供更合适的输入。
    

---

## 2. 经典归一化流程

Shifted logarithm 可以拆成两步：

```text
raw counts
    ↓
normalize_total
    ↓
library-size normalized counts
    ↓
log1p
    ↓
log-normalized expression
```

对应 Scanpy：

```python
sc.pp.normalize_total(adata, target_sum=None)
sc.pp.log1p(adata)
```

---

## 3. `normalize_total`

### 3.1 作用

`normalize_total` 主要校正不同 cell 之间 **total counts（总 counts）**的尺度差异。

假设：

```text
Cell 1 total counts = 40
Cell 2 total counts = 200
Cell 3 total counts = 200
```

为了方便举例，设：

```text
target_sum = 100
```

目标是将每个 cell 整体缩放到相同的总量：

```text
Cell 1:  40 → 100
Cell 2: 200 → 100
Cell 3: 200 → 100
```

### 3.2 Size factor

**size factor（尺度因子）**表示一个 cell 相对于目标尺度需要缩放多少。

```text
size factor = cell total counts / target_sum
```

例如：

```text
Cell 1: 40 / 100  = 0.4
Cell 2: 200 / 100 = 2
Cell 3: 200 / 100 = 2
```

每个 cell 内的所有 gene 都除以同一个 size factor：

```text
normalized count = raw count / size factor
```

---

## 4. 手算示例

原始矩阵：

|Cell|Gene A|Gene B|Gene C|Gene D|Total|
|---|--:|--:|--:|--:|--:|
|Cell 1|20|10|10|0|40|
|Cell 2|120|60|20|0|200|
|Cell 3|80|40|20|60|200|

设：

```text
target_sum = 100
```

### Cell 1

```text
size factor = 40 / 100 = 0.4

Gene A = 20 / 0.4 = 50
Gene B = 10 / 0.4 = 25
Gene C = 10 / 0.4 = 25
Gene D =  0 / 0.4 = 0
```

### Cell 2

```text
size factor = 200 / 100 = 2

Gene A = 120 / 2 = 60
Gene B =  60 / 2 = 30
Gene C =  20 / 2 = 10
Gene D =   0 / 2 = 0
```

### Cell 3

```text
size factor = 200 / 100 = 2

Gene A = 80 / 2 = 40
Gene B = 40 / 2 = 20
Gene C = 20 / 2 = 10
Gene D = 60 / 2 = 30
```

归一化后：

|Cell|Gene A|Gene B|Gene C|Gene D|Total|
|---|--:|--:|--:|--:|--:|
|Cell 1|50|25|25|0|100|
|Cell 2|60|30|10|0|100|
|Cell 3|40|20|10|30|100|

需要注意：

> `target_sum` 并不是平均分配给每个 gene。

同一个 cell 内所有 gene 使用相同的缩放系数，因此原来的相对表达比例不会被改变。

---

## 5. `target_sum`

`target_sum` 表示归一化后每个 cell 的目标总量。

例如：

```python
sc.pp.normalize_total(adata, target_sum=1e4)
```

表示将每个 cell 的总表达量缩放到约：

```text
10,000
```

在本教程使用：

```python
sc.pp.normalize_total(
    adata,
    target_sum=None
)
```

时，目标尺度由数据本身决定，使用细胞 raw total counts 的中位数作为归一化尺度。

---

## 6. `log1p`

完成 `normalize_total` 后，还需要进行 logarithmic transformation（对数变换）：

```python
sc.pp.log1p(adata)
```

其计算为：

```text
log(1 + x)
```

Scanpy 默认使用 natural logarithm（自然对数）。

其中的：

```text
+1
```

称为 **pseudo-count（伪计数）**。

它没有真实的 UMI 生物学意义，主要用于避免：

```text
log(0)
```

无法计算的问题。

因此：

```text
x = 0

log(1 + 0)
= log(1)
= 0
```

---

## 7. 为什么需要 `log1p`？

完成 `normalize_total` 后，不同 gene 的表达值跨度仍然可能很大。

例如：

```text
1
10
100
1000
```

经过 `log1p`：

```text
log(2)    ≈ 0.69
log(11)   ≈ 2.40
log(101)  ≈ 4.62
log(1001) ≈ 6.91
```

因此 `log1p` 会压缩较大的表达值，改善表达数据的分布和 mean-variance relationship（均值-方差关系），避免高表达 gene 在后续 PCA 等分析中过度主导结果。

---

## 8. 两步的区别

### `normalize_total`

主要解决：

```text
不同 cell 的 total counts 不同
```

属于：

**library-size normalization（文库大小归一化）**

可以理解为 cell-level scaling（细胞层面的尺度调整）。

### `log1p`

主要解决：

```text
gene expression 数值跨度过大
```

属于：

**log transformation（对数变换）**

作用于每一个 cell × gene 的表达值。

因此经典路线可以理解为：

```text
normalize_total
    ↓
调整 cell-level count depth

log1p
    ↓
压缩 gene-expression dynamic range
```

---

## 9. 在 AnnData 中保留 raw counts

归一化前应保留原始 count matrix，避免后续需要 raw counts 时无法恢复。

例如：

```python
adata.layers["counts"] = adata.X.copy()
```

然后再进行：

```python
sc.pp.normalize_total(adata, target_sum=None)
sc.pp.log1p(adata)
```

这样：

```text
adata.layers["counts"]
```

保存 raw counts，而：

```text
adata.X
```

可以用于存储当前的 log-normalized expression matrix。

---

## 10. 当前流程中的位置

```text
FASTQ
  ↓
mapping / quantification
  ↓
count matrix
  ↓
ambient RNA correction
  ↓
QC
  ↓
normalization
  │
  ├─ normalize_total
  └─ log1p
  ↓
feature selection / HVG
  ↓
PCA
  ↓
neighbors
  ↓
UMAP / clustering
```

Normalization 属于 **preprocessing（预处理）**。

---

## 11. 当前需要记住的术语

|Term|Full name / 中文|含义|
|---|---|---|
|normalization|归一化|减少技术尺度差异|
|total counts|总 counts / 总计数|一个 cell 中所有 gene counts 的总和|
|target sum|目标总量 / 目标尺度|一个 cell 被整体缩放到的总量|
|size factor|尺度因子|描述一个 cell 需要整体缩放多少|
|library-size normalization|文库大小归一化|根据每个 cell 的总 counts 调整表达尺度|
|log transformation|对数变换|压缩高表达值的数值范围|
|pseudo-count|伪计数|`log` 前人为加入的数值，如 1|
|sparse|稀疏|矩阵中大部分值为 0|
|technical sampling effect|技术采样效应|RNA 捕获、逆转录、测序等过程带来的随机技术差异|
|overdispersion|过度离散|实际数据的波动程度超过统计模型的预期|

---

## 12. 当前采用的经典路线

现阶段只使用：

```python
sc.pp.normalize_total(adata, target_sum=None)
sc.pp.log1p(adata)
```

暂不加入：

```text
Scran normalization
Analytic Pearson residuals
```

后续再单独学习和比较不同 normalization 方法。