## Feature Selection

### 1. 特征选择是什么

在 scRNA-seq 里：

> **feature 基本就是 gene。**

Feature selection 的目的就是：

> 从上万甚至两万多个基因里，挑出更有信息量、能够反映细胞之间差异的基因，供后面的 PCA、邻居图、聚类等步骤使用。

也就是：

```text
所有 genes
↓
挑出 informative genes
↓
再做 PCA
```

这里的 **informative gene** 不是指“这个基因生物学上一定很重要”，而是：

> 它的表达模式能够帮助区分不同细胞。

---

## 2. 为什么不能所有 gene 都直接拿去 PCA

因为很多基因：

- 几乎全是 0；
    
- 在所有细胞中表达都差不多；
    
- 只有很弱的随机波动。

这些基因对区分细胞帮助很小，反而可能引入噪声。

所以一般先：

```text
约 20,000 genes
↓
Feature selection
↓
例如 4,000 informative genes
↓
PCA
```

---

## 3. HVG 是什么

HVG 全称：

> **Highly Variable Genes，高可变基因**

它是最常见的 feature selection 方法之一。

核心思想是：

> 找那些在不同细胞之间表达变化比较明显的基因。

传统 HVG 一般会基于处理后的表达值，看：

```text
mean
+
dispersion / variability
```

来判断哪些 genes 比预期更 variable。

所以：

> **HVG 是 feature selection 的一种方法，不等于 feature selection 本身。**

---

## 4. 教程主要使用的方法

这篇教程重点介绍的是：

> **deviance-based feature selection**
> 
> 基于偏差度的特征选择

它和 HVG 的最终目的差不多：

> 都是在找 informative genes。

但判断方法不同。

### HVG

大致是：

```text
normalized / transformed expression
↓
mean + dispersion
↓
判断 gene 是否高度 variable
```

### Deviance-based 方法

是：

```text
raw counts
↓
建立简单的 null model
↓
看实际表达模式偏离模型有多严重
↓
deviance 高 → 更 informative
```

---

## 5. 为什么教程喜欢 deviance 方法

一个重要优点是：

> **它可以直接基于 raw counts 做 feature selection。**

因此不像传统 HVG 那么依赖前面的处理：

```text
normalization
log transformation
pseudo-count
target_sum
```

> deviance-based feature selection 可以减少前面 preprocessing 方法对“哪些基因被选中”的影响。

---

## 6. Deviance 应该怎么直观理解

可以简单理解为：

> **某个 gene 的真实表达模式，和“普通、相对稳定的表达模式”差得有多远。**

例如：

```text
Gene A:
5, 5, 5, 5, 5
```

所有细胞都一样。

所以：

```text
很符合简单模型
↓
deviance 小
↓
信息量低
```

而：

```text
Gene B:
0, 0, 0, 20, 25
```

明显存在：

```text
Cell1–3：不表达
Cell4–5：高表达
```

所以：

```text
实际模式明显偏离简单模型
↓
deviance 大
↓
更 informative
```

因此最后按 deviance 排序，选择最高的一部分 genes。

---

## 7. 怎么选 gene（ deviance 方法）

每个 gene 算完 deviance 后，大概会得到：

```text
Gene A → deviance 小
Gene B → deviance 大
Gene C → deviance 中等
...
```

然后：

```text
按 deviance 排序
↓
选择 top genes
```

教程示例里选择：

```text
top 4000 genes
```

然后在：

```python
adata.var
```

中加入类似：

```text
binomial_deviance
highly_deviant
```

其中：

```text
highly_deviant = True
```

表示这个 gene 被选中了。

注意：

> 这里主要是在给 gene 打标签，并不一定马上把其他 genes 从 AnnData 中删除。

---

## 8. dispersion 是什么

- variance = 方差
    
- dispersion = 离散度 

两者都描述数据有多分散，但 dispersion 更强调：

> **相对于基因自身平均表达量，它到底有多大的变异。**

因为一个高表达基因天然可能有更大的绝对 variance，所以只看 variance 不够公平。

---

## 9. 教程最后又运行 `sc.pp.highly_variable_genes()`

这里并不代表教程突然改回 HVG 作为主要 feature selection 方法。

教程主要是想利用这个函数计算：

```text
means
dispersions
dispersions_norm
```

然后画：

```text
x = mean
y = dispersion
```

再把：

```text
highly_deviant = True
```

的 genes 标出来。

所以真正选 gene 的依据仍然是：

> **deviance**

不是  `highly_variable`

---

## 10.  上一步画的 mean–dispersion 图在检查什么

它主要有两个作用。
### 第一：观察 highly deviant genes 的统计特征

比如看：

> deviance 选中的 genes，是不是大部分也有比较高的 dispersion？

也就是把它们放到传统 mean–dispersion 框架中看看分布。

### 第二：直观比较 deviance 和传统 HVG 思路

如果很多 highly deviant genes 都落在 high-dispersion 区域，说明：

> 两种方法虽然统计思想不同，但经常会挑中相似类型的 informative genes。

不过这只是：

> **直观比较**

如果真的想计算两种方法的 gene 重合度，应该直接比较：

```python
adata.var["highly_deviant"]
adata.var["highly_variable"]
```

例如计算两者同时为 `True` 的 gene 数量。

---

## 11. Feature selection 和 PCA 的关系
### Feature selection

是在原来的 gene 中做筛选：

```text
20,000 genes
↓
4,000 genes
```

这些仍然是真实 genes。

### PCA

是进一步重新构造新的维度：

```text
4,000 genes
↓
PC1
PC2
PC3
...
```

所以：

> **Feature selection 是挑 gene；PCA 是把选出的 gene 压缩成新的 principal components。**

---

