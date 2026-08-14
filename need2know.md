
## 1. scRNA-seq overall workflow

一个典型的  **scRNA-seq** 流程可以简化为：

```text
Cells
↓
mRNA
↓
cDNA library preparation
↓
Sequencing
↓
FASTQ
↓
Barcode / UMI processing
↓
Read mapping
↓
Gene assignment
↓
UMI deduplication
↓
Cell calling
↓
Gene × cell count matrix
↓
QC
↓
Normalization
↓
Dimensionality reduction
↓
Clustering
↓
Marker gene identification
↓
Cell type annotation
```

可以粗略分成两部分：

```text
上游：
FASTQ → count matrix

下游：
count matrix → QC → clustering → annotation → biological interpretation
```

---

## 2. Sequencing library

测序文库不是某一条 DNA，而是：

> 一大批已经经过加工、能够被测序仪读取的 DNA 分子的集合。

在 scRNA-seq 中，大致过程是：

```text
mRNA
↓
reverse transcription
↓
cDNA
↓
加入 barcode、UMI、adapter 等结构
↓
大量加工好的 DNA molecules
↓
sequencing library
```

原始生物分子是 RNA，但真正送进测序仪读取的是经过处理后的 DNA 分子。

---

## 3. cDNA （**complementary DNA**）

cDNA 是以 RNA 为模板，通过：

逆转录合成得到的 DNA。

在 scRNA-seq 中，cDNA 的序列来源于原来的 mRNA，因此可以用来判断：

> 这条 RNA 来自哪个 gene。

---

## 4. Cell barcode

它回答的问题是：

> 这条 RNA 来自哪个 cell？

例如：

```text
Cell A → barcode AAA...
Cell B → barcode CCC...
```

同一个细胞中的 RNA molecule 会被赋予相同的 cell barcode。

不要和 **sample index，样本索引**混淆。

```text
Sample index
→ 区分不同样本

Cell barcode
→ 区分同一个样本里的不同细胞
```

---

## 5. UMI

**UMI = Unique Molecular Identifier，唯一分子标识符**

它回答的问题是：

> 这是这个细胞里的哪个原始 RNA molecule？

例如：

```text
Cell A + CD3D + UMI U1
Cell A + CD3D + UMI U1
Cell A + CD3D + UMI U2
```

虽然一共有 3 条 reads，但不同的 UMI 只有：

```text
U1
U2
```

所以真正代表的原始 RNA molecule 数量更接近：

```text
2
```

因此：

```text
Cell A 的 CD3D count = 2
```

为什么要这样做？

因为一个原始 cDNA molecule ，经过PCR 之后可能被扩增成很多 copies。

如果直接数 reads，就可能把 PCR 复制产生的重复 reads 当成多个原始 RNA molecule。

UMI 可以帮助进行：

**UMI deduplication，UMI 去重**

---

## 6. R1 and R2

R1 和 R2 是：

**paired-end reads，双端测序 reads**

它们来自同一个 sequencing library molecule，但读取的是不同区域。

在典型的 10x 3' Gene Expression 数据中，可以先这样理解：

```text
[Cell barcode][UMI] -------- [cDNA]
       ↑                       ↑
      R1                      R2
```

---

### R1

R1 主要包含：

```text
Cell barcode
+
UMI
```

它回答：

```text
Which cell?
→ 哪个细胞？

Which RNA molecule?
→ 哪个原始 RNA 分子？
```

对于常见的 **10x Chromium 3' v3 chemistry**，可以粗略理解为：

```text
R1
├── 16 bp cell barcode
└── 12 bp UMI
```

所以总长度约为：

```text
28 bp
```

---

### R2

R2 主要读取：

**cDNA sequence，cDNA 序列**

它的作用是：

> 判断这条 RNA 来自哪个 gene。

例如：

```text
R2:
ACTGCTAGCT...
↓
mapping
↓
reference
↓
gene assignment
↓
CD3D
```

注意：

R2 本身并不会直接写着 “CD3D”。

它只有一串 DNA sequence，需要经过 mapping 才能确定来源。

---

## 7. 为什么 R2 能读到基因信息

R2 不是“碰巧”读到了 gene sequence。

这是由：**文库设计**提前决定的。

建库时，会人为把最终分子设计成类似：

```text
[barcode][UMI] -------- [cDNA insert]
     ↑                        ↑
    R1                       R2
```

R2 sequencing primer 的位置被设计成：

> 从这里开始读，就会进入 cDNA insert。

因此 R2 读取的是来源于原始 mRNA 的序列。

---

## 8. FASTQ

FASTQ 文件用于保存 sequencing reads。

每条 FASTQ record 有 4 行：

```text
1. Read identifier
2. Sequence
3. +
4. Quality scores
```

例如：

```text
@read001
ACTGCTAGCT...
+
FFFFFFFFFF...
```

---

## 9. I1

**I1 = Index Read 1

I1 通常包含：**sample index，样本索引**

它回答的问题是：

> 这条 read 属于哪个 sample？

例如一次 sequencing run 中混合多个样本：

```text
Sample A → index AAAA
Sample B → index CCCC
Sample C → index GGGG
```

测序完成后，再根据这些 index 将 reads 拆回对应样本。

这个过程叫**样本拆分**

---

## 10. Sample index vs Cell barcode

这是很容易混淆的一组概念。

```text
I1 / sample index
↓
Which sample?
→ 哪个样本？

R1 / cell barcode
↓
Which cell?
→ 哪个细胞？

R1 / UMI
↓
Which RNA molecule?
→ 哪个原始 RNA molecule？

R2
↓
Which gene?
→ 哪个 gene？
```

可以理解成层级关系：

```text
Sample
└── Cell
    └── RNA molecule
        └── Gene identity
```

---

## 11. Sequencing lane

**Sequencing lane = 测序通道 / 测序泳道**

例如：

```text
L001 = Lane 1
L002 = Lane 2
```

同一个 biological sample 可以分布在多个 lane 中进行 sequencing。

分析时，可以把不同 lane 产生的 reads 再合并起来。

---

## 12.FASTQ 文件名

例如：

```text
pbmc_1k_v3_S1_L001_R2_001.fastq.gz
```

可以拆成：

```text
pbmc_1k_v3
→ sample name，样本名

S1
→ Sample 1

L001
→ Lane 1，第 1 测序通道

R2
→ Read 2

001
→ file/chunk number，文件分块编号

.fastq.gz
→ gzip-compressed FASTQ，gzip 压缩的 FASTQ
```

---

## 13. Reference genome

假设 R2 得到：

```text
ACTGCTAGCTGAC...
```

电脑无法只看这一串字符就知道：

> 这是 CD3D。

必须进行：**mapping / alignment**

大致流程：

```text
R2 sequence
↓
mapping
↓
reference genome
↓
genomic location
↓
gene annotation
↓
gene identity
```

---

## 14. GRCh38｜人类参考基因组版本

**GRCh38 = Genome Reference Consortium Human Build 38，人类参考基因组第 38 版**

它是一个 human reference genome assembly。

GRCh38 不是单独某一个具体文件。

基于 GRCh38 可以生成不同软件需要的 reference/index：

```text
GRCh38 FASTA + annotation
│
├── Cell Ranger reference
├── STAR index
├── Salmon index
└── splici reference
```

---

## 15. FASTA

FASTA 用于保存 nucleotide 或 protein sequence。

在人类 reference genome 中，例如：

```text
>chr1
ACTGCTGACT...
```

它告诉软件：

> chr1 实际上是什么 DNA sequence。

---

## 16. GTF (Gene Transfer Format)

GTF 是常见的**基因注释文件格式**

它会描述：

```text
gene
transcript
exon
chromosome
start position
end position
strand
gene_id
transcript_id
gene_name
```

例如它可以告诉软件：

```text
chr1 的某个区域
↓
属于某个 gene
↓
这个区域还是 exon
↓
这个 exon 属于某个 transcript
```

---

## 17. Multi-mapping

如果一条 read sequence 可以匹配多个不同位置，就不能简单确定它到底来自哪里。

例如：

```text
Read X
├── matches Gene A
└── matches Gene B
```

这种 read 的 gene assignment 会更加模糊。

所以：

> 并不是每一条 R2 最后都一定能可靠得到一个 gene。

---

## 18. Gene annotation｜基因注释

**Gene annotation = 基因注释**

它告诉软件：

```text
哪个位置是什么 gene
哪些区域是 exon
哪些区域属于 transcript
gene 和 transcript 的对应关系
```

FASTA 和 annotation 的作用不同：

```text
FASTA
→ 这里是什么 sequence？

Annotation
→ 这里是什么 biological feature？
```

---

## 19. Gene × cell count matrix｜基因 × 细胞计数矩阵

**Gene × cell count matrix = 基因 × 细胞计数矩阵**

在经过：

```text
cell barcode
+
gene assignment
+
UMI deduplication
```

之后，可以得到：

```text
          Cell A   Cell B   Cell C
CD3D         3        0        2
MS4A1        0        4        0
LYZ          1        0        8
```

矩阵中的数字代表：

> 某个 cell 中某个 gene 检测到的 UMI-based molecule count。

之后这个矩阵才会进入 Scanpy。

例如后续可能存入：

```python
adata.X
```

---

## 20. Cell calling｜细胞识别

**Cell calling = 细胞识别 / 细胞判定**

不是每一个观察到的 barcode 都一定对应一个真实细胞。

因此软件需要判断：

```text
哪些 barcode
→ 很可能是真正的 cell

哪些 barcode
→ 更像 empty droplet / background
```

这个过程叫 cell calling。

---

## 21. Empty droplet｜空液滴

**Empty droplet = 空液滴**

在 droplet-based scRNA-seq 中，不是每个 droplet 都一定装有一个 cell。

有的 droplet 是空的。

但是空液滴里仍然可能出现少量 RNA。

---

## 22. Ambient RNA｜环境 RNA

**Ambient RNA = 环境 RNA / 游离 RNA**

样本处理中，如果有细胞破裂：

```text
cell lysis
↓
RNA leaked into solution
↓
ambient RNA
```

这些 RNA 可能进入 empty droplets。

于是：

```text
empty droplet
+
ambient RNA
↓
仍然产生少量 UMI / reads
```

因此：

> 有 barcode + RNA signal，并不一定就代表那里真的有一个完整细胞。

---

## 23. Barcode rank plot｜条形码排序图

**Barcode rank plot = 条形码排序图**

通常把 barcode 按 UMI 数量从高到低排列。

大致形状可能像：

```text
UMI count
↑
│\
│ \
│  \
│   \__
│      \________
└────────────────→ barcode rank
```

前面的高-count barcodes 通常更像真正的 cells。

后面低-count 的长尾中通常包含大量 empty droplets。

---

## 24. Raw matrix vs Filtered matrix｜原始矩阵与过滤矩阵

### Raw matrix

**Raw matrix = 原始计数矩阵**

通常包含比较广泛的 barcode，包括很多 background / empty droplet barcodes。

### Filtered matrix

**Filtered matrix = 过滤后的计数矩阵**

通常保留经过 cell calling 后，被认为更可能对应真实细胞的 barcodes。

下游分析常常从 filtered matrix 开始。

---

## 25. Salmon｜RNA-seq 定量 / mapping 工具

在我们当前流程中，它主要负责上游 mapping 等工作。

大致：

```text
FASTQ
↓
Salmon
↓
mapping
↓
产生供 alevin-fry 后续处理的信息
```

---

## 26. alevin-fry｜单细胞 barcode / UMI 定量工具

**alevin-fry = single-cell quantification framework，单细胞定量处理框架**

当前安装版本：

```text
alevin-fry 0.17.1
```

它主要参与：

```text
barcode processing
UMI processing
permit-list generation
collation
quantification
```

最终目标还是得到：

```text
gene × cell count matrix
```

---

## 27. splici reference｜spliced + intronic reference

**splici = spliced + intronic**

可以理解成：

**剪接转录本序列 + 内含子序列构成的 reference**

它包含：

```text
spliced transcript sequences
+
intronic sequences
```

为什么需要 intronic sequence？

因为真实 scRNA-seq reads 并不全部来自已经完全成熟、完全剪接的 mRNA。

有些 reads 可能来自：

```text
unspliced / partially processed RNA
```

所以加入 intronic sequence 可以帮助更合理地解释这些 reads。

---

# Core mental model｜核心思维模型

目前最重要的是记住这条身份链：

```text
                       Sequencing library
                               │
                 ┌─────────────┼─────────────┐
                 │             │             │
                I1            R1            R2
                 │             │             │
                 ▼        ┌────┴────┐        ▼
           Sample index   │         │      cDNA
                 │        ▼         ▼        │
                 │   Cell barcode   UMI      │
                 │        │         │        ▼
                 │        │         │     Mapping
                 │        │         │        │
                 │        │         │        ▼
                 │        │         │       Gene
                 │        │         │        │
                 └────────┴─────────┴────────┘
                               │
                               ▼
                    Cell + Gene + UMI
                               │
                       UMI deduplication
                               │
                               ▼
                    Gene × Cell matrix
```

