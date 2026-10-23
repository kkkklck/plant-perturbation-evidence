# PlantVC v1 — GSE153103 Arabidopsis salt 6h 前瞻测试冻结规则

冻结时间：2026-10-20，在读取 6h 表达数值之前。

## 数据
GSE153103 experiment 1 raw count workbook。

样本：
- rep1: 0h RAS048, 1h RAS049, 3h RAS050, 6h RAS051
- rep2: 0h RAS056, 1h RAS057, 3h RAS058, 6h RAS059

冻结阶段只允许读取：RAS048, RAS049, RAS050, RAS056, RAS057, RAS058。

目标 RAS051 / RAS059 在 freeze 阶段禁止读取数值。

## 归一化
每个样本：CPM = count / library_size * 1e6。
表达值：log2(CPM + 1)。

对每个基因：
- r1 = mean(log2CPM_1h) - mean(log2CPM_0h)
- r3 = mean(log2CPM_3h) - mean(log2CPM_0h)

## PlantVC v1 持续型历史分数

若 sign(r1) == sign(r3) 且二者非零：

`persistent_bonus = sqrt(abs(r1 * r3))`

否则：

`persistent_bonus = 0`

最终预测强度：

`magnitude = 0.5 * abs(r3) + 0.5 * persistent_bonus`

预测方向：

`direction = sign(r3)`

有符号预测：

`pred = direction * magnitude`

不做 3h→6h 数值线性外推，避免两个早期时间点造成过度外推。

## 简单基线
1. 3h baseline：pred = r3
2. 1h baseline：pred = r1

PlantVC v1 的价值必须和这些简单基线比较，而不是只和随机比较。

## 基因范围
主分析技术范围：在六个 pretarget 样本中，CPM >= 1 的样本数至少为 2。

同时保留 all-gene sensitivity，避免技术过滤成为人为生杀门槛。

## 置信度（只报告，不参与主排名）
对每个基因计算：
- 1h 两个 matched replicate response 是否同方向；
- 3h 两个 matched replicate response 是否同方向；
- 平均 r1 与 r3 是否同方向。

confidence = 三项中一致项比例，范围 0–1。

## 目标揭晓后的评价
目标：6h 相对 0h 的平均 log2(CPM+1) response r6。

连续报告：
- full-vector Spearman；
- top100 / 250 / 500 / 1000 absolute-response overlap 与随机富集倍数；
- overlap 中方向一致率；
- target top-k 在预测排名中的百分位；
- PlantVC v1 vs 3h baseline vs 1h baseline；
- 两个 6h biological replicates 自己之间的 top-k 重复性，作为实验上限参考。

不设项目级 PASS/FAIL 阈值。

## 解释边界
这是持续盐胁迫、bulk RNA-seq 的 v1 单专家测试，不等于完整单细胞路由系统验证。
它测试的是：持续型历史证据是否能提高未来强响应候选的排序。
