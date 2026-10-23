# PlantVC v1 — GSE158761 20h 外部盲测冻结规则（2026-10-20）

## 目的
这是统一 PlantVC v1 入口的外部持续/进展型测试。

目标不是重新为该数据集发明规则，而是直接使用已经前瞻验证过的持续型规则：

> 最近一次可靠的同状态响应，作为后续阶段候选排序的默认证据。

## 数据与时间
作者标注对象：GSE158761_all_cells_labeled.rds

metadata 映射：
- 0h：sample 1、5
- 8h：sample 2、6
- 20h：sample 3、7

生物重复按编号轨迹配对：
- replicate A：1 -> 2 -> 3
- replicate B：5 -> 6 -> 7

## 目标隔离
冻结预测阶段：
- 允许读取并计算 0h、8h 表达；
- 20h counts 不参与任何 pseudobulk、技术基因筛选、候选排序、置信描述或参数选择；
- 20h 只在预测排名和哈希保存完成后才允许用于评分。

RDS 文件整体会被反序列化并包含全部细胞，但冻结代码在任何表达统计前通过 metadata 索引选择 0h/8h 细胞；20h 表达值不用于预测。

## 主纳入规则（仅由 metadata 决定）
一个作者 `celltype` 只有在 6 个样本（1,2,3,5,6,7）中每个样本都至少有 20 个细胞，才进入主盲测。

冻结得到的 9 个 eligible cell types：
1. Ambiguous Stele Cells
2. Columella/ Root Cap
3. Epidermis
4. Lateral Root Endodermis
5. Mature Pericycle
6. Phloem
7. Phloem Pole Pericycle
8. Xylem
9. Xylem Pole Pericycle

以下状态在看表达前即排除，不进入主结果：
- Cortex
- Endodermis
- Lateral Root Primordia

## 表达汇总与技术基因范围
对每个 eligible cell type 独立计算 pseudobulk。

每个样本内：
- 对该 cell type 的细胞 counts 求和；
- 转成 CPM；
- 使用 `log2(CPM + 1)`。

冻结阶段技术基因范围（每个 cell type 独立）：
- 在 4 个 pretarget pseudobulk（sample 1,2,5,6）中，至少 2 个样本 CPM >= 1。

技术基因筛选完全不使用 20h。

## PlantVC v1 持续型预测
对每个 eligible cell type：

- replicate A 的 8h 响应：`r8_A = logCPM(sample2) - logCPM(sample1)`
- replicate B 的 8h 响应：`r8_B = logCPM(sample6) - logCPM(sample5)`
- 最近历史：`r8 = (r8_A + r8_B) / 2`

统一 PlantVC v1 的主预测：
- `predicted_response = r8`
- `priority_score = abs(r8)`
- 方向 = sign(r8)

`r8_A` / `r8_B` 只用于描述重复一致性，不改变候选排名。

## 冻结对照：跨细胞共同响应
为了直接检验“类型优先”的证据选择是否合理，对每个目标 cell type 同时冻结一个不参与主模型调参的对照：

- 使用另外 8 个 eligible cell types 的 `r8`；
- 对每个基因做等权有符号平均；
- 得到 `broad_shared_r8`；
- 按 `abs(broad_shared_r8)` 排名。

这个对照对应急性型 PlantVC 使用的“跨状态共同响应”思路。

因此本次测试会同时回答：
1. 最近同状态历史本身能否预测 20h；
2. 在这个持续/进展型问题里，它是否比跨状态 shared predictor 更合适。

该对照完全由 0h/8h 构建，不读取 20h。

## 20h 目标定义（揭晓后）
评分阶段才允许计算：
- `r20_A = logCPM(sample3) - logCPM(sample1)`
- `r20_B = logCPM(sample7) - logCPM(sample5)`
- `r20 = (r20_A + r20_B) / 2`

## 评价
每个 cell type 独立报告：
- 完整向量 Spearman rho（辅助指标）；
- 全基因方向一致率（辅助指标）；
- top100 / 250 / 500 / 1000 候选重叠；
- 相对随机的富集倍数；
- top-k 重叠基因方向一致率；
- 真实 top-k 在预测排名中的中位百分位；
- 20h 两个生物重复之间的可重复性，作为目标可靠性参考。

主解释关注候选排序，不要求完整未来转录组精确重建。

## 禁止事项
- 不因 20h 结果修改 eligible cell type；
- 不因 20h 结果修改技术基因阈值；
- 不因 20h 结果修改时间配对；
- 不因 20h 结果给 r8 增加额外权重、平滑或模型；
- 不因 20h 结果修改或删除 frozen broad-shared 对照；
- 若结果较弱，记录适用边界，不事后发明新公式挽救。
