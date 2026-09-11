# 中文使用说明

这是为投稿/公开归档整理的最终代码仓库版本。

**最终 RWR restart probability 已锁定为 `alpha = 0.5`。** 旧版 `ALPHA = 0.7` 文件不放入公开仓库，避免审稿人误认为主分析存在参数冲突。

最重要的改动不是重写统计方法，而是：

- 删除旧版/无关代码；
- 将所有主分析路径改为相对路径；
- 将 STRING 400/700/900、CAP meta vector、core/full target lists 作为可直接运行的 derived inputs 整理；
- 新增一个统一的 downstream reproduction 入口；
- 新增自动数值核对，确保重新计算结果与最终 Supplementary Table 中的值一致；
- 明确区分“可以完全复现的 downstream network analysis”和“因 CNKI/WanFang/原始处方文件缺失而需要外部输入的 upstream GAT stage”。

完整运行：

```bash
python scripts/run_downstream_reproducibility.py
```

快速测试：

```bash
python scripts/run_downstream_reproducibility.py --smoke
```

如果准备上传 GitHub，建议先创建私有仓库测试一次；确认 README 页面、文件名和运行结果无误后再公开，并在 Zenodo 连接 GitHub 后生成 DOI。


## 发布前状态

本次整理已经实际运行并通过 downstream smoke test；自动核对脚本能够恢复网络规模、target 数量、CAP-up/CAP-down proximity、top-200 enrichment、positive-D Spearman correlation 与 PCI 等确定性主结果。完整 1,000 次 permutation robustness 计算量明显更大，建议在本地或服务器完成。

仓库已经加入 GitHub Actions smoke test。正式公开前请阅读 `docs/RELEASE_CHECKLIST.md`，并在所有共同作者确认后再确定最终 LICENSE、作者顺序及 Zenodo DOI。
