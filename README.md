<div align="center">
  <img src="helix_logo.svg" alt="HelixWDL Logo" width="200">
  <h1>HelixWDL</h1>
  <p><strong>基于WDL语言的NGS分析流程，专注于临床肿瘤分析和遗传WES分析</strong></p>
</div>

## ⚠️ 注意 ⚠️

项目处于快速迭代阶段，部分功能尚未完善，技术文档尚未更新，不建议用于临床分析。

## 项目概述

HelixWDL是一个完全基于WDL构建的下一代测序(NGS)数据分析流程。项目采用容器化技术，使用公开数据库，确保所有用户都能够轻松复现和快速部署。

**设计理念**: 本项目专注于从原始测序数据(Fastq)到变异位点识别(VCF)的核心分析流程，提供高质量的变异检测和基础质量过滤。我们有意将变异注释作为可选功能，这样的模块化设计使得分析结果可以灵活对接各种下游分析系统，用户可以根据具体的临床或研究需求选择最适合的注释和解读工具。

### 核心特性

- 🧬 **全面的NGS分析**: 支持WES、肿瘤配对/单样本、CNV和MSI分析
- 🐳 **Docker容器化**: 所有工具和依赖都通过Docker封装，确保环境一致性
- 📊 **公开数据库**: 使用公开可获取的参考数据库
- 🔄 **可重现性**: 标准化流程确保结果的可重现性
- ⚡ **快速部署**: 简化的配置和部署流程
- 🎯 **临床导向**: 专为临床应用设计的分析流程

## 分析流程

### 1. WES分析流程
- **单样本WES分析** (`wes_single.wdl`): 从Fastq到VCF的完整单样本外显子组分析
- **家系WES分析** (`wes_trio.wdl`): 从gVCF到VCF的多样本家系的遗传变异分析

### 2. 肿瘤分析流程
- **肿瘤单样本分析** (`cancer_single.wdl`): 从Fastq到VCF的肿瘤组织单独分析
- **肿瘤配对分析** (`cancer_pair.wdl`): 从Fastq到VCF的肿瘤-正常组织配对分析

## 技术架构

### 工作流引擎
- **Cromwell**: 支持本地、集群(LSF/Slurm)等多种执行环境
- **WDL 2.0**: 使用最新的WDL语言规范

### 容器化部署
- **Docker**: 所有生物信息学工具都封装在Docker镜像中
- **Ubuntu 24.04**: 基于最新LTS版本构建
- **Singularity**: 同时支持Singularity
- **Podman**: 同时支持Podman

### 执行环境支持
- **单机模式**: 适合小规模数据处理
- **LSF集群**: 支持IBM LSF作业调度系统
- **Slurm集群**: 支持Slurm作业调度系统

## 目录结构

```
HelixWDL/
├── workflows/                 # WDL工作流文件
│   └── tasks/                 # 任务模块
├── containers/                # Docker容器配置
├── cromwell_config/           # Cromwell配置文件
├── scripts/                   # 辅助脚本
├── inputs/                    # 输入文件目录
├── README.md                  # 项目说明文档
└── LICENSE                    # 开源许可证
```

## 使用说明

详细的部署和使用说明请参考项目[Wiki](https://github.com/pzweuj/HelixWDL/wiki)文档。

## 许可证

本项目采用MIT许可证，详见[LICENSE](LICENSE)文件。
