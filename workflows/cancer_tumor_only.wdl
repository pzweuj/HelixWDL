version 2.0
# 肿瘤单样本分析流程

import "./tasks/struct.wdl"
import "./tasks/qc.wdl" as qc
import "./tasks/mapping.wdl" as mapping
import "./tasks/variant.wdl" as variant
import "./tasks/mitochondrial.wdl" as mito
import "./tasks/cnv.wdl" as cnv
import "./tasks/cancer_advance.wdl" as advance

workflow Cancer_Tumor_Only {
    input {
        String sample_id
        String output_dir
        File raw_read1
        File raw_read2
        File bed
        IndexBundle reference
        CnvBundle cnv_reference
        File msi_reference
        Int threads = 16
        Int flank = 20
        String sex
        Boolean report = false
        Boolean use_bwamem2 = true  # 默认使用BwaMem2，设为false使用Bwa
        String tmb_formula
    }





}

