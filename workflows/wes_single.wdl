version 2.0
# WES分析流程

import "./utils/struct.wdl"
import "./utils/qc.wdl" as utils_qc
import "./utils/mapping.wdl" as utils_mapping

workflow WES_Single {
    input {
        String sample_id
        String output_dir
        File raw_read1
        File raw_read2
        Int threads
        File bed
        File cnv_baseline
        BwaIndex reference
    }

    # 质控
    String qc_output_dir = output_dir + "/QC"
    call utils_qc.Fastp as Fastp {input: sample_id=sample_id, output_dir=qc_output_dir, raw_read1=raw_read1, raw_read2=raw_read2, threads=threads}

    # 比对
    call utils_mapping.BwaMem2 as BwaMem2 {}




}



