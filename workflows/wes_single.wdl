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





}



