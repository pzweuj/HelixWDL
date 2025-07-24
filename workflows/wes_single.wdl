version 2.0
# WES分析流程

import "./tasks/struct.wdl"
import "./tasks/qc.wdl" as qc
import "./tasks/mapping.wdl" as mapping

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
    call qc.Fastp as Fastp {input: sample_id=sample_id, output_dir=qc_output_dir, raw_read1=raw_read1, raw_read2=raw_read2, threads=threads}

    # 比对
    String bam_output_dir = output_dir + "/Bam"
    call mapping.BwaMem2 as BwaMem2 {input: sample_id=sample_id, output_dir=bam_output_dir, read1=Fastp.clean_read1, read2=Fastp.clean_read2, threads=threads, reference=reference}
    call mapping.MarkDuplicates as MarkDuplicates {input: sample_id=sample_id, output_dir=bam_output_dir, bam=BwaMem2.sort_bam, bai=BwaMem2.bai, threads=threads}

    # 质量报告


    # SNP/InDel


    # CNV


    # 线粒体


    # 




}



