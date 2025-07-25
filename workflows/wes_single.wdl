version 2.0
# WES分析流程

import "./tasks/struct.wdl"
import "./tasks/qc.wdl" as qc
import "./tasks/mapping.wdl" as mapping
import "./tasks/variant.wdl" as variant
import "./tasks/mitochondrial.wdl" as mito

workflow WES_Single {
    input {
        String sample_id
        String output_dir
        File raw_read1
        File raw_read2
        Int threads
        File bed
        File cnv_baseline
        IndexBundle reference
        Int flank
    }

    # 质控
    String qc_output_dir = output_dir + "/QC"
    call qc.Fastp as Fastp {input: sample_id=sample_id, output_dir=qc_output_dir, raw_read1=raw_read1, raw_read2=raw_read2, threads=threads}

    # 比对
    String bam_output_dir = output_dir + "/Bam"
    call mapping.BwaMem2 as BwaMem2 {input: sample_id=sample_id, output_dir=bam_output_dir, read1=Fastp.clean_read1, read2=Fastp.clean_read2, threads=threads, reference=reference}
    call mapping.MarkDuplicates as MarkDuplicates {input: sample_id=sample_id, output_dir=bam_output_dir, bam=BwaMem2.sort_bam, bai=BwaMem2.bai, threads=threads}

    # 质控2
    call qc.Bamdst as Bamdst {input: sample_id=sample_id, output_dir=qc_output_dir, bam=MarkDuplicates.mark_bam, bai=MarkDuplicates.mark_bai, bed=bed}
    call qc.Mosdepth as Mosdepth {input: sample_id=sample_id, output_dir=qc_output_dir, bam=MarkDuplicates.mark_bam, bai=MarkDuplicates.mark_bai, bed=bed, threads=4}
    call qc.CollectQCMetrics as CollectQCMetrics {input: sample_id=sample_id, output_dir=qc_output_dir, reference=reference, bam=MarkDuplicates.mark_bam, bai=MarkDuplicates.mark_bai, bed=bed}
    call qc.SRYCount as SRYCount {input: sample_id=sample_id, output_dir=qc_output_dir, fai=reference.fai, bam=MarkDuplicates.mark_bam, bai=MarkDuplicates.mark_bai}

    # SNP/InDel
    String vcf_output_dir = output_dir + "/Vcf"
    call variant.DeepVariant as DeepVariant {input: sample_id=sample_id, output_dir=vcf_output_dir, bam=MarkDuplicates.mark_bam, bai=MarkDuplicates.mark_bai, bed=bed, reference=reference, threads=threads, flank=flank}
    call variant.LeftAlignAndTrimVariants as LeftAlignAndTrimVariants {input: sample_id=sample_id, output_dir=vcf_output_dir, vcf=DeepVariant.vcf, idx=DeepVariant.tbi, reference=reference}
    call variant.WhatsHap as WhatsHap {input: sample_id=sample_id, output_dir=vcf_output_dir, vcf=LeftAlignAndTrimVariants.left_vcf, bam=MarkDuplicates.mark_bam, bai=MarkDuplicates.mark_bai, reference=reference}
    
    # CNV


    # 线粒体
    String mito_output_dir = output_dir + "/Mitochondrial"
    call mito.Mutect2Mito as Mutect2Mito {input: sample_id=sample_id, output_dir=mito_output_dir, bam=MarkDuplicates.mark_bam, bai=MarkDuplicates.mark_bai, reference=reference, threads=threads}

    # ROH


    # STR

    # 输出结果文件
    output {

    }


}



