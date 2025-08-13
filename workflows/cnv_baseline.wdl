version 2.0
# CNV基线建立流程
# 不使用refflat进行注释
# https://github.com/pzweuj/ManeSelectBed
# 定义为外显子级别的CNV分析

import "./tasks/cnv.wdl" as cnv

workflow CNV_Baseline {
    input {
        String prefix
        String output_dir
        File bed
        File mane_bed
        File exclude_bed
        File reference
        String reference_bam_dir
        Boolean use_chr = true
    }

    call cnv.ManeBedIntersect as ManeBedIntersect {input: prefix=prefix, output_dir=output_dir, bed=bed, mane_bed=mane_bed, use_chr=use_chr}
    call cnv.CNVkitTargetAndAntitarget as CNVkitTargetAndAntitarget {input: prefix=prefix, output_dir=output_dir, bed=ManeBedIntersect.cnv_bed, exclude_bed=exclude_bed, reference=reference}

    String coverage_dir = output_dir + "/Coverage"
    
    # 获取reference_bam_dir下的所有BAM文件
    Array[File] bam_files = glob(reference_bam_dir + "/*.bam")
    
    # 对每个BAM文件执行CNVkitCoverage任务
    scatter (bam_file in bam_files) {
        # 从BAM文件路径提取样本ID（去掉路径和.bam扩展名）
        String sample_id = basename(bam_file, ".bam")
        
        # 构造对应的BAI文件路径
        File bai_file = bam_file + ".bai"
        
        call cnv.CNVkitCoverage as CNVkitCoverage {
            input: 
                sample_id = sample_id,
                output_dir = coverage_dir,
                bam = bam_file,
                bai = bai_file,
                target_bed = CNVkitTargetAndAntitarget.target_bed,
                antitarget_bed = CNVkitTargetAndAntitarget.antitarget_bed
        }
    }
    
    # 建立CNV基线
    call cnv.CNVkitReference as CNVkitReference {
        input:
            prefix = prefix,
            coverage_dir = coverage_dir,
            output_dir = output_dir,
            reference = reference
    }
    
    output {
        File target_bed = CNVkitTargetAndAntitarget.target_bed
        File antitarget_bed = CNVkitTargetAndAntitarget.antitarget_bed
        File reference_cnn = CNVkitReference.ref_cnn
    }
}
