version 2.0
# 通用质控脚本

## Fastp
### https://github.com/OpenGene/fastp
task Fastp {
    input {
        String sample_id
        String output_dir
        File raw_read1
        File raw_read2
        Int threads
    }

    # fastp只支持最多16线程
    Int fastp_threads = if threads > 16 then 16 else threads

    command <<<
        if [ ! -d ~{output_dir} ]; then
            mkdir -p ~{output_dir}
        fi

        fastp \
            -i ~{raw_read1} \
            -I ~{raw_read2} \
            -o ~{output_dir}/~{sample_id}.clean_1.fq.gz \
            -O ~{output_dir}/~{sample_id}.clean_2.fq.gz \
            -w ~{fastp_threads} \
            -j ~{output_dir}/~{sample_id}.fastp_stats.json \
            -h ~{output_dir}/~{sample_id}.fastp_stats.html \
            --detect_adapter_for_pe
    >>>

    output {
        File clean_read1 = "~{output_dir}/~{sample_id}.clean_1.fq.gz"
        File clean_read2 = "~{output_dir}/~{sample_id}.clean_2.fq.gz"
        File json_report = "~{output_dir}/~{sample_id}.fastp_stats.json"
        File html_report = "~{output_dir}/~{sample_id}.fastp_stats.html"
    }

    runtime {
        container: "ghcr.io/pzweuj/mapping:2025aug"
        binding: "~{output_dir}:~{output_dir}"
        cpus: ~{fastp_threads}
    }
}

## Bamdst
### https://github.com/shiquan/bamdst
task Bamdst {
    input {
        String sample_id
        String output_dir
        File bam
        File bai
        File bed
    }

    command <<<
        if [ ! -d ~{output_dir} ]; then
            mkdir -p ~{output_dir}
        fi

        bamdst -p ~{target_bed} \
            --cutoffdepth 20 \
            -o ~{output_dir} ~{bam}
    >>>

    output {
        File bamdst_cov_file = "~{output_dir}/coverage.report"
        File bamdst_region_file = "~{output_dir}/region.tsv.gz"
    }

    runtime {
        container: "ghcr.io/pzweuj/mapping:2025aug"
        binding: "~{output_dir}:~{output_dir}"
    }

}

## Mosdepth
### https://github.com/brentp/mosdepth
### 一般是WGS才用
task Mosdepth {
    input {
        String sample_id
        String output_dir
        File bam
        File bai
        File? bed
        Int threads = 4
    }

    command <<<
        if [ ! -d ~{output_dir} ]; then
            mkdir -p ~{output_dir}
        fi

        # 基本的深度统计
        mosdepth_d4 \
            --threads ~{threads} \
            ~{if defined(bed) then "--by " + bed else ""} \
            --no-per-base \
            ~{output_dir}/~{sample_id} \
            ~{bam}

        # 生成全基因组深度统计（如果没有提供bed文件）
        ~{if !defined(bed) then "mosdepth --threads " + threads + " --by 500 " + output_dir + "/" + sample_id + ".windows " + bam else ""}
    >>>

    output {
        File mosdepth_summary = "~{output_dir}/~{sample_id}.mosdepth.summary.txt"
        File? mosdepth_regions = if defined(bed) then "~{output_dir}/~{sample_id}.regions.bed.gz" else None
        File? mosdepth_windows = if !defined(bed) then "~{output_dir}/~{sample_id}.windows.regions.bed.gz" else None
        File mosdepth_dist = "~{output_dir}/~{sample_id}.mosdepth.global.dist.txt"
    }

    runtime {
        container: "ghcr.io/pzweuj/mapping:2025aug"
        binding: "~{output_dir}:~{output_dir}"
        cpus: ~{threads}
    }
}

## CollectQCMetrics
task CollectQCMetrics {
    input {
        String sample_id
        String output_dir
        File bam
        File bai
        File bed
        IndexBundle reference
    }

    command <<<
        if [ ! -d ~{output_dir} ]; then
            mkdir -p ~{output_dir}
        fi

        gatk CollectInsertSizeMetrics \
            -I ~{bam} \
            -O ~{output_dir}/~{sample_id}.insertsize.txt \
            -H ~{output_dir}/~{sample_id}.histogram.pdf

        gatk CollectAlignmentSummaryMetrics \
            -I ~{bam} \
            -R ~{reference.fasta} \
            -O ~{output_dir}/~{sample_id}.align_summary.txt

        gatk BedToIntervalList \
            -I ~{bed} \
            -O ~{sample_id}.interval_list \
            -SD ~{reference.fasta}

        gatk CollectHsMetrics \
            -BI ~{sample_id}.interval_list \
            -TI ~{sample_id}.interval_list \
            -I ~{bam} \
            -O ~{output_dir}/~{sample_id}.hs_metric.txt
    >>>

    output {
        File insertsize_file = "~{output_dir}/~{sample_id}.insertsize.txt"
        File insertsize_pdf = "~{output_dir}/~{sample_id}.histogram.pdf"
        File align_summary = "~{output_dir}/~{sample_id}.align_summary.txt"
        File hs_metric = "~{output_dir}/~{sample_id}.hs_metric.txt"
    }

    runtime {
        container: "broadinstitute/gatk:4.6.2.0"
        binding: "~{output_dir}:~{output_dir}"
    }
}

### SRY计数
task SRYCount {
    input {
        String sample_id
        String output_dir
        File fai
        File bam
        File bai
    }

    command <<<
        if [ ! -d ~{output_dir} ]; then
            mkdir -p ~{output_dir}
        fi

        # 检测Y染色体信息
        y_info=$(awk '$1 ~ /^(chr)?Y$/ {print $1, $2}' ~{fai})
        
        if [ -z "$y_info" ]; then
            echo "Error: Y chromosome not found in reference" >&2
            exit 1
        fi
        
        y_chr=$(echo $y_info | cut -d' ' -f1)
        y_length=$(echo $y_info | cut -d' ' -f2)
        
        # 根据Y染色体长度判断基因组版本并设置SRY坐标
        case $y_length in
            57227415)
                # hg38/GRCh38
                sry_start=2786855
                sry_end=2787682
                genome_version="hg38"
                ;;
            59373566)
                # hg19/GRCh37
                sry_start=2654896
                sry_end=2655723
                genome_version="hg19"
                ;;
            *)
                echo "Warning: Unknown genome version (Y length: $y_length), using hg19 coordinates" >&2
                sry_start=2654896
                sry_end=2655723
                genome_version="unknown"
                ;;
        esac
        
        sry_region="${y_chr}:${sry_start}-${sry_end}"
        
        echo "Detected: $genome_version, Y chromosome: $y_chr, SRY region: $sry_region" >&2
        
        # 统计SRY区域的reads数量
        sry_count=$(samtools view -F 2052 ~{bam} ${sry_region} | wc -l)
        
        # 输出结果，包含更多信息
        echo -e "Sample\tGenome_Version\tY_Chromosome\tSRY_Region\tSRY_Count" > ~{output_dir}/~{sample_id}.SRY.txt
        echo -e "~{sample_id}\t${genome_version}\t${y_chr}\t${sry_region}\t${sry_count}" >> ~{output_dir}/~{sample_id}.SRY.txt
    >>>

    output {
        File sry_result = "~{output_dir}/~{sample_id}.SRY.txt"
    }

    runtime {
        container: "ghcr.io/pzweuj/mapping:2025aug"
        binding: "~{output_dir}:~{output_dir}"
    }
}


