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



