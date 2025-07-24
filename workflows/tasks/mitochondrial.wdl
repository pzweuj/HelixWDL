version 2.0
# 线粒体分析脚本

import "struct.wdl"

## 线粒体Calling
### 直接从bam得到
### 线粒体是细胞质遗传，不遵循孟德尔定律，这里直接输出为vcf
### 参考 https://gatk.broadinstitute.org/hc/en-us/articles/4403870837275-Mitochondrial-short-variant-discovery-SNVs-Indels
task MitoCalling {
    input {
        String sample_id
        String output_dir
        File bam
        File bai
        GATKIndex reference
        Int threads
    }

    command <<<
        if [ ! -d ~{output_dir} ]; then
            mkdir -p ~{output_dir}
        fi

        # 同时指定MT和chrM，GATK会自动选择存在的染色体
        gatk Mutect2 \
            -R ~{reference.fasta} \
            -L MT -L chrM \
            --mitochondria-mode \
            -I ~{bam} \
            -O ~{output_dir}/~{sample_id}.mt.vcf.gz

        gatk FilterMutectCalls \
            -V ~{output_dir}/~{sample_id}.mt.vcf.gz \
            -R ~{reference.fasta} \
            -O ~{output_dir}/~{sample_id}.mt.filtered.vcf.gz
        
        gatk SelectVariants \
            -V ~{output_dir}/~{sample_id}.mt.filtered.vcf.gz \
            --exclude-filtered true \
            -O ~{output_dir}/~{sample_id}.mt.vcf
    >>>

    output {
        File vcf = "~{output_dir}/Mitochondrion/~{sample_id}.mt.vcf"
        File idx = "~{output_dir}/Mitochondrion/~{sample_id}.mt.vcf.idx"
    }

    runtime {
        container: "broadinstitute/gatk:4.6.2.0"
        binding: "~{output_dir}:~{output_dir}"
    }
}
