version 2.0
# Variant Call

import "struct.wdl"

## deepVariant
### https://github.com/google/deepvariant
task DeepVariant {
    input {
        String sample_id
        String output_dir
        File bam
        File bai
        File bed
        IndexBundle reference
        Int threads
        Int flank
    }

    command <<<
        if [ ! -d ~{output_dir} ]; then
            mkdir -p ~{output_dir}
        fi

        awk 'BEGIN {OFS="\t"} {start=$2; end=$3; new_start=start-~{flank}; new_end=end+~{flank}; if(new_start<0) new_start=0; print $1,new_start,new_end}' ~{bed} > extended.bed

        /opt/deepvariant/bin/run_deepvariant \
            --model_type WES \
            --ref ~{reference.fasta} \
            --reads ~{bam} \
            --output_vcf ~{output_dir}/~{sample_id}.vcf.gz \
            --output_gvcf ~{output_dir}/~{sample_id}.g.vcf.gz \
            --num_shards ~{threads} \
            --regions extended.bed
    >>>

    output {
        File vcf = "~{output_dir}/~{sample_id}.vcf.gz"
        File tbi = "~{output_dir}/~{sample_id}.vcf.gz.tbi"
        File g_vcf = "~{output_dir}/~{sample_id}.g.vcf.gz"
        File g_tbi = "~{output_dir}/~{sample_id}.g.vcf.gz.tbi"
    }

    runtime {
        cpus: threads
        container: "google/deepvariant:1.9.0"
        binding: "~{output_dir}:~{output_dir}"
    }
}

## deepSomatic
### https://github.com/google/deepsomatic
task DeepSomaticTumorOnly {
    input {
        String sample_id
        String output_dir
        File bam
        File bai
        File bed
        IndexBundle reference
        Int threads
    }

    command <<<
        if [ ! -d ~{output_dir} ]; then
            mkdir -p ~{output_dir}
        fi

        run_deepsomatic \
            --model_type=FFPE_WES_TUMOR_ONLY \
            --ref=~{reference.fasta} \
            --reads_tumor=~{bam} \
            --output_vcf=~{output_dir}/~{sample_id}.ds.vcf.gz \
            --sample_name_tumor=~{sample_id} \
            --num_shards=~{threads} \
            --logging_dir=logs \
            --intermediate_results_dir=intermediate_results_dir \
            --use_default_pon_filtering=true \
            --regions ~{bed}
    >>>

    output {
        File vcf = "~{output_dir}/~{sample_id}.ds.vcf.gz"
        File tbi = "~{output_dir}/~{sample_id}.ds.vcf.gz.tbi"
    }

    runtime {
        cpus: threads
        container: "google/deepsomatic:1.9.0"
        binding: "~{output_dir}:~{output_dir}"
    }
}

task DeepSomaticPair {
    input {
        String sample_id
        String normal_id
        String output_dir
        File bam
        File bai
        File normal_bam
        File normal_bai
        File bed
        IndexBundle reference
        Int threads
    }

    command <<<
        if [ ! -d ~{output_dir} ]; then
            mkdir -p ~{output_dir}
        fi

        run_deepsomatic \
            --model_type=FFPE_WES \
            --ref=~{reference.fasta} \
            --reads_tumor=~{bam} \
            --reads_normal=~{normal_bam} \
            --output_vcf=~{output_dir}/~{sample_id}.ds.vcf.gz \
            --sample_name_tumor=~{sample_id} \
            --sample_name_normal=~{normal_id} \
            --num_shards=~{threads} \
            --logging_dir=logs \
            --intermediate_results_dir=intermediate_results_dir \
            --use_default_pon_filtering=true \
            --regions ~{bed}
    >>>

    output {
        File vcf = "~{output_dir}/~{sample_id}.ds.vcf.gz"
        File tbi = "~{output_dir}/~{sample_id}.ds.vcf.gz.tbi"
    }

    runtime {
        cpus: threads
        container: "google/deepsomatic:1.9.0"
        binding: "~{output_dir}:~{output_dir}"
    }
}

## Mutect2
### https://github.com/broadinstitute/gatk
task Mutect2TumorOnly {

}


task Mutect2Pair {

}

## 左对齐
### 用于将位点按照基因组坐标进行左对齐，并分拆多突变方向
task LeftAlignAndTrimVariants {
    input {
        String sample_id
        String output_dir
        File vcf
        File idx
        IndexBundle reference
    }

    command <<<
        if [ ! -d ~{output_dir} ]; then
            mkdir -p ~{output_dir}
        fi

        gatk LeftAlignAndTrimVariants \
            -R ~{reference.fasta} \
            -V ~{vcf} \
            -O ~{output_dir}/~{sample_id}.left.vcf.gz \
            --split-multi-allelics \
            --create-output-variant-index true
    >>>

    output {
        File left_vcf = "~{output_dir}/~{sample_id}.left.vcf.gz"
        File left_vcf_tbi = "~{output_dir}/~{sample_id}.left.vcf.gz.tbi"
    }

    runtime {
        container: "broadinstitute/gatk:4.6.2.0"
        binding: "~{output_dir}:~{output_dir}"
    }
}

## 单倍型检查
## WhatsHap
### https://whatshap.readthedocs.io/en/latest/
task WhatsHap {
    input {
        String sample_id
        String output_dir
        File vcf
        File bam
        File bai
        IndexBundle reference
    }

    command <<<
        if [ ! -d ~{output_dir} ]; then
            mkdir -p ~{output_dir}
        fi

        whatshap phase \
            --indels \
            --reference=~{reference.fasta} \
            -o ~{output_dir}/~{sample_id}.phase.vcf \
            ~{vcf} \
            ~{bam}
        
        bgzip ~{output_dir}/~{sample_id}.phase.vcf
        tabix -f -p vcf ~{output_dir}/~{sample_id}.phase.vcf.gz
    >>>

    output {
        File phase_vcf = "~{output_dir}/~{sample_id}.phase.vcf.gz"
        File phase_vcf_tbi = "~{output_dir}/~{sample_id}.phase.vcf.gz.tbi"
    }

    runtime {
        container: "fellen31/whatshap-tabix:2.2"
        singularity_binding: "~{output_dir}:~{output_dir}"
    }

}


