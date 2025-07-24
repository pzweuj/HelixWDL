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
        BwaIndex reference
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
        BwaIndex reference
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
        BwaIndex reference
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





