version 2.0
# 通用比对脚本

import "struct.wdl"

## Bwa-mem2
### https://github.com/bwa-mem2/bwa-mem2
task BwaMem2 {
    input {
        String sample_id
        String output_dir
        File read1
        File read2
        Int threads
        BwaIndex reference
    }

    command <<<
        if [ ! -d ~{output_dir} ]; then
            mkdir -p ~{output_dir}
        fi

        bwa-mem2 mem -M -R "@RG\tID:mapping\tPL:helix\tLB:helix\tSM:~{sample_id}" \
            -t ~{threads} ~{reference.fasta} \
            ~{read1} ~{read2} | \
            samtools view -bSh - > ~{output_dir}/~{sample_id}.bam
        samtools sort -@ ~{threads} ~{output_dir}/~{sample_id}.bam -o ~{output_dir}/~{sample_id}.sorted.bam
        samtools index ~{output_dir}/~{sample_id}.sorted.bam
        rm ~{output_dir}/~{sample_id}.bam
    >>>

    output {
        File sort_bam = "~{output_dir}/~{sample_id}.sorted.bam"
        File sort_bam_bai = "~{output_dir}/~{sample_id}.sorted.bam.bai"
    }

    runtime {
        container: "ghcr.io/pzweuj/mapping:2025aug"
        binding: "~{output_dir}:~{output_dir}"
        cpu: ~{threads}
    }
}

## MarkDuplicates
### https://github.com/broadinstitute/gatk
task MarkDuplicates {
    




}


