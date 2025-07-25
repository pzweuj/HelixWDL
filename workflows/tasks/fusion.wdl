version 2.0
# fusion分析脚本

## Manta
### https://github.com/Illumina/manta
task MantaTumorOnly {
    input {
        String sample_id
        String output_dir
        File bam
        File bai
        File reference
        Int threads
        File bed
    }

    command <<<
        if [ ! -d ~{output_dir} ]; then
            mkdir -p ~{output_dir}
        fi

        bgzip ~{bed}
        tabix -p bed ~{bed}.gz

        configManta.py \
            --tumorBam ~{bam} \
            --referenceFasta ~{reference} \
            --exome \
            --generateEvidenceBam \
            --runDir . \
            --callRegions ~{bed}.gz
        
        ./runWorkflow.py -j ~{threads}

        zcat ./results/variants/tumorSV.vcf.gz > ~{output_dir}/~{sample_id}.manta.vcf
    >>>

    output {
        File vcf = "~{output_dir}/~{sample_id}.manta.vcf"
    }

    runtime {
        cpus: threads
        container: "dceoy/manta:latest"
        binding: "~{output_dir}:~{output_dir}"
    }
}

task MantaPair {
    input {
        String sample_id
        String normal_id
        String output_dir
        File bam
        File bai
        File normal_bam
        File normal_bai
        File reference
        Int threads
        File bed
    }

    command <<<
        if [ ! -d ~{output_dir} ]; then
            mkdir -p ~{output_dir}
        fi

        bgzip ~{bed}
        tabix -p bed ~{bed}.gz

        configManta.py \
            --tumorBam ~{bam} \
            --normalBam ~{normal_bam} \
            --referenceFasta ~{reference} \
            --exome \
            --generateEvidenceBam \
            --runDir . \
            --callRegions ~{bed}.gz
        
        ./runWorkflow.py -j ~{threads}

        zcat ./results/variants/somaticSV.vcf.gz > ~{output_dir}/~{sample_id}.manta.vcf
    >>>

    output {
        File vcf = "~{output_dir}/~{sample_id}.manta.vcf"
    }

    runtime {
        cpus: threads
        container: "dceoy/manta:latest"
        binding: "~{output_dir}:~{output_dir}"
    }
}

## TIDDIT
### https://github.com/SciLifeLab/TIDDIT
task TIDDIT {
    input {
        String sample_id
        String output_dir
        File bam
        File bai
        File reference
        Int threads
    }

    command <<<
        if [ ! -d ~{output_dir} ]; then
            mkdir -p ~{output_dir}
        fi

        python tiddit --sv \
            -o ~{output_dir}/~{sample_id}.tiddit \
            --threads ~{threads} \
            --bam ~{bam} --ref ~{reference} -z 15
    >>>

    output {
        File vcf = "~{output_dir}/~{sample_id}.tiddit.vcf"
    }

    runtime {
        cpus: threads
        container: "quay.io/biocontainers/tiddit:3.9.3--py310h20b60a1_0"
        binding: "~{output_dir}:~{output_dir}"
    }
}

## SVDB
### https://github.com/J35P312/SVDB
task SVDB {
    input {
        String sample_id
        String output_dir
        File vcf
        File db
    }

    command <<<
        if [ ! -d ~{output_dir} ]; then
            mkdir -p ~{output_dir}
        fi

        svdb --query \
            --query_vcf ~{vcf} \
            --out_occ GNOMAD_AC \
            --out_frq GNOMAD_AF \
            --in_occ AN \
            --out_frq AF \
            --db ~{db} \
            --prefix ~{output_dir}/~{sample_id}.svdb
    >>>

    output {
        File svdb_vcf = "~{output_dir}/~{sample_id}.svdb.vcf"
    }

    runtime {
        container: "clinicalgenomics/svdb:2.8.2"
        binding: "~{output_dir}:~{output_dir}"
    }
}


