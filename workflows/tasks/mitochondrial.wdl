version 2.0
# 线粒体分析脚本
# 注意VEP的版本格式不符合标准，应在拉取时通过tag重订格式

import "struct.wdl"

## 线粒体Calling
### 直接从bam得到
### 线粒体是细胞质遗传，不遵循孟德尔定律，这里直接输出为vcf
### 参考 https://gatk.broadinstitute.org/hc/en-us/articles/4403870837275-Mitochondrial-short-variant-discovery-SNVs-Indels
task Mutect2Mito {
    input {
        String sample_id
        String output_dir
        File bam
        File bai
        IndexBundle reference
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
            -O ~{sample_id}.mt.raw.vcf.gz

        gatk FilterMutectCalls \
            -V ~{sample_id}.mt.raw.vcf.gz \
            -R ~{reference.fasta} \
            -O ~{sample_id}.mt.filtered.vcf.gz
        
        gatk SelectVariants \
            -V ~{sample_id}.mt.filtered.vcf.gz \
            --exclude-filtered true \
            -O ~{output_dir}/~{sample_id}.mt.vcf.gz
    >>>

    output {
        File vcf = "~{output_dir}/~{sample_id}.mt.vcf.gz"
        File tbi = "~{output_dir}/~{sample_id}.mt.vcf.gz.tbi"
    }

    runtime {
        container: "docker.io/broadinstitute/gatk:4.6.2.0"
        binding: "~{output_dir}:~{output_dir}"
        cpus: ~{threads}
    }
}

task MitoVEP {
    input {
        String sample_id
        File vcf
        Int threads
        File reference
        String output_dir
        String vep_database
    }

    command <<<
        if [ ! -d ~{output_dir} ]; then
            mkdir -p ~{output_dir}
        fi

        cache_str="Uploaded_variation,Location,REF_ALLELE,Allele,Consequence,IMPACT,DOMAINS,SYMBOL,HGNC_ID,HGVSg,Protein_position,Amino_acids,Codons,PUBMED,Existing_variation"
        custom_str="ClinVar_CLNSIG,ClinVar_CLNREVSTAT,ClinVar_CLNDN,ClinVar_CLNHGVS,MitoMap_aachange,MitoMap_homoplasmy,MitoMap_heteroplasmy,MitoMap_PubmedIDs,MitoMap_Disease,MitoMap_DiseaseStatus,MT_Gene"

        vep \
            --offline --cache \
            --dir_cache ~{vep_database} --merged \
            --force_overwrite --fork ~{threads} \
            -i ~{vcf} -o ~{output_dir}/~{sample_id}.mt.anno.vcf \
            --format vcf --vcf \
            --fa ~{reference} \
            --shift_3prime 1 --assembly GRCh37 --no_escape --check_existing -exclude_predicted --uploaded_allele --show_ref_allele --numbers --domains \
            --total_length --hgvs --hgvsg --symbol --ccds --uniprot --max_af --pubmed --pick \
            --custom file=~{vep_database}/clinvar/clinvar.vcf.gz,short_name=ClinVar,format=vcf,type=exact,coords=0,fields=CLNSIG%CLNREVSTAT%CLNDN%CLNHGVS \
            --custom file=~{vep_database}/MitoMap/mitomap.20230621.r8.disease.vcf.gz,short_name=MitoMap,format=vcf,type=exact,coords=0,fields=aachange%homoplasmy%heteroplasmy%PubmedIDs%Disease%DiseaseStatus \
            --custom file=~{vep_database}/MitoMap/NC_012920.1.bed.gz,short_name=MT_Gene,format=bed,type=overlap,coords=0 \
            --fields "${cache_str},${custom_str}"
    >>>

    output {
        File annoVcf = "~{output_dir}/~{sample_id}.mt.anno.vcf"
    }

    runtime {
        cpus: threads
        container: "docker.io/ensemblorg/ensembl-vep:release_114.2"
        binding: "~{output_dir}:~{output_dir},~{vep_database}:~{vep_database},~{plugin_dir}:~{plugin_dir}"
    }
}
