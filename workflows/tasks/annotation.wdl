version 2.0
# 注释流程

# Exomiser
# --spring.config.location必须是最后一个参数
task Exomiser {
    input {
        String sample_id
        String output_dir
        String hpo_list
        File vcf
        File config_yaml
        File reference
        String genome
        String data_version
        String db_path
    }

    command <<<
        if [ ! -d ~{output_dir} ]; then
            mkdir -p ~{output_dir}
        fi

        hpo_list_string=`awk '{printf "%s,", $0}' ~{hpo_list} | sed 's/,$//'`

        /deploy/run_exomiser.sh \
            --sample_id ~{sample_id} \
            --vcf_file ~{vcf} \
            --output_dir ~{output_dir} \
            --genome ~{genome} \
            --hpo_string $hpo_list_string \
            --data_dir /data/exomiser \
            --data_version ~{data_version}
    >>>

    output {
        File result_file = "~{output_dir}/~{sample_id}.variants.tsv"
    }

    runtime {
        container: "ghcr.io/pzweuj/exomiser:v14.1.0"
        binding: "~{output_dir}:~{output_dir},~{db_path}:/data/exomiser"
        cpus: 4
    }
}

# VEP
## 需要确保所需的字段都有
## --shift_hgvs 1 用于在转录本水平启动3'原则，默认的开启的
## --shift_3prime 1 用于在基因组水平启动3'原则，默认是关闭的
## 注，VEP插件可注释spiceAI、dbNSFP、dbscsnv 需要自行从相关网站下载数据库
## 注释主要需求：1，rsID；2，Clinvar；3，有害性预测；4，人群频率；5，cytoBand；6，基因；7，HGVS
## 指定--pick_allele来对每个方向的变异选择一个注释结果，优先顺序查看
## https://asia.ensembl.org/info/docs/tools/vep/script/vep_other.html#pick
## exac数据库已弃用，已经并入GnomAD项目了
## cytoband下载后转成bed格式，按以下连接操作
## https://asia.ensembl.org/info/docs/tools/vep/script/vep_custom.html
## 输出为vcf，便于整理深度等信息
## 注释时即过滤EMSEMBL的转录本 --transcript_filter "stable_id match N[MR]_"
task VEP {
    input {
        String sample_id
        File vcf
        Int threads
        File reference
        String output_dir
        String vep_database
        String plugin_dir
    }

    command <<<
        if [ ! -d ~{output_dir} ]; then
            mkdir -p ~{output_dir}
        fi
 
        cache_str="Uploaded_variation,Location,REF_ALLELE,Allele,Consequence,IMPACT,DOMAINS,Feature,DISTANCE,EXON,INTRON,SYMBOL,STRAND,HGNC_ID,HGVSc,HGVSp,HGVSg,MAX_AF,Protein_position,Amino_acids,Codons,PUBMED,Existing_variation"
        custom_str="cytoBand,ClinVar_CLNSIG,ClinVar_CLNREVSTAT,ClinVar_CLNDN,ClinVar_CLNHGVS,clinPath,InterVar_Intervar,GnomHemi_Gnomad_AC_XY"
        dbnsfp_str="rs_dbSNP,REVEL_score,REVEL_rankscore,M-CAP_score,M-CAP_rankscore,M-CAP_pred,GERP++_RS,GERP++_RS_rankscore,MVP_score,MVP_rankscore,phyloP100way_vertebrate,phyloP100way_vertebrate_rankscore,phyloP470way_mammalian,phyloP470way_mammalian_rankscore"
        dbnsfp_gnomad_str="gnomAD_exomes_AF,gnomAD_exomes_AN,gnomAD_exomes_AC,gnomAD_exomes_EAS_AF,gnomAD_exomes_EAS_AN,gnomAD_exomes_EAS_AC,gnomAD_exomes_POPMAX_AF,gnomAD_exomes_POPMAX_AN,gnomAD_exomes_POPMAX_AC,gnomAD_exomes_nhomalt,gnomAD_genomes_AF,gnomAD_genomes_AN,gnomAD_genomes_AC,gnomAD_genomes_EAS_AF,gnomAD_genomes_EAS_AN,gnomAD_genomes_EAS_AC,gnomAD_genomes_POPMAX_AF,gnomAD_genomes_POPMAX_AN,gnomAD_genomes_POPMAX_AC,gnomAD_genomes_nhomalt"
        spliceai_str="SpliceAI_pred_DP_AG,SpliceAI_pred_DP_AL,SpliceAI_pred_DP_DG,SpliceAI_pred_DP_DL,SpliceAI_pred_DS_AG,SpliceAI_pred_DS_AL,SpliceAI_pred_DS_DG,SpliceAI_pred_DS_DL"
        self_plugin_str="FlankingSequence,MissenseZscore"

        vep \
            --offline --cache \
            --dir_cache ~{vep_database} --merged \
            --dir_plugins ~{plugin_dir} \
            --force_overwrite --fork ~{threads} \
            -i ~{vcf} -o ~{output_dir}/~{sample_id}.vep.anno.vcf \
            --format vcf --vcf \
            --fa ~{reference} \
            --shift_3prime 1 --assembly GRCh37 --no_escape --check_existing -exclude_predicted --uploaded_allele --show_ref_allele --numbers --domains \
            --total_length --hgvs --hgvsg --symbol --ccds --uniprot --max_af --pubmed \
            --transcript_filter "stable_id match N[MR]_" \
            --custom file=~{vep_database}/clinvar/clinvar.vcf.gz,short_name=ClinVar,format=vcf,type=exact,coords=0,fields=CLNSIG%CLNREVSTAT%CLNDN%CLNHGVS \
            --custom file=~{vep_database}/intervar/intervar.vcf.gz,short_name=InterVar,format=vcf,type=exact,coords=0,fields=Intervar \
            --custom file=~{vep_database}/cytoband/cytoBand.bed.gz,short_name=cytoBand,format=bed,type=overlap,coords=0 \
            --custom file=~{vep_database}/clinvar/clinvar_pathogenic.merge.bed.gz,short_name=clinPath,format=bed,type=overlap,coords=0 \
            --custom file=~{vep_database}/GnomAD/gnomad.v4.1.sites.combined_AC_XY.b37.vcf.gz,short_name=GnomHemi,format=vcf,type=exact,coords=0,fields=Exomes_AC_XY%Genomes_AC_XY%Gnomad_AC_XY \
            --plugin dbNSFP,~{vep_database}/dbNSFP/dbNSFP5.1a_grch37.gz,"${dbnsfp_str},${dbnsfp_gnomad_str}" \
            --plugin SpliceAI,snv=~{vep_database}/spliceAI/spliceai_scores.raw.snv.hg19.vcf.gz,indel=~{vep_database}/spliceAI/spliceai_scores.raw.indel.hg19.vcf.gz,cutoff=0.5 \
            --plugin FlankingSequence,10 \
            --plugin MissenseZscoreTranscript,~{vep_database}/GnomAD/missenseByTranscript.hg38.v4.1.bed \
            --fields "${cache_str},${custom_str},${dbnsfp_str},${dbnsfp_gnomad_str},${spliceai_str},${self_plugin_str}"
    >>>

    output {
        File annoVcf = "~{output_dir}/~{sample_id}.vep.anno.vcf"
    }

    runtime {
        cpus: threads
        container: "docker.io/ensemblorg/ensembl-vep:release_114.2"
        binding: "~{output_dir}:~{output_dir},~{vep_database}:~{vep_database},~{plugin_dir}:~{plugin_dir}"
    }
}

## SNPEff
### 用于融合注释
task SNPEff {
    input {
        String sample_id
        File vcf
        Int threads
        String genome
        String output_dir
        String snpeff_database
    }

    command <<<
        if [ ! -d ~{output_dir} ]; then
            mkdir -p ~{output_dir}
        fi

        java -jar /snpeff/snpEff.jar \
            -c /snpeff/snpEff.config \
            -v ~{genome} \
            ~{vcf} > ~{output_dir}/~{sample_id}.snpEff.anno.vcf
    >>>

    output {
        File annoVcf = "~{output_dir}/~{sample_id}.snpEff.anno.vcf"
    }

    runtime {
        cpus: threads
        container: "docker.io/staphb/snpeff:5.2f"
        binding: "~{output_dir}:~{output_dir},~{snpeff_database}:/snpeff/data"
    }
}

