version 2.0
# WES相关高级分析

## ExpansionHunter
## https://github.com/Illumina/ExpansionHunter
task ExpansionHunter {
    input {
        String sample_id
        String output_dir
        File bam
        File bai
        File reference
        File fai
        Int threads
        String sex = "female"
    }

    command <<<
        if [ ! -d ~{output_dir} ]; then
            mkdir -p ~{output_dir}
        fi

        y_info=$(awk '$1 ~ /^(chr)?Y$/ {print $1, $2}' ~{fai})
        if [ -z "$y_info" ]; then
            echo "Error: Y chromosome not found in reference" >&2
            exit 1
        fi
        
        y_chr=$(echo $y_info | cut -d' ' -f1)
        y_length=$(echo $y_info | cut -d' ' -f2)

        # 根据Y染色体长度和命名方式判断基因组版本
        case $y_length in
            57227415)
                # GRCh38/hg38
                if [[ $y_chr == "chrY" ]]; then
                    genome_version="hg38"
                else
                    genome_version="grch38"
                fi
                ;;
            59373566)
                # GRCh37/hg19
                if [[ $y_chr == "chrY" ]]; then
                    genome_version="hg19"
                else
                    genome_version="grch37"
                fi
                ;;
            *)
                echo "Warning: Unknown genome version (Y length: $y_length), using default coordinates" >&2
                if [[ $y_chr == "chrY" ]]; then
                    genome_version="hg19"
                else
                    genome_version="grch37"
                fi
                ;;
        esac

        # 开发中

        ExpansionHunter \
            --reads ~{bam} \
            --reference ~{reference} \
            --variant-catalog /opt/ExpansionHunter/variant_catalog/${genome_version}/variant_catalog.json \
            --output-prefix ~{output_dir}/~{sample_id}.eh \
            -n ~{threads} \
            --sex ~{sex}
     >>>

    output {
        File str_json = "~{output_dir}/~{sample_id}.eh.json"
        File str_vcf = "~{output_dir}/~{sample_id}.eh.vcf"
        File str_bam = "~{output_dir}/~{sample_id}.eh_realigned.bam"
        File str_result = "~{output_dir}/~{sample_id}.ExpansionHunter.txt"
    }

    runtime {
        cpus: threads
        container: "docker.io/zihhuafang/expansionhunter:v5.0.0"
        binding: "~{output_dir}:~{output_dir}"
    }
}

## ROH分析
### https://github.com/mquinodo/AutoMap
task ROH {
    input {
        String sample_id
        File vcf
        String output_dir
        File fai
    }

    command <<<
        if [ ! -d ~{output_dir} ]; then
            mkdir -p ~{output_dir}
        fi

        y_info=$(awk '$1 ~ /^(chr)?Y$/ {print $1, $2}' ~{fai})
        if [ -z "$y_info" ]; then
            echo "Error: Y chromosome not found in reference" >&2
            exit 1
        fi
        
        y_length=$(echo $y_info | cut -d' ' -f2)

        # 根据Y染色体长度判断基因组版本（ROH只需要hg19和hg38）
        case $y_length in
            57227415)
                # hg38
                genome_version="hg38"
                ;;
            59373566)
                # hg19
                genome_version="hg19"
                ;;
            *)
                echo "Warning: Unknown genome version (Y length: $y_length), using hg19 as default" >&2
                genome_version="hg19"
                ;;
        esac

        echo -e "#Chr\tBegin\tEnd\tSize(Mb)\tNb_variants\tPercentage_homozygosity\n" > ~{output_dir}/~{sample_id}.ROH.txt

        bash /opt/AutoMap/AutoMap_v1.3.sh \
            --vcf ~{vcf} \
            --genome ${genome_version} \
            --out ~{output_dir} \
            --id ~{sample_id}
        
        if [ -f ~{output_dir}/~{sample_id}.HomRegions.tsv ]; then
            cat ~{output_dir}/~{sample_id}.HomRegions.tsv | grep -v '##' > ~{output_dir}/~{sample_id}.ROH.txt
        fi
    >>>

    output {
        File roh_result = "~{output_dir}/~{sample_id}.ROH.txt"
    }

    runtime {
        container: "ghcr.io/pzweuj/automap:v1.3.p7"
        binding: "~{output_dir}:~{output_dir}"
        cpus: 2
    }
}