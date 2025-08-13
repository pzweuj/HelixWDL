version 2.0
# CNV分析脚本

## 使用基于Mane转录本的bed || 注，脚本未编写
task ManeBedIntersect {
    input {
        String prefix
        String output_dir
        File bed
        File mane_bed
        Boolean use_chr = true
    }

    command <<<
        if [ ! -d ~{output_dir} ]; then
            mkdir -p ~{output_dir}
        fi

        python3 /opt/HelixWDL/scripts/mane_bed_intersect.py \
            -i ~{bed} \
            -r ~{mane_bed} \
            -o ~{output_dir}/~{prefix}.cnv.target.bed \
            ~{if use_chr then "--chr" else ""}
    >>>

    output {
        File cnv_bed = "~{output_dir}/~{prefix}.cnv.target.bed"
    }

    runtime {
        container: "ghcr.io/pzweuj/helixwdl:v0.0.1"
        binding: "~{output_dir}:~{output_dir}"
    }
}

## CNVkit路线
### Target文件及AntiTarget文件生成
### exclude_bed_hg19: http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeMapability/wgEncodeDukeMapabilityRegionsExcludable.bed.gz
### exclude_bed_hg38: https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz
### refFlat_hg19: https://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/refFlat.txt.gz
### refFlat_hg38: https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refFlat.txt.gz
task CNVkitTargetAndAntitarget {
    input {
        String prefix
        String output_dir
        File reference
        File bed
        File exclude_bed
        File? refflat
    }

    command <<<
        if [ ! -d ~{output_dir} ]; then
            mkdir -p ~{output_dir}
        fi

        # 生成access文件
        cnvkit.py access ~{reference} \
            -x ~{exclude_bed} \
            -o ~{output_dir}/~{prefix}.access.bed

        # 如果提供了refflat文件，则执行target步骤；否则直接使用输入的bed文件作为target
        if [ -n "~{refflat}" ]; then
            cnvkit.py target ~{bed} --short-names --annotate ~{refflat} -o ~{output_dir}/~{prefix}.target.bed
        else
            # 当没有refflat文件时，直接将输入的bed文件重命名为target.bed供antitarget步骤使用
            cp ~{bed} ~{output_dir}/~{prefix}.target.bed
        fi

        # 生成antitarget文件
        cnvkit.py antitarget ~{output_dir}/~{prefix}.target.bed \
            -g ~{output_dir}/~{prefix}.access.bed \
            -o ~{output_dir}/~{prefix}.antitarget.bed
    >>>

    output {
        File access_bed = "~{output_dir}/~{prefix}.access.bed"
        File target_bed = "~{output_dir}/~{prefix}.target.bed"
        File antitarget_bed = "~{output_dir}/~{prefix}.antitarget.bed"
    }

    runtime {
        container: "ghcr.io/pzweuj/cnvkit:v0.9.11.p4"
        binding: "~{output_dir}:~{output_dir}"
    }
}

## 获得单个样本的coverage
task CNVkitCoverage {
    input {
        String sample_id
        String output_dir
        File bam
        File bai
        File target_bed
        File? antitarget_bed
    }

    command <<<
        if [ ! -d ~{output_dir} ]; then
            mkdir -p ~{output_dir}
        fi

        cnvkit.py coverage ~{bam} ~{target_bed} -o ~{output_dir}/~{sample_id}.targetcoverage.cnn
        
        ~{if defined(antitarget_bed) then 
            "cnvkit.py coverage " + bam + " " + antitarget_bed + " -o " + output_dir + "/" + sample_id + ".antitargetcoverage.cnn"
        else 
            "touch " + output_dir + "/" + sample_id + ".antitargetcoverage.cnn"
        }
    >>>

    output {
        File target_coverage = "~{output_dir}/~{sample_id}.targetcoverage.cnn"
        File antitarget_coverage = "~{output_dir}/~{sample_id}.antitargetcoverage.cnn"
    }

    runtime {
        container: "ghcr.io/pzweuj/cnvkit:v0.9.11.p4"
        binding: "~{output_dir}:~{output_dir}"
    }
}

## 建立基线
task CNVkitReference {
    input {
        String prefix
        String coverage_dir
        String output_dir
        File reference
    }

    command <<<
        if [ ! -d ~{output_dir} ]; then
            mkdir -p ~{output_dir}
        fi

        cnvkit.py reference ~{coverage_dir}/*.{,anti}targetcoverage.cnn \
            --fasta ~{reference} -o ~{output_dir}/~{prefix}.ref.cnn
    >>>

    output {
        File ref_cnn = "~{output_dir}/~{prefix}.ref.cnn"
    }

    runtime {
        container: "ghcr.io/pzweuj/cnvkit:v0.9.11.p4"
        binding: "~{output_dir}:~{output_dir},~{coverage_dir}:~{coverage_dir}"
    }
}

## 单样本的CNV分析
task CNVkit {
    input {
        String sample_id
        String output_dir
        File target_cnn
        File antitarget_cnn
        File cnv_baseline
        String method
    }

    command <<<
        if [ ! -d ~{output_dir} ]; then
            mkdir -p ~{output_dir}
        fi

        cnvkit.py fix ~{target_cnn} ~{antitarget_cnn} ~{cnv_baseline} -o ~{output_dir}/~{sample_id}.cnr
        cnvkit.py segment ~{output_dir}/~{sample_id}.cnr -o ~{output_dir}/~{sample_id}.segment.cnr -m ~{method}
        cnvkit.py call ~{output_dir}/~{sample_id}.cnr -o ~{output_dir}/~{sample_id}.cns
        cnvkit.py call ~{output_dir}/~{sample_id}.segment.cnr -o ~{output_dir}/~{sample_id}.segment.cns
        cnvkit.py scatter ~{output_dir}/~{sample_id}.cnr -s ~{output_dir}/~{sample_id}.cns -o ~{output_dir}/~{sample_id}.scatter.pdf
        cnvkit.py diagram ~{output_dir}/~{sample_id}.cnr -s ~{output_dir}/~{sample_id}.cns -o ~{output_dir}/~{sample_id}.diagram.pdf
    >>>

    output {
        File cns = "~{output_dir}/~{sample_id}.cns"
        File seg_cns = "~{output_dir}/~{sample_id}.segment.cns"
        File scatter_pdf = "~{output_dir}/~{sample_id}.scatter.pdf"
        File diagram_pdf = "~{output_dir}/~{sample_id}.diagram.pdf"
    }

    runtime {
        container: "ghcr.io/pzweuj/cnvkit:v0.9.11.p4"
        binding: "~{output_dir}:~{output_dir}"
    }
}

## CNV结果整理
task CNVFormat {
    input {
        String sample_id
        String output_dir
        File cnvkit_cnr
    }
}

