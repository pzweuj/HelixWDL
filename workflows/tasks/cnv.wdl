version 2.0
# CNV分析脚本

## 使用基于Mane转录本的bed || 注，脚本未编写
task ManeBedIntersect {
    input {
        String prefix
        String output_dir
        File bed
        File mane_bed
    }

    command <<<
        if [ ! -d ~{output_dir} ]; then
            mkdir -p ~{output_dir}
        fi

        python3 /opt/HelixWDL/scripts/mane_bed_intersect.py -i ~{bed} -m ~{mane_bed} -o ~{output_dir}/~{prefix}.cnv.bed
    >>>

    output {
        File cnv_bed = "~{output_dir}/~{prefix}.cnv.bed"
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

        # 如果提供了refflat文件，则执行target步骤
        ~{if defined(refflat) then "cnvkit.py target " + bed + " --annotate " + refflat + " -o " + output_dir + "/" + prefix + ".target.bed" else "cp " + bed + " " + output_dir + "/" + prefix + ".target.bed"}

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










