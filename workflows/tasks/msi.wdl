version 2.0
# 微卫星分析

## MSISensor-Pro
### https://github.com/xjtu-omics/msisensor-pro
### 建立库
task MSISensorProList {
    input {
        String prefix
        String output_dir
        File reference
    }

    command <<<
        if [ ! -d ~{output_dir} ]; then
            mkdir -p ~{output_dir}
        fi

        msisensor-pro scan -d ~{reference} -o ~{output_dir}/~{prefix}.msi.list
    >>>

    output {
        File msi_list = "~{output_dir}/~{prefix}.msi.list"
    }

    runtime {
        container: "ghcr.io/pzweuj/msisensor-pro:v1.3.0"
        binding: "~{output_dir}:~{output_dir}"
    }
}

### 建立基线
### 需要把所有的正常样本bam放置于一个文件夹里
### 建议用20个样本以上建立
task MSISensorProBaseline {
    input {
        String prefix
        String output_dir
        File msi_list
        String input_dir
    }

    command <<<
        if [ ! -d ~{output_dir} ]; then
            mkdir -p ~{output_dir}
        fi

        find ~{input_dir} -name "*.bam" -exec realpath {} \; | awk -F/ '{file=$NF; gsub(/\.bam$/,"",file); print file "\t" $0}' > sample.txt

        msisensor-pro baseline -d ~{msi_list} -i sample.txt -o ~{output_dir}/~{prefix}.msi.ref.txt
    >>>

    output {
        File msi_ref = "~{output_dir}/~{prefix}.msi.ref.txt"
    }

    runtime {
        container: "ghcr.io/pzweuj/msisensor-pro:v1.3.0"
        binding: "~{output_dir}:~{output_dir},~{input_dir}:~{input_dir}"
    }
}

### 单样本模式
task MSISensorProTumorOnly {
    input {
        String sample_id
        String output_dir
        File bam
        File bai
        File baseline
    }

    command <<<
        if [ ! -d ~{output_dir} ]; then
            mkdir -p ~{output_dir}
        fi

        msisensor-pro pro -d ~{baseline} -t ~{bam} -o ~{output_dir}/~{sample_id}.msi.txt
    >>>

    output {
        File msi_result = "~{output_dir}/~{sample_id}.msi.txt"
    }

    runtime {
        container: "ghcr.io/pzweuj/msisensor-pro:v1.3.0"
        binding: "~{output_dir}:~{output_dir}"
    }
}

### 配对模式
task MSISensorProPair {
    input {
        String sample_id
        String output_dir
        File bam
        File bai
        File normal_bam
        File normal_bai
        File msi_list
    }

    command <<<
        if [ ! -d ~{output_dir} ]; then
            mkdir -p ~{output_dir}
        fi

        msisensor-pro msi -d ~{msi_list} -t ~{bam} -n ~{normal_bam} -o ~{output_dir}/~{sample_id}.msi.pair.txt
    >>>

    output {
        File msi_result = "~{output_dir}/~{sample_id}.msi.pair.txt"
    }

    runtime {
        container: "ghcr.io/pzweuj/msisensor-pro:v1.3.0"
        binding: "~{output_dir}:~{output_dir}"
    }
}
