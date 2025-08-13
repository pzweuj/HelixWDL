version 2.0
# 微卫星基线建立流程

import "./tasks/msi.wdl" as msi

workflow MSISensorProBaseline {
    input {
        String prefix
        String output_dir
        File reference
        String input_dir
    }

    call msi.MSISensorProList as MSISensorProList {input: prefix=prefix, output_dir=output_dir, reference=reference}
    call msi.MSISensorProBaseline as MSISensorProBaseline {input: prefix=prefix, output_dir=output_dir, msi_list=MSISensorProList.msi_list, input_dir=input_dir}

    output {
        File msi_ref = MSISensorProBaseline.msi_ref
    }
}
