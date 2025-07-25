version 2.0
# 结构定义脚本

struct BwaIndex {
    File fasta
    File fai
    Array[File] index_file
}

struct GATKIndex {
    File fasta
    File fai
    File dict
}

struct IndexBundle {
    File fasta
    Array[File] index_file
}
