# coding=utf-8
# pzw
# 20250813
# 使用MANE Bed和输入的Bed文件，取得交集，作为CNVkit的Target Bed
# https://github.com/pzweuj/ManeSelectBed

import argparse
import sys
from collections import defaultdict

class BedRecord:
    """BED文件记录类"""
    def __init__(self, chrom, start, end, **kwargs):
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
        self.extra_fields = kwargs
    
    def __str__(self):
        fields = [self.chrom, str(self.start), str(self.end)]
        for key in sorted(self.extra_fields.keys()):
            fields.append(str(self.extra_fields[key]))
        return '\t'.join(fields)


class BedIntersector:
    """BED文件交集处理类"""
    
    def __init__(self, use_chr=True):
        self.reference_intervals = defaultdict(list)
        self.use_chr = use_chr
    
    def normalize_chrom(self, chrom):
        """标准化染色体名称 - 参考文件总是有chr，根据use_chr参数决定输出格式"""
        if self.use_chr:
            # 输出需要chr前缀
            if not chrom.startswith('chr'):
                return f'chr{chrom}'
            return chrom
        else:
            # 输出不需要chr前缀
            if chrom.startswith('chr'):
                return chrom[3:]
            return chrom
    
    def chrom_to_ref_format(self, chrom):
        """将输入染色体名称转换为参考文件格式（总是有chr）"""
        if not chrom.startswith('chr'):
            return f'chr{chrom}'
        return chrom
    
    def load_reference_bed(self, reference_file):
        """加载参考BED文件"""
        print(f"Loading reference BED file: {reference_file}")
        
        with open(reference_file, 'r', encoding='utf-8') as f:
            header = f.readline().strip()
            if header.startswith('#'):
                # 解析表头
                header_fields = header[1:].split('\t')
                print(f"Reference BED header: {header_fields}")
            
            for line_num, line in enumerate(f, 2):
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                
                try:
                    fields = line.split('\t')
                    if len(fields) < 3:
                        print(f"Warning: Line {line_num} has insufficient columns, skipping")
                        continue
                    
                    chrom = fields[0]
                    start = int(fields[1])
                    end = int(fields[2])
                    
                    # 存储额外字段
                    extra_fields = {}
                    if len(fields) > 3:
                        field_names = ['location', 'symbol', 'refseq', 'ensembl', 'strand']
                        for i, field_name in enumerate(field_names):
                            if i + 3 < len(fields):
                                extra_fields[field_name] = fields[i + 3]
                    
                    record = BedRecord(chrom, start, end, **extra_fields)
                    self.reference_intervals[chrom].append(record)
                    
                except (ValueError, IndexError) as e:
                    print(f"Warning: Error parsing line {line_num}: {e}")
                    continue
        
        # 按起始位置排序每个染色体的区间
        for chrom in self.reference_intervals:
            self.reference_intervals[chrom].sort(key=lambda x: x.start)
        
        total_intervals = sum(len(intervals) for intervals in self.reference_intervals.values())
        print(f"Loaded {total_intervals} intervals from reference BED file")
    
    def find_overlaps(self, chrom, start, end):
        """查找与给定区间重叠的参考区间"""
        overlaps = []
        
        # 将输入染色体名称转换为参考文件格式（总是有chr）
        ref_chrom = self.chrom_to_ref_format(chrom)
        
        if ref_chrom not in self.reference_intervals:
            return overlaps
        
        for ref_record in self.reference_intervals[ref_chrom]:
            # 检查是否有重叠
            if ref_record.end <= start or ref_record.start >= end:
                continue
            
            # 计算重叠区间
            overlap_start = max(start, ref_record.start)
            overlap_end = min(end, ref_record.end)
            overlap_length = overlap_end - overlap_start
            
            overlaps.append({
                'ref_record': ref_record,
                'overlap_start': overlap_start,
                'overlap_end': overlap_end,
                'overlap_length': overlap_length
            })
        
        return overlaps
    
    def intersect_beds(self, input_file, output_file):
        """执行BED文件交集操作"""
        print(f"Processing input BED file: {input_file}")
        print(f"Output will be written to: {output_file}")
        
        intersected_count = 0
        total_input_count = 0
        
        with open(input_file, 'r', encoding='utf-8') as infile, \
             open(output_file, 'w', encoding='utf-8') as outfile:
            
            # 写入输出文件表头
            outfile.write('#chrom\tstart\tend\tname\n')
            
            for line_num, line in enumerate(infile, 1):
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                
                total_input_count += 1
                
                try:
                    fields = line.split('\t')
                    if len(fields) < 3:
                        print(f"Warning: Input line {line_num} has insufficient columns, skipping")
                        continue
                    
                    chrom = fields[0]
                    start = int(fields[1])
                    end = int(fields[2])
                    
                    # 保存输入文件的额外信息
                    input_info = '\t'.join(fields[3:]) if len(fields) > 3 else ''
                    
                    # 查找重叠
                    overlaps = self.find_overlaps(chrom, start, end)
                    
                    if overlaps:
                        intersected_count += 1
                        
                        # 对于每个重叠区间，输出注释后的结果
                        for overlap in overlaps:
                            ref_record = overlap['ref_record']
                            overlap_start = overlap['overlap_start']
                            overlap_end = overlap['overlap_end']
                            
                            # 构建name字段 (gene|location)
                            gene = ref_record.extra_fields.get('symbol', '')
                            location = ref_record.extra_fields.get('location', '')
                            name = f"{gene}|{location}"
                            
                            # 构建输出行 - 使用标准化的染色体名称
                            output_chrom = self.normalize_chrom(chrom)
                            output_fields = [
                                output_chrom,
                                str(overlap_start),
                                str(overlap_end),
                                name
                            ]
                            
                            outfile.write('\t'.join(output_fields) + '\n')
                    
                except (ValueError, IndexError) as e:
                    print(f"Warning: Error processing input line {line_num}: {e}")
                    continue
        
        print(f"Intersection complete:")
        print(f"  Total input intervals: {total_input_count}")
        print(f"  Intervals with overlaps: {intersected_count}")
        print(f"  Success rate: {intersected_count/total_input_count*100:.1f}%" if total_input_count > 0 else "  No input intervals processed")


def main():
    parser = argparse.ArgumentParser(
        description='BED文件交集工具 - 将输入BED文件与参考BED文件取交集并注释基因信息',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
使用示例:
  python mane_bed_intersect.py -r reference.bed -i input.bed -o output.bed
  python mane_bed_intersect.py -r reference.bed -i input.bed -o output.bed --chr
  
输入文件格式:
  - 参考BED文件: #chrom start end location symbol refseq ensembl strand (总是包含chr)
  - 输入BED文件: 标准BED格式 (chrom start end [其他字段]) (可能包含或不包含chr)
  
输出文件格式:
  - 简化的BED文件: chrom start end name (name格式为gene|location)
  - 使用--chr参数控制输出染色体是否包含chr前缀
        '''
    )
    
    parser.add_argument('-r', '--reference', required=True,
                        help='参考BED文件路径 (包含基因注释信息)')
    parser.add_argument('-i', '--input', required=True,
                        help='输入BED文件路径')
    parser.add_argument('-o', '--output', required=True,
                        help='输出BED文件路径')
    parser.add_argument('--chr', action='store_true', default=False,
                        help='输入和输出BED文件是否包含chr前缀 (参考文件总是包含chr)')
    
    args = parser.parse_args()
    
    try:
        # 创建交集处理器
        intersector = BedIntersector(use_chr=args.chr)
        
        # 加载参考BED文件
        intersector.load_reference_bed(args.reference)
        
        # 执行交集操作
        intersector.intersect_beds(args.input, args.output)
        
        print("BED文件交集操作完成!")
        
    except FileNotFoundError as e:
        print(f"错误: 文件未找到 - {e}")
        sys.exit(1)
    except Exception as e:
        print(f"错误: {e}")
        sys.exit(1)


if __name__ == '__main__':
    main()
