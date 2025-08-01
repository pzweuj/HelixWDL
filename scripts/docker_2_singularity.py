#!/usr/bin/env python3
# coding=utf-8
"""
Docker to Singularity Container Converter for WDL Scripts

This script converts container references in WDL files between Docker and Singularity formats.
It recursively processes all .wdl files in a specified directory.

Author: pzw
Date: 20250801
"""

import os
import argparse
import sys
from pathlib import Path
from typing import List, Optional


def docker_to_singularity(wdl_file: Path) -> bool:
    """
    将WDL文件中的Docker容器引用转换为Singularity格式
    使用URL编码风格的转换方案来避免字符冲突
    
    Args:
        wdl_file: WDL文件路径
        
    Returns:
        bool: 转换是否成功
    """
    try:
        temp_file = wdl_file.with_suffix('.tmp')
        
        with open(wdl_file, 'r', encoding='utf-8') as input_file, \
             open(temp_file, 'w', encoding='utf-8') as output_file:
            
            for line_num, line in enumerate(input_file, 1):
                if 'container:' in line:
                    # 处理容器行: container: "ghcr.io/pzweuj/mapping_tool:2025aug"
                    try:
                        # 提取容器名称
                        container_part = line.split('container:')[1].strip()
                        if container_part.startswith('"') and container_part.endswith('"'):
                            container_name = container_part[1:-1]
                        else:
                            container_name = container_part.strip('"')
                        
                        # 使用简洁的编码方案转换为Singularity格式
                        # 先保护原有下划线，然后转换路径和标签分隔符
                        # _ -> -U-
                        # / -> -S-  
                        # : -> -C-
                        sif_name = (container_name
                                   .replace('_', '-U-')
                                   .replace('/', '-S-')
                                   .replace(':', '-C-') + '.sif')
                        
                        new_line = line.split('container:')[0] + f'container: "{sif_name}"\n'
                        output_file.write(new_line)
                        
                    except Exception:
                        output_file.write(line)
                else:
                    output_file.write(line)
        
        # 替换原文件
        temp_file.replace(wdl_file)
        return True
        
    except Exception:
        if temp_file.exists():
            temp_file.unlink()
        return False


def singularity_to_docker(wdl_file: Path) -> bool:
    """
    将WDL文件中的Singularity容器引用转换为Docker格式
    使用URL编码风格的解码方案
    
    Args:
        wdl_file: WDL文件路径
        
    Returns:
        bool: 转换是否成功
    """
    try:
        temp_file = wdl_file.with_suffix('.tmp')
        
        with open(wdl_file, 'r', encoding='utf-8') as input_file, \
             open(temp_file, 'w', encoding='utf-8') as output_file:
            
            for line_num, line in enumerate(input_file, 1):
                if 'container:' in line and '.sif' in line:
                    try:
                        # 提取SIF文件名
                        container_part = line.split('container:')[1].strip()
                        if container_part.startswith('"') and container_part.endswith('"'):
                            sif_name = container_part[1:-1]
                        else:
                            sif_name = container_part.strip('"')
                        
                        # 转换为Docker格式
                        if sif_name.endswith('.sif'):
                            base_name = sif_name[:-4]  # 移除.sif后缀
                            
                            # 解码回原始Docker镜像名
                            # -C- -> :
                            # -S- -> /
                            # -U- -> _
                            docker_name = (base_name
                                         .replace('-C-', ':')
                                         .replace('-S-', '/')
                                         .replace('-U-', '_'))
                            
                            new_line = line.split('container:')[0] + f'container: "{docker_name}"\n'
                            output_file.write(new_line)
                        else:
                            output_file.write(line)
                            
                    except Exception:
                        output_file.write(line)
                else:
                    output_file.write(line)
        
        # 替换原文件
        temp_file.replace(wdl_file)
        return True
        
    except Exception:
        if temp_file.exists():
            temp_file.unlink()
        return False


def find_wdl_files(directory: Path) -> List[Path]:
    """
    递归查找目录中的所有WDL文件
    
    Args:
        directory: 搜索目录
        
    Returns:
        List[Path]: WDL文件路径列表
    """
    wdl_files = []
    for wdl_file in directory.rglob('*.wdl'):
        if wdl_file.is_file():
            wdl_files.append(wdl_file)
    
    return wdl_files


def process_directory(input_dir: Path, mode: str, file_pattern: Optional[str] = None) -> int:
    """
    处理目录中的WDL文件
    
    Args:
        input_dir: 输入目录
        mode: 转换模式 ('d2s' 或 's2d')
        file_pattern: 文件名模式过滤
        
    Returns:
        int: 成功处理的文件数量
    """
    if not input_dir.exists() or not input_dir.is_dir():
        return 0
    
    # 查找WDL文件
    wdl_files = find_wdl_files(input_dir)
    
    # 应用文件名过滤
    if file_pattern:
        wdl_files = [f for f in wdl_files if file_pattern in f.name]
    
    if not wdl_files:
        return 0
    
    success_count = 0
    for wdl_file in wdl_files:
        if mode == 'd2s':
            success = docker_to_singularity(wdl_file)
        elif mode == 's2d':
            success = singularity_to_docker(wdl_file)
        else:
            continue
            
        if success:
            success_count += 1
    
    return success_count

def main():
    """主函数"""
    parser = argparse.ArgumentParser(
        description='Convert container references in WDL files between Docker and Singularity formats',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Convert Docker to Singularity in workflows directory
  python docker_2_singularity.py -d workflows -m d2s
  
  # Convert Singularity to Docker in current directory
  python docker_2_singularity.py -d . -m s2d
  
  # Process only files containing 'cancer' in filename
  python docker_2_singularity.py -d workflows -m d2s -p cancer
  
  # Process with dry run
  python docker_2_singularity.py -d workflows -m d2s --dry-run
        """
    )
    
    parser.add_argument(
        '-d', '--directory',
        type=str,
        default='workflows',
        help='Input directory containing WDL files (default: workflows)'
    )
    
    parser.add_argument(
        '-m', '--mode',
        choices=['d2s', 's2d'],
        required=True,
        help='Conversion mode: d2s (Docker to Singularity) or s2d (Singularity to Docker)'
    )
    
    parser.add_argument(
        '-p', '--pattern',
        type=str,
        help='Filter files by name pattern (optional)'
    )
    

    
    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='Show what would be processed without making changes'
    )
    
    args = parser.parse_args()
    
    input_dir = Path(args.directory)
    
    if args.dry_run:
        wdl_files = find_wdl_files(input_dir)
        if args.pattern:
            wdl_files = [f for f in wdl_files if args.pattern in f.name]
        
        print(f"Would process {len(wdl_files)} files:")
        for wdl_file in wdl_files:
            print(f"  {wdl_file}")
        return 0
    
    # 处理文件
    try:
        success_count = process_directory(input_dir, args.mode, args.pattern)
        return 0 if success_count > 0 else 1
            
    except KeyboardInterrupt:
        return 1
    except Exception:
        return 1


if __name__ == '__main__':
    sys.exit(main())
