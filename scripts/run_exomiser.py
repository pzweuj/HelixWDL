# coding=utf-8
# Exomiser的运行器
# pzw
# 20250728
# 一并封装到exomiser的镜像里，因此调用的是镜像中的路径

import os
import yaml

# 读取默认的yaml文件，形成字典
def read_default_yaml_file(yaml_file_path="config/default.yaml"):
    """
    读取默认的YAML配置文件
    
    Args:
        yaml_file_path (str): YAML文件路径，默认为 config/default.yaml
    
    Returns:
        dict: 解析后的YAML内容字典
    
    Raises:
        FileNotFoundError: 当YAML文件不存在时
        yaml.YAMLError: 当YAML文件格式错误时
    """
    try:
        with open(yaml_file_path, 'r', encoding='utf-8') as file:
            config_data = yaml.safe_load(file)
            return config_data if config_data is not None else {}
    except FileNotFoundError:
        print(f"警告: 配置文件 {yaml_file_path} 不存在")
        return {}
    except yaml.YAMLError as e:
        print(f"错误: YAML文件解析失败 - {e}")
        return {}
    except Exception as e:
        print(f"错误: 读取配置文件时发生未知错误 - {e}")
        return {}

# 更新
def update_config_dict(input_dict, sample_id, vcf_file, genome=None, hpo_string=None, ped_file=None):
    input_dict['analysis']['vcf'] = os.path.abspath(vcf_file)
    if genome is not None:
        input_dict['analysis']['genomeAssembly'] = genome
    
    if hpo_string is not None:
        hpo_list = hpo_string.split(',')
        input_dict['analysis']['hpoIds'] = hpo_list

    if ped_file is not None:
        input_dict['analysis']['ped'] = ped_file_basename
        input_dict['analysis']['proband'] = sample_id
    
    running_dir = '/run_data'
    if not os.path.exists(running_dir):
        os.makedirs(running_dir)
    
    # 将字典输出为一个yaml文件
    output_yaml_path = os.path.join(running_dir, f"{sample_id}_exomiser_config.yaml")
    
    try:
        with open(output_yaml_path, 'w', encoding='utf-8') as yaml_file:
            yaml.dump(input_dict, yaml_file, default_flow_style=False, 
                     allow_unicode=True, indent=2, sort_keys=False)
        print(f"配置文件已生成: {output_yaml_path}")
        return output_yaml_path
    except Exception as e:
        print(f"错误: 生成YAML配置文件失败 - {e}")
        return None









