#!/bin/bash
# 1. 可修改：指定样本文件夹所在的父目录
target_dir="/histor/kang/houli/sc_RNAseq/locust_brain/GDR25030573-中科院动物所12例动物脑10x单细胞转录组测序"
# 输出配置文件
ExpInformation="config.csv"

# 清空并写入表头
> "$ExpInformation"
echo -e "sample_name\tfile_path\tgroup" >> "$ExpInformation"

# 切换到目标目录遍历文件夹
for dir in "$target_dir"/*/; do
    # 去除路径前缀与末尾斜杠，得到纯样本名
    dir_name="${dir%/}"
    dir_name="${dir_name##*/}"

    sample_name="$dir_name"
    file_path="$target_dir/$dir_name"
    # 前两位字符作为分组
    group="${dir_name:0:2}"

    echo -e "$sample_name\t$file_path\t$group" >> "$ExpInformation"
done

echo "配置文件 $ExpInformation 生成成功！目标目录：$target_dir"
