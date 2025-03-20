#!/bin/bash

# ディレクトリの指定
input_dir="../../data/Experiment2_group/1_rawdata/tracks"
output_dir="../../data/Experiment2_group/2_moddata/tracks"
file=$1
# sh 2_2_modify_tracking_group.sh 20250121

# 出力ディレクトリが存在しない場合は作成
mkdir -p "$output_dir"

# 各ファイルを処理
for input_file in "$input_dir"/$file*.tsv; do
    # ファイル名の取得と拡張子の削除
    filename=$(basename "$input_file" .tsv)
    
    # 出力ファイル名の作成
    output_file="$output_dir/${filename}_modified.tsv"
    
    # Pythonスクリプトの実行
    python3 2_2_modify_tracking_group.py $input_file $output_file 360 300
done