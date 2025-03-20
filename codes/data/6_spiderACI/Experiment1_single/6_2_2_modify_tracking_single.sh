#!/bin/bash

# ディレクトリの指定
input_dir="../../../../data/6_spiderACI/Experiment1_single/rawdata/tracks"
output_dir="../../../../data/6_spiderACI/Experiment1_single/moddata/tracks"
file=$1
# sh 2_modify_tracking.sh 20241102

# 出力ディレクトリが存在しない場合は作成
mkdir -p "$output_dir"

# 各ファイルを処理
for input_file in "$input_dir"/$file*.tsv; do
    # ファイル名の取得と拡張子の削除
    filename=$(basename "$input_file" .tsv)
    
    # 出力ファイル名の作成
    output_file="$output_dir/${filename}_modified.tsv"
    
    # Pythonスクリプトの実行
    python3 6_2_1_modify_tracking_single.py "$input_file" "$output_file" 360 300
done