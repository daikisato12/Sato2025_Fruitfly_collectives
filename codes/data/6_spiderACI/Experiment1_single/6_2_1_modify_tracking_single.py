import sys
import pandas as pd
import numpy as np

# データを読み込む
infile = sys.argv[1]  # 入力ファイル例: "m1_1_1_12.0s_100px_6.0mmps_track.tsv"
outfile = sys.argv[2]  # 出力ファイル例: "m1_1_1_12.0s_100px_6.0mmps_track_modified.tsv"
total_duration = float(sys.argv[3])  # トラッキングデータ全体の時間（秒）
cut_duration = float(sys.argv[4])    # 出力したい後半部分の時間（秒）

# トラッキングデータを読み込む
data = pd.read_csv(infile, sep='\t')

# `Stopped` 列を3行下にシフト
if 'Stopped' in data.columns:
    data['Stopped'] = data['Stopped'].shift(3)

# 下にずらした結果として追加された最後の3行は無視
data = data.iloc[:-3]

# `Time` 列と `Time_interval` 列があるかを確認して適切なフレームでトリミング
if 'Time' in data.columns and 'Time_interval' in data.columns:
    # データの後半 `cut_duration` 秒のフレームを取得
    end_time = data['Time'].iloc[-1]
    start_time = end_time - cut_duration
    
    # `Time` 列の値が `start_time` 以上のフレームを抽出
    data = data[data['Time'] >= start_time].reset_index(drop=True)
else:
    # `Time` 列と `Time_interval` 列がない場合は、フレームレートに基づいて計算
    fps = len(data) / total_duration  # フレームレートを計算
    frames_to_keep = int(cut_duration * fps)  # 切り取りたい秒数に基づく必要なフレーム数
    start_frame = len(data) - frames_to_keep  # 後半部分の開始フレーム位置
    data = data.iloc[start_frame:].reset_index(drop=True)
    data['Time_interval'] = cut_duration / len(data)

data['Time'] = data['Time_interval'].cumsum() - data['Time_interval'].iloc[0]
# 必ず `Frame`, `Time`, `Time_interval` が最初の3列になるようにカラム順を指定
columns_order = ['Frame', 'Time', 'Time_interval'] + [col for col in data.columns if col not in ['Frame', 'Time', 'Time_interval']]
data = data[columns_order]

# `Spider_Direction_Radians` 列を削除
if 'Spider_Direction_Radians' in data.columns:
    data = data.drop(columns=['Spider_Direction_Radians'])

# 前回の有効なクモの位置を保持する変数を初期化
last_valid_spider_x, last_valid_spider_y = None, None

# 入れ替わりと見なすための閾値
threshold = 10  # 同一視する位置差（調整可能）

# 各行を処理して座標を修正し、Distanceを再計算
for i in range(len(data)):
    # 現在のハエとクモの位置
    fly_x, fly_y = data.loc[i, 'Fly_X'], data.loc[i, 'Fly_Y']
    spider_x, spider_y = data.loc[i, 'Spider_X'], data.loc[i, 'Spider_Y']
    
    # Step 1: クモの位置が (500, 0) の場合、直前の位置で置き換える
    if spider_x == 500 and spider_y == 0:
        if last_valid_spider_x is not None and last_valid_spider_y is not None:
            data.loc[i, 'Spider_X'] = last_valid_spider_x
            data.loc[i, 'Spider_Y'] = last_valid_spider_y

    # Step 2: ハエとクモの位置が閾値内で入れ替わっているか確認
    if i > 0:
        prev_fly_x, prev_fly_y = data.loc[i - 1, 'Fly_X'], data.loc[i - 1, 'Fly_Y']
        prev_spider_x, prev_spider_y = data.loc[i - 1, 'Spider_X'], data.loc[i - 1, 'Spider_Y']
        
        # 位置が閾値内で入れ替わっているか確認
        if (abs(prev_fly_x - prev_spider_x) > threshold and abs(prev_fly_y - prev_spider_y) > threshold) and \
            (abs(spider_x - fly_x) <= threshold and abs(spider_y - fly_y) <= threshold):
            # 位置を元に戻す
            data.loc[i, 'Spider_X'], data.loc[i, 'Spider_Y'] = prev_spider_x, prev_spider_y
    
    # Step 3: クモの位置が (500, 0) の場合、NAにする
    if data.loc[i, 'Spider_X'] == 500 and data.loc[i, 'Spider_Y'] == 0:
        data.loc[i, 'Spider_X'] = np.nan
        data.loc[i, 'Spider_Y'] = np.nan
 
    last_valid_spider_x, last_valid_spider_y = data.loc[i, 'Spider_X'], data.loc[i, 'Spider_Y']

    # ハエとクモの距離を計算してDistance列を更新
    data.loc[i, 'Distance'] = np.sqrt((data.loc[i, 'Fly_X'] - data.loc[i, 'Spider_X'])**2 +
                                      (data.loc[i, 'Fly_Y'] - data.loc[i, 'Spider_Y'])**2)

# 修正後のデータを保存
data.to_csv(outfile, sep='\t', index=False)