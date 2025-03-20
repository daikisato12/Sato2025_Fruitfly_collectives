import sys
import pandas as pd
import numpy as np

# データを読み込む
infile = sys.argv[1]  # 入力ファイル例: "m1_1_1_6-6-6-6-6-6s_50px_24.0mmps_weak_track.tsv"
outfile = sys.argv[2]  # 出力ファイル例: "m1_1_1_6-6-6-6-6-6s_50px_24.0mmps_weak_track_modified.tsv"
total_duration = float(sys.argv[3])  # トラッキングデータ全体の時間（秒）
cut_duration = float(sys.argv[4])    # 出力したい後半部分の時間（秒）

# トラッキングデータを読み込む
data = pd.read_csv(infile, sep='\t')

# `Stopped` 列を各ハエごとに3行下にシフト
if 'Stopped' in data.columns:
    data['Stopped'] = data.groupby('Fly_ID')['Stopped'].shift(3)

# 下にずらした結果として追加された最後の3行を削除
data = data.dropna(subset=['Stopped'], how='all')

# `Time` 列と `Time_interval` 列があるかを確認して適切なフレームでトリミング
if 'Time' in data.columns and 'Time_interval' in data.columns:
    # 各ハエごとに処理
    def process_fly_group(group):
        # データの後半 `cut_duration` 秒のフレームを取得
        end_time = group['Time'].iloc[-1]
        start_time = end_time - cut_duration
        
        # `Time` 列の値が `start_time` 以上のフレームを抽出
        trimmed = group[group['Time'] >= start_time].copy()
        trimmed.reset_index(drop=True, inplace=True)
        
        # `Time` 列を再計算
        trimmed['Time'] = trimmed['Time_interval'].cumsum() - trimmed['Time_interval'].iloc[0]
        
        return trimmed

    data = data.groupby('Fly_ID').apply(process_fly_group).reset_index(drop=True)
else:
    # `Time` 列と `Time_interval` 列がない場合は、フレームレートに基づいて計算
    # 各ハエごとに処理
    def process_fly_group_no_time(group):
        fps = len(group) / total_duration  # フレームレートを計算
        frames_to_keep = int(cut_duration * fps)  # 切り取りたい秒数に基づく必要なフレーム数
        start_frame = max(len(group) - frames_to_keep, 0)  # 後半部分の開始フレーム位置
        trimmed = group.iloc[start_frame:].copy()
        trimmed.reset_index(drop=True, inplace=True)
        trimmed['Time_interval'] = cut_duration / len(trimmed)
        trimmed['Time'] = trimmed['Time_interval'].cumsum() - trimmed['Time_interval'].iloc[0]
        return trimmed

    data = data.groupby('Fly_ID').apply(process_fly_group_no_time).reset_index(drop=True)

# `Spider_Direction_Radians` 列を削除
if 'Spider_Direction_Radians' in data.columns:
    data = data.drop(columns=['Spider_Direction_Radians'])

# クモの位置をフレームごとに修正
# フレームをソート
data = data.sort_values(by=['Frame', 'Fly_ID']).reset_index(drop=True)

# クモの位置を保持する変数を初期化
last_valid_spider_x, last_valid_spider_y = None, None

# 入れ替わりと見なすための閾値
threshold = 10  # 同一視する位置差（調整可能）

# フレームごとに処理
for frame in sorted(data['Frame'].unique()):
    frame_data = data[data['Frame'] == frame]
    
    # クモの位置が (500, 0) の行が存在するか確認
    invalid_spider = frame_data[(frame_data['Spider_X'] == 500) & (frame_data['Spider_Y'] == 0)]
    
    if not invalid_spider.empty:
        # クモの位置を前回の有効な位置に置き換える
        if last_valid_spider_x is not None and last_valid_spider_y is not None:
            data.loc[data['Frame'] == frame, 'Spider_X'] = last_valid_spider_x
            data.loc[data['Frame'] == frame, 'Spider_Y'] = last_valid_spider_y
    else:
        # クモの位置を更新
        # クモの位置が全てNAの場合はスキップ
        valid_spider = frame_data.dropna(subset=['Spider_X', 'Spider_Y'])
        if not valid_spider.empty:
            # ここでは最初の有効なクモの位置を使用
            current_spider_x = valid_spider['Spider_X'].iloc[0]
            current_spider_y = valid_spider['Spider_Y'].iloc[0]
            last_valid_spider_x, last_valid_spider_y = current_spider_x, current_spider_y
        else:
            # クモの位置がNAの場合は何もしない
            pass
    
    # フレーム内のハエとクモの位置が入れ替わっているか確認
    if frame > data['Frame'].min():
        prev_frame_data = data[data['Frame'] == frame - 1]
        current_frame_data = data[data['Frame'] == frame]
        
        # クモの前の位置
        prev_spider = prev_frame_data[['Spider_X', 'Spider_Y']].dropna().drop_duplicates()
        if not prev_spider.empty:
            prev_spider_x, prev_spider_y = prev_spider.iloc[0]
            
            # 現在のクモの位置
            current_spider = current_frame_data[['Spider_X', 'Spider_Y']].dropna().drop_duplicates()
            if not current_spider.empty:
                current_spider_x, current_spider_y = current_spider.iloc[0]
                
                # ハエの位置との入れ替わりをチェック
                for _, row in current_frame_data.iterrows():
                    fly_x, fly_y = row['Fly_X'], row['Fly_Y']
                    distance_to_prev_spider = np.sqrt((fly_x - prev_spider_x)**2 + (fly_y - prev_spider_y)**2)
                    distance_to_current_spider = np.sqrt((fly_x - current_spider_x)**2 + (fly_y - current_spider_y)**2)
                    
                    if distance_to_prev_spider <= threshold and distance_to_current_spider <= threshold:
                        # 位置を元に戻す
                        data.loc[(data['Frame'] == frame) & (data['Fly_ID'] == row['Fly_ID']), ['Spider_X', 'Spider_Y']] = [prev_spider_x, prev_spider_y]
    
    # クモの位置が (500, 0) の場合、NAにする
    data.loc[(data['Frame'] == frame) & (data['Spider_X'] == 500) & (data['Spider_Y'] == 0), ['Spider_X', 'Spider_Y']] = [np.nan, np.nan]
    
    # 有効なクモの位置を更新
    frame_valid_spider = data[(data['Frame'] == frame) & (~data['Spider_X'].isna()) & (~data['Spider_Y'].isna())]
    if not frame_valid_spider.empty:
        # 最初の有効なクモの位置を保持
        last_valid_spider_x = frame_valid_spider['Spider_X'].iloc[0]
        last_valid_spider_y = frame_valid_spider['Spider_Y'].iloc[0]

    # ハエとクモの距離を計算してDistance列を更新
    data.loc[data['Frame'] == frame, 'Distance'] = data[data['Frame'] == frame].apply(
        lambda row: np.sqrt((row['Fly_X'] - row['Spider_X'])**2 + (row['Fly_Y'] - row['Spider_Y'])**2) 
                     if not pd.isna(row['Spider_X']) and not pd.isna(row['Spider_Y']) else np.nan,
        axis=1
    )

# 最終的なデータを保存
data.to_csv(outfile, sep='\t', index=False)
