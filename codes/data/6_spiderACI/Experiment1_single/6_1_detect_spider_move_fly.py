import sys
sys.path.append("/home/lab-macbook/LabExperiment/Codes/TIS/SnapAnImageAndConverttoOpenCV/")

import os
import cv2
import numpy as np
import math
import random
import time
import TIS
import csv
from datetime import datetime
from screeninfo import get_monitors
from collections import deque

# コマンドライン引数の処理
if len(sys.argv) < 9:
    print("Usage: python 6_1_detect_spider_move_fly.py <output_directory> <duration_seconds> <stop_duration> <distance_threshold> <fly_speed> <individual_id> <sex> <trial_id>")
    sys.exit(1)

# 引数からパラメータの取得
output_base_dir = sys.argv[1]
duration_seconds = int(sys.argv[2])
stop_duration = float(sys.argv[3])
distance_threshold = int(sys.argv[4])
average_speed_mm = float(sys.argv[5])
individual_id = sys.argv[6]
sex = sys.argv[7]
trial_id = sys.argv[8]

# 出力ディレクトリの生成（ディレクトリ名に引数を含める）
start_time_str = datetime.now().strftime("%Y%m%d%H%M%S")
output_dir_arena = os.path.join(output_base_dir, f"images/{start_time_str}_{duration_seconds}s_{stop_duration}s_{distance_threshold}px_{average_speed_mm}mmps_{individual_id}_{sex}_{trial_id}/arena")
output_dir_tracking = os.path.join(output_base_dir, f"images/{start_time_str}_{duration_seconds}s_{stop_duration}s_{distance_threshold}px_{average_speed_mm}mmps_{individual_id}_{sex}_{trial_id}/tracking")
os.makedirs(output_dir_arena, exist_ok=True)
os.makedirs(output_dir_tracking, exist_ok=True)

# TSVファイルのパス設定
output_file = os.path.join(output_base_dir, f"tracks/{start_time_str}_{duration_seconds}s_{stop_duration}s_{distance_threshold}px_{average_speed_mm}mmps_{individual_id}_{sex}_{trial_id}_track.tsv")

# ディスプレイ設定とその他の基本設定
monitors = get_monitors()
main_display = monitors[0]
sub_display = monitors[1] if len(monitors) > 1 else monitors[0]
camera_width, camera_height = 640, 480
arena_diameter_mm = 84
desired_arena_diameter_px = 600
scale_px_per_mm = desired_arena_diameter_px / arena_diameter_mm
arena_radius_px = desired_arena_diameter_px // 2 - 10
camera_arena_scale = camera_height / (desired_arena_diameter_px)
horizontal_margin = 0
vertical_margin = 0
spider_center_x, spider_center_y = 500, 0
fps = 60
pixel_speed = (average_speed_mm * scale_px_per_mm) / fps
must_move_until = 0

# アリーナの中心位置（ウィンドウ内での固定位置）を指定
arena_center_x = 617
arena_center_y = 556

# スケール補正係数
vertical_scale_correction = 0.9
horizontal_scale_correction = 0.9

# TISカメラの設定
Tis = TIS.TIS()
Tis.openDevice("28020355", camera_width, camera_height, "60/1", TIS.SinkFormats.GRAY8, True)
Tis.Start_pipeline()
Tis.Set_Property("Gain Auto", True)
Tis.Set_Property("Exposure Auto", False)
Tis.Set_Property("Exposure Time (us)", 100)

# 背景の初期化
if Tis.Snap_image(0.1):
    background = Tis.Get_image()
    background = cv2.resize(background, (camera_width, camera_height))
else:
    print("Failed to capture background image.")
    background = np.zeros((camera_height, camera_width), dtype=np.uint8)

# 出力用TSVファイルの初期化
with open(output_file, mode='w', newline='') as file:
    writer = csv.writer(file, delimiter='\t')
    writer.writerow(["Frame", "Time_interval", "Time", "Fly_X", "Fly_Y", "Spider_X", "Spider_Y", "Distance", "Stopped", "Fly_Direction_Radians", "Spider_Direction_Radians"])

# 距離を計算する関数
def Calc_distance(center1, center2):
    return np.sqrt((center2[0] - center1[0]) ** 2 + (center2[1] - center1[1]) ** 2)

# クモの向きを計算する関数
def calculate_spider_direction(contour):
    approx = cv2.approxPolyDP(contour, 0.04 * cv2.arcLength(contour, True), True)
    edges = [(approx[i][0], approx[j][0]) for i in range(len(approx)) for j in range(i+1, len(approx))]
    sorted_edges = sorted(edges, key=lambda pair: Calc_distance(pair[0], pair[1]), reverse=True)
    
    long_edge = sorted_edges[0]
    long_edge_midpoint = ((long_edge[0][0] + long_edge[1][0]) / 2, (long_edge[0][1] + long_edge[1][1]) / 2)
    
    if len(approx) == 3:
        opposite_vertex = [pt[0] for pt in approx if not np.array_equal(pt[0], long_edge[0]) and not np.array_equal(pt[0], long_edge[1])][0]
        angle = math.atan2(opposite_vertex[1] - long_edge_midpoint[1], opposite_vertex[0] - long_edge_midpoint[0])
    elif len(approx) == 4:
        short_edge = sorted_edges[-1]
        short_edge_midpoint = ((short_edge[0][0] + short_edge[1][0]) / 2, (short_edge[0][1] + short_edge[1][1]) / 2)
        angle = math.atan2(short_edge_midpoint[1] - long_edge_midpoint[1], short_edge_midpoint[0] - long_edge_midpoint[0])
    else:
        angle = None
    
    return angle

# レヴィウォークの関数
def generate_levy_step(scale=1.0, alpha=1.5):
    return np.random.pareto(alpha) * scale * 1.5

def fly_motion_levy(center_bfr, theta, pixel_speed, large_step=False):
    step_size = generate_levy_step(scale=pixel_speed * (2 if large_step else 1))
    delta = [step_size * np.cos(theta), step_size * np.sin(theta)]
    center_now = [center_bfr[0] + delta[0], center_bfr[1] + delta[1]]
    theta += np.random.normal(0, np.pi / 12)
    return center_now, theta

def keep_within_arena(position, radius_px):
    dist_from_center = np.sqrt((position[0] - arena_center_x) ** 2 + (position[1] - arena_center_y) ** 2)
    if dist_from_center > radius_px - 5: #アリーナの縁に当たらないように
        angle = math.atan2(position[1] - arena_center_y, position[0] - arena_center_x)
        position[0] = arena_center_x + (radius_px - 5) * np.cos(angle)
        position[1] = arena_center_y + (radius_px - 5) * np.sin(angle)
    return position

# ウィンドウ設定
window_name_arena = 'Arena'
window_name_tracking = 'Tracking'
cv2.namedWindow(window_name_arena, cv2.WINDOW_NORMAL)
cv2.namedWindow(window_name_tracking, cv2.WINDOW_NORMAL)

# サブディスプレイの解像度にウィンドウサイズを設定
sub_display_width = sub_display.width
sub_display_height = sub_display.height
cv2.resizeWindow(window_name_arena, sub_display_width, sub_display_height)

# アリーナの表示位置（サブディスプレイ内の固定位置を指定）
arena_offset_x = (sub_display_width - desired_arena_diameter_px) // 2 - 300
arena_offset_y = (sub_display_height - desired_arena_diameter_px) // 2

# サブディスプレイの左上にウィンドウを移動
cv2.moveWindow(window_name_arena, sub_display.x, sub_display.y)

# メインディスプレイ上のウィンドウ配置
cv2.resizeWindow(window_name_tracking, camera_width, camera_height)
cv2.moveWindow(window_name_tracking, main_display.x + 1000, main_display.y)

# 3コマ前の位置を保持するためのキューを初期化
fly_position_queue = deque([(int(camera_width / 2), int(camera_height / 2))] * 3, maxlen=3)
# 3コマ前の位置を取得
three_frames_ago_position = fly_position_queue[0]

# 初期設定
center_fly_now = [arena_center_x, arena_center_y]
theta = 0
frame_count = 0
end_time = time.time() + duration_seconds
stop_duration_tmp = 0
start_time = time.time()
current_time = start_time
time_passed = 0

# 撮影ループ
while time.time() < end_time:
    current_time_prev = current_time
    current_time = time.time()
    time_interval = current_time - current_time_prev
    time_passed = time_passed  + time_interval
    if Tis.Snap_image(0.001):
        current_frame_tmp = Tis.Get_image()
        current_frame = cv2.resize(current_frame_tmp, (camera_width, camera_height))

        # グレースケール画像の閾値処理
        _, thresh = cv2.threshold(current_frame, 70, 255, cv2.THRESH_BINARY_INV)
        thresh = cv2.GaussianBlur(thresh, (5, 5), 0)
        thresh = cv2.erode(thresh, None, iterations=2)
        thresh = cv2.dilate(thresh, None, iterations=3)
        
        # マスクを設定し、枠線の検出を抑制
        mask = np.zeros((camera_height, camera_width), dtype=np.uint8)
        cv2.circle(mask, (camera_width // 2, camera_height // 2), int(arena_radius_px * 0.69), (255), -1)
        thresh = cv2.bitwise_and(thresh, thresh, mask=mask)

        # 輪郭の検出とクモの位置・向きの計算
        contours, _ = cv2.findContours(thresh, cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE)
        center_spider = (spider_center_x, spider_center_y)
        spider_direction = None

        if contours:
            for i, contour in enumerate(contours):
                area = cv2.contourArea(contour)
                perimeter = cv2.arcLength(contour, True)
                M = cv2.moments(contour)
                if M["m00"] > 0:
                    contour_center = (int(M["m10"] / M["m00"]), int(M["m01"] / M["m00"]))
                    distance_from_arena_center = Calc_distance(contour_center, (camera_width // 2, camera_height // 2))
                    if distance_from_arena_center <= arena_radius_px and 300 < area < 5000 and perimeter > 0:
                        center_spider = contour_center
                        spider_direction = calculate_spider_direction(contour)
                        cv2.drawContours(current_frame, [contour], -1, (255, 0, 0), 2)  # クモの輪郭を描画
                        break

        # サブディスプレイにハエの位置を表示
        arena_frame = np.ones((sub_display_height, sub_display_width, 3), dtype=np.uint8) * 255
        cv2.circle(arena_frame, (arena_center_x, arena_center_y), arena_radius_px, (0, 0, 0), 2)
        fly_position = (int(center_fly_now[0]),
                        int(center_fly_now[1]))
        cv2.ellipse(arena_frame, fly_position, (10, 5), math.degrees(theta), 0, 360, (0, 0, 0), -1)
        cv2.imshow(window_name_arena, arena_frame)

        # ハエの位置をカメラ入力画像におけるハエの位置を取得
        relative_a = (center_fly_now[0] - arena_center_x) * camera_arena_scale * vertical_scale_correction
        relative_b = (center_fly_now[1] - arena_center_y) * camera_arena_scale * horizontal_scale_correction
        fly_position_on_camera = (
            int(camera_width / 2 + relative_b + horizontal_margin),
            int(camera_height / 2 - relative_a - vertical_margin)
        )

        # デバッグ用トラッキングウィンドウにマスクの円を描画
        debug_frame = cv2.cvtColor(current_frame, cv2.COLOR_GRAY2BGR)
        cv2.circle(debug_frame, (camera_width // 2, camera_height // 2), int(arena_radius_px * 0.69), (255, 255, 0), 2)  # マスク範囲を黄色で表示
        cv2.circle(debug_frame, center_spider, 5, (0, 0, 255), -1)  # クモの位置を赤で表示
        cv2.circle(debug_frame, three_frames_ago_position, distance_threshold, (0, 255, 0), 1)  # ハエの枠を緑で表示
        cv2.imshow(window_name_tracking, debug_frame)

        # フレーム画像の保存
        arena_frame_filename = os.path.join(output_dir_arena, f"frame_{frame_count}.png")
        tracking_frame_filename = os.path.join(output_dir_tracking, f"frame_tracking_{frame_count}.png")
        cv2.imwrite(arena_frame_filename, current_frame_tmp) #save original TIS image
        cv2.imwrite(tracking_frame_filename, debug_frame)

        theta_on_camera = theta - (np.pi/2)

        # ハエとクモの距離を計算
        distance = Calc_distance(fly_position_on_camera, center_spider)
        print(f'Fly position: {fly_position_on_camera}, Spider position: {center_spider}')

        # ハエの動きと停止を制御
        stopped = False  # stopped変数を初期化
        if stop_duration_tmp > 0:
            stopped = True
            if current_time - start_time >= stop_duration_tmp:
                stop_duration_tmp = 0
                must_move_until = current_time + 3
                big_move_next = True
                stopped = False
        else:
            if current_time < must_move_until:
                center_fly_bfr = center_fly_now[:]
                center_fly_now, theta = fly_motion_levy(center_fly_bfr, theta, pixel_speed, large_step=big_move_next)
                center_fly_now = keep_within_arena(center_fly_now, arena_radius_px)
                big_move_next = False
            else:
                if stop_duration > 0 and distance <= distance_threshold:
                    stop_duration_tmp = stop_duration
                    start_time = current_time
                    print("The fly meets a spider.")
                    stopped = True
                else:
                    center_fly_bfr = center_fly_now[:]
                    center_fly_now, theta = fly_motion_levy(center_fly_bfr, theta, pixel_speed)
                    center_fly_now = keep_within_arena(center_fly_now, arena_radius_px)
        
        # TSVファイルへの記録
        with open(output_file, mode='a', newline='') as file:
            writer = csv.writer(file, delimiter='\t')
            writer.writerow([
                frame_count, time_interval, time_passed,
                three_frames_ago_position[0], three_frames_ago_position[1],
                center_spider[0], center_spider[1], distance, int(stopped),
                theta_on_camera, spider_direction if spider_direction is not None else "NA"
                ])
        
        # ハエの位置を4フレーム分保存するキューに追加
        fly_position_queue.append(fly_position_on_camera)
        # 3コマ前の位置を取得
        three_frames_ago_position = fly_position_queue[0]

        frame_count += 1
        if cv2.waitKey(1) & 0xFF == 27:
            break

# カメラのパイプラインを停止してウィンドウを閉じる
Tis.Stop_pipeline()
cv2.destroyAllWindows()
