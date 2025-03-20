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
#python3 6_1_detect_spider_move_6flies_each_duration_interaction.py ./test 360 6 50 24 m3 male 1 6 1,1,1,3,3,3 0.2


# Command-line argument processing
if len(sys.argv) < 9:
    print("Usage: python 6_1_detect_spider_move_6flies_each_duration_interaction.py <output_directory> <duration_seconds> <stop_duration> <distance_threshold> <fly_speed> <individual_id> <sex> <trial_id> [<num_flies>] [<stop_durations_comma_separated>] [<influence_factor>]")
    sys.exit(1)

# Extract parameters from arguments
output_base_dir = sys.argv[1]
duration_seconds = int(sys.argv[2])
default_stop_duration = float(sys.argv[3])
distance_threshold = int(sys.argv[4])
average_speed_mm = float(sys.argv[5])
individual_id = sys.argv[6]
sex = sys.argv[7]
trial_id = sys.argv[8]

# Optional: Number of flies (default to 1 if not provided)
num_flies = int(sys.argv[9]) if len(sys.argv) >= 10 else 1
if num_flies < 1 or num_flies > 6:
    print("Number of flies must be between 1 and 6.")
    sys.exit(1)

# Optional: Per-fly stop durations (comma-separated)
per_fly_stop_durations = []
if len(sys.argv) >= 11:
    # sys.argv[10] should contain comma-separated stop durations
    stop_durations_str = sys.argv[10]
    per_fly_stop_durations = [float(s) for s in stop_durations_str.split(',')]
    if len(per_fly_stop_durations) != num_flies:
        print(f"Number of stop durations ({len(per_fly_stop_durations)}) does not match number of flies ({num_flies}).")
        sys.exit(1)
else:
    # Use default stop_duration for all flies
    per_fly_stop_durations = [default_stop_duration] * num_flies

# Optional: Influence factor (default to 0.1 if not provided)
influence_factor = 0
if len(sys.argv) >= 12:
    try:
        influence_factor = float(sys.argv[11])
    except ValueError:
        print("Influence factor must be a float.")
        sys.exit(1)

# Create output directories with argument details
start_time_str = datetime.now().strftime("%Y%m%d%H%M%S")
stop_durations_formatted = sys.argv[10].replace(",", "-") if len(sys.argv) >= 11 else f"{default_stop_duration}"
arena_dir_name = f"images/{start_time_str}_{duration_seconds}s_{stop_durations_formatted}s_{distance_threshold}px_{average_speed_mm}mmps_{influence_factor}inf_{individual_id}_{sex}_{trial_id}/arena"
tracking_dir_name = f"images/{start_time_str}_{duration_seconds}s_{stop_durations_formatted}s_{distance_threshold}px_{average_speed_mm}mmps_{influence_factor}inf_{individual_id}_{sex}_{trial_id}/tracking"
output_dir_arena = os.path.join(output_base_dir, arena_dir_name)
output_dir_tracking = os.path.join(output_base_dir, tracking_dir_name)
os.makedirs(output_dir_arena, exist_ok=True)
os.makedirs(output_dir_tracking, exist_ok=True)

# Set TSV file path
output_file = os.path.join(output_base_dir, f"tracks/{start_time_str}_{duration_seconds}s_{stop_durations_formatted}s_{distance_threshold}px_{average_speed_mm}mmps_{influence_factor}inf_{individual_id}_{sex}_{trial_id}_track.tsv")
os.makedirs(os.path.dirname(output_file), exist_ok=True)  # Ensure the tracks directory exists

# Display and camera settings
monitors = get_monitors()
main_display = monitors[0]
sub_display = monitors[1] if len(monitors) > 1 else monitors[0]
camera_width, camera_height = 640, 480
arena_diameter_mm = 84
desired_arena_diameter_px = 600
scale_px_per_mm = desired_arena_diameter_px / arena_diameter_mm
arena_radius_px = desired_arena_diameter_px // 2 - 10
camera_arena_scale = camera_height / desired_arena_diameter_px
horizontal_margin = 0
vertical_margin = 0
fps = 60
pixel_speed = (average_speed_mm * scale_px_per_mm) / fps  # Base pixel speed
spider_center_x, spider_center_y = 500, 0  # Adjust as needed

# Arena center position in the window (fixed)
arena_center_x = 617
arena_center_y = 556

# Scale correction factors
vertical_scale_correction = 0.9  # 1920:1080 640:480
horizontal_scale_correction = 0.9  # 1920:1080 640:480


# Define the Fly class to manage individual fly states
class Fly:
    def __init__(self, fly_id, initial_position, theta=0, stop_duration=0, influence_factor=0.1, pixel_speed = pixel_speed):
        self.fly_id = fly_id
        self.position = initial_position[:]  # Current position [x, y] in arena coordinates
        self.theta = theta  # Current direction in radians
        self.must_move_until = 0  # Timestamp until which the fly must move
        self.stopped = False  # Whether the fly is currently stopped
        self.big_move_next = False  # Flag for a large movement in the next step
        self.fly_position_queue = deque([self.calculate_camera_position(initial_position)] * 3, maxlen=3)  # History of camera positions
        self.three_frames_ago_position = self.fly_position_queue[0]  # Position three frames ago in camera coordinates
        self.stop_duration_tmp = 0  # Temporary stop duration
        self.stop_start_time = 0  # Time when the fly started to stop
        self.stop_duration = stop_duration  # Individual stop duration
        self.influence_factor = influence_factor  # Influence factor for motion cue
        self.fly_position_on_camera = self.calculate_camera_position(initial_position)  # Current position in camera coordinates
        self.motion_cue = 0  # Motion cue based on other flies
        self.prev_motion_cue = 0
        self.prev2_motion_cue = 0
        self.prev_speed = pixel_speed

    # Method to calculate fly_position_on_camera based on current arena position
    def calculate_camera_position(self, current_position, arena_center_x=617, arena_center_y=556,
                                  camera_arena_scale=1.0, vertical_scale_correction=0.9,
                                  horizontal_scale_correction=0.9, camera_width=640, camera_height=480,
                                  horizontal_margin=0, vertical_margin=0):
        relative_a = (current_position[0] - arena_center_x) * camera_arena_scale * vertical_scale_correction
        relative_b = (current_position[1] - arena_center_y) * camera_arena_scale * horizontal_scale_correction
        fly_position_on_camera = (
            int(camera_width / 2 + relative_b + horizontal_margin),
            int(camera_height / 2 - relative_a - vertical_margin)
        )
        return fly_position_on_camera

    # Method to update camera position
    def update_camera_position(self, arena_center_x, arena_center_y, camera_arena_scale,
                               vertical_scale_correction, horizontal_scale_correction,
                               camera_width, camera_height, horizontal_margin, vertical_margin):
        self.fly_position_on_camera = self.calculate_camera_position(
            self.position,
            arena_center_x,
            arena_center_y,
            camera_arena_scale,
            vertical_scale_correction,
            horizontal_scale_correction,
            camera_width,
            camera_height,
            horizontal_margin,
            vertical_margin
        )
        self.prev_speed = self.motion_speed
        self.fly_position_queue.append(self.fly_position_on_camera)
        self.three_frames_ago_position = self.fly_position_queue[0]

    # Method to calculate motion cue based on other flies
    def calculate_motion_cue(self, other_flies):
        total_cue = 0
        for other_fly in other_flies:
            if other_fly.fly_id == self.fly_id:
                continue  # Skip self
            # Calculate relative position
            dx = other_fly.position[0] - self.position[0]
            dy = other_fly.position[1] - self.position[1]
            distance = math.hypot(dx, dy)
            if distance == 0:
                continue  # Avoid division by zero
            # Influence decreases with distance
            # dividing by distance approximates multiplying (2*math.atan(1/(2*distance)))
            influence = other_fly.motion_speed / (distance + 1e-5)
            total_cue += influence
        self.prev2_motion_cue = self.prev_motion_cue
        self.prev_motion_cue = self.motion_cue
        self.motion_cue = total_cue

    @property
    def motion_speed(self):
        # Base speed plus motion cue adjusted by influence factor
        # You can adjust the influence factor as needed
        delta_cue = self.prev_motion_cue - self.prev2_motion_cue
        return self.prev_speed + delta_cue * self.influence_factor

# Function to calculate distance between two points
def Calc_distance(center1, center2):
    return np.sqrt((center2[0] - center1[0]) ** 2 + (center2[1] - center1[1]) ** 2)

# Function to calculate spider's direction based on contour
def calculate_spider_direction(contour):
    approx = cv2.approxPolyDP(contour, 0.04 * cv2.arcLength(contour, True), True)
    edges = [(approx[i][0], approx[j][0]) for i in range(len(approx)) for j in range(i+1, len(approx))]
    sorted_edges = sorted(edges, key=lambda pair: Calc_distance(pair[0], pair[1]), reverse=True)
    
    if not sorted_edges:
        return None
    
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

# Function to generate a Levy flight step
def generate_levy_step(scale=1.0, alpha=1.5):
    return np.random.pareto(alpha) * scale * 1.5

# Function to update fly's position based on Levy flight
def fly_motion_levy(center_bfr, theta, pixel_speed, large_step=False):
    step_size = generate_levy_step(scale=pixel_speed * (2 if large_step else 1))
    delta = [step_size * np.cos(theta), step_size * np.sin(theta)]
    center_now = [center_bfr[0] + delta[0], center_bfr[1] + delta[1]]
    theta += np.random.normal(0, np.pi / 12)
    return center_now, theta

# Function to ensure fly stays within the arena
def keep_within_arena(position, radius_px, arena_center_x, arena_center_y):
    dist_from_center = np.sqrt((position[0] - arena_center_x) ** 2 + (position[1] - arena_center_y) ** 2)
    if dist_from_center > radius_px - 5:  # Prevent touching the arena edge
        angle = math.atan2(position[1] - arena_center_y, position[0] - arena_center_x)
        position[0] = arena_center_x + (radius_px - 5) * math.cos(angle)
        position[1] = arena_center_y + (radius_px - 5) * math.sin(angle)
    return position

# Initialize TIS camera
Tis = TIS.TIS()
Tis.openDevice("28020355", camera_width, camera_height, "60/1", TIS.SinkFormats.GRAY8, True)
Tis.Start_pipeline()
Tis.Set_Property("Gain Auto", True)
Tis.Set_Property("Exposure Auto", False)
Tis.Set_Property("Exposure Time (us)", 100)

# Capture background
if Tis.Snap_image(0.1):
    background = Tis.Get_image()
    background = cv2.resize(background, (camera_width, camera_height))
else:
    print("Failed to capture background image.")
    background = np.zeros((camera_height, camera_width), dtype=np.uint8)

# Initialize TSV file with header including Fly_ID
with open(output_file, mode='w', newline='') as file:
    writer = csv.writer(file, delimiter='\t')
    # Header with multiple flies
    header = ["Frame", "Fly_ID", "Time_interval", "Time", "Fly_X", "Fly_Y",
              "Spider_X", "Spider_Y", "Distance", "Stopped", "Fly_Direction_Radians",
              "Spider_Direction_Radians"]
    writer.writerow(header)

# Initialize window settings
window_name_arena = 'Arena'
window_name_tracking = 'Tracking'
cv2.namedWindow(window_name_arena, cv2.WINDOW_NORMAL)
cv2.namedWindow(window_name_tracking, cv2.WINDOW_NORMAL)

# Set window sizes based on sub-display
sub_display_width = sub_display.width
sub_display_height = sub_display.height
cv2.resizeWindow(window_name_arena, sub_display_width, sub_display_height)

# Move arena window to sub-display
cv2.moveWindow(window_name_arena, sub_display.x, sub_display.y)

# Set tracking window on main display
cv2.resizeWindow(window_name_tracking, camera_width, camera_height)
cv2.moveWindow(window_name_tracking, main_display.x + 1000, main_display.y)

# Initialize multiple flies with random positions within the arena
flies = []
for i in range(num_flies):
    angle = random.uniform(0, 2 * math.pi)
    radius = random.uniform(0, arena_radius_px - 10)
    x = arena_center_x + int(radius * math.cos(angle))
    y = arena_center_y + int(radius * math.sin(angle))
    initial_theta = random.uniform(0, 2 * math.pi)
    fly = Fly(
        fly_id=i,
        initial_position=[x, y],
        theta=initial_theta,
        stop_duration=per_fly_stop_durations[i],
        influence_factor=influence_factor  # Pass influence_factor to each Fly instance
    )
    # Initialize fly_position_queue with initial camera position
    fly.fly_position_on_camera = fly.calculate_camera_position(
        fly.position,
        arena_center_x,
        arena_center_y,
        camera_arena_scale,
        vertical_scale_correction,
        horizontal_scale_correction,
        camera_width,
        camera_height,
        horizontal_margin,
        vertical_margin
    )
    fly.fly_position_queue = deque([fly.fly_position_on_camera] * 3, maxlen=3)
    fly.three_frames_ago_position = fly.fly_position_queue[0]
    flies.append(fly)

frame_count = 0
end_time = time.time() + duration_seconds
start_time_global = time.time()
current_time = start_time_global
time_passed = 0

# Main capture and processing loop
while time.time() < end_time:
    current_time_prev = current_time
    current_time = time.time()
    time_interval = current_time - current_time_prev
    time_passed += time_interval

    if Tis.Snap_image(0.001):
        current_frame_tmp = Tis.Get_image()
        current_frame = cv2.resize(current_frame_tmp, (camera_width, camera_height))

        # Thresholding and preprocessing
        _, thresh = cv2.threshold(current_frame, 70, 255, cv2.THRESH_BINARY_INV)
        thresh = cv2.GaussianBlur(thresh, (5, 5), 0)
        thresh = cv2.erode(thresh, None, iterations=2)
        thresh = cv2.dilate(thresh, None, iterations=3)
        
        # Apply mask to suppress frame edges
        mask = np.zeros((camera_height, camera_width), dtype=np.uint8)
        cv2.circle(mask, (camera_width // 2, camera_height // 2), int(arena_radius_px * 0.69), (255), -1)
        thresh = cv2.bitwise_and(thresh, thresh, mask=mask)

        # Detect contours to find the spider
        contours, _ = cv2.findContours(thresh, cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE)
        center_spider = (spider_center_x, spider_center_y)
        spider_direction = None

        if contours:
            for contour in contours:
                area = cv2.contourArea(contour)
                perimeter = cv2.arcLength(contour, True)
                M = cv2.moments(contour)
                if M["m00"] > 0:
                    contour_center = (int(M["m10"] / M["m00"]), int(M["m01"] / M["m00"]))
                    distance_from_arena_center = Calc_distance(contour_center, (camera_width // 2, camera_height // 2))
                    # 修正ポイント: areaの下限を引き上げ
                    if distance_from_arena_center <= arena_radius_px and 1000 < area < 5000 and perimeter > 0:
                        center_spider = contour_center
                        spider_direction = calculate_spider_direction(contour)
                        cv2.drawContours(current_frame, [contour], -1, (255, 0, 0), 2)  # Draw spider contour
                        break

        # Create a blank arena frame and draw the arena boundary
        arena_frame = np.ones((sub_display_height, sub_display_width, 3), dtype=np.uint8) * 255
        cv2.circle(arena_frame, (arena_center_x, arena_center_y), arena_radius_px, (0, 0, 0), 2)

        # Calculate motion cues for all flies
        for fly in flies:
            fly.calculate_motion_cue(flies)

        # Update and draw each fly
        for fly in flies:
            # Handle stopping logic
            if fly.stop_duration_tmp > 0:
                fly.stopped = True
                if current_time - fly.stop_start_time >= fly.stop_duration_tmp:
                    fly.stop_duration_tmp = 0
                    fly.must_move_until = current_time + 3
                    fly.big_move_next = True
                    fly.stopped = False
            else:
                if current_time < fly.must_move_until:
                    # Perform movement with possible large_step
                    center_bfr = fly.position[:]
                    fly.position, fly.theta = fly_motion_levy(
                        center_bfr, fly.theta, fly.motion_speed, large_step=fly.big_move_next
                    )
                    fly.position = keep_within_arena(fly.position, arena_radius_px, arena_center_x, arena_center_y)
                    fly.big_move_next = False
                else:
                    # Calculate distance using arena coordinates
                    distance_to_spider = Calc_distance(fly.fly_position_on_camera, center_spider)
                    
                    if fly.stop_duration > 0 and distance_to_spider <= distance_threshold:
                        fly.stop_duration_tmp = fly.stop_duration
                        fly.stop_start_time = current_time
                        fly.stopped = True
                        print(f"Fly {fly.fly_id} meets a spider.")
                    else:
                        # Perform normal movement
                        center_bfr = fly.position[:]
                        fly.position, fly.theta = fly_motion_levy(center_bfr, fly.theta, fly.motion_speed)
                        fly.position = keep_within_arena(fly.position, arena_radius_px, arena_center_x, arena_center_y)

            # Update camera position and position queue
            fly.update_camera_position(
                arena_center_x, arena_center_y, camera_arena_scale,
                vertical_scale_correction, horizontal_scale_correction,
                camera_width, camera_height, horizontal_margin, vertical_margin
            )

            # Calculate distance to spider using camera coordinates
            distance = Calc_distance(fly.fly_position_on_camera, center_spider)

            # Log data for the fly
            with open(output_file, mode='a', newline='') as file:
                writer = csv.writer(file, delimiter='\t')
                writer.writerow([
                    frame_count,
                    fly.fly_id,
                    round(time_interval, 6),
                    round(time_passed, 6),
                    fly.three_frames_ago_position[0],
                    fly.three_frames_ago_position[1],
                    center_spider[0],
                    center_spider[1],
                    round(distance, 2),
                    int(fly.stopped),
                    round(fly.theta - (np.pi / 2), 6),
                    round(spider_direction, 6) if spider_direction is not None else "NA"
                ])
            
            # Draw the fly on the arena frame (arena coordinates)
            fly_position = tuple(map(int, fly.position))
            cv2.ellipse(arena_frame, fly_position, (10, 5), math.degrees(fly.theta), 0, 360, (0, 0, 0), -1)
        
        # Display the updated arena with all flies
        cv2.imshow(window_name_arena, arena_frame)

        # Debug/tracking window visuals
        debug_frame = cv2.cvtColor(current_frame, cv2.COLOR_GRAY2BGR)
        cv2.circle(debug_frame, (camera_width // 2, camera_height // 2), int(arena_radius_px * 0.69), (255, 255, 0), 2)  # Mask boundary
        cv2.circle(debug_frame, center_spider, 5, (0, 0, 255), -1)  # Spider position
        for fly in flies:
            # Use fly_position_on_camera from the queue for drawing
            fly_third_pos = tuple(map(int, fly.three_frames_ago_position))
            cv2.circle(debug_frame, fly_third_pos, distance_threshold, (0, 255, 0), 1)  # Green frame around fly's past position

        # Display the debug/tracking frame
        cv2.imshow(window_name_tracking, debug_frame)

        # Save frames to output directories
        arena_frame_filename = os.path.join(output_dir_arena, f"frame_{frame_count}.png")
        tracking_frame_filename = os.path.join(output_dir_tracking, f"frame_tracking_{frame_count}.png")
        cv2.imwrite(arena_frame_filename, arena_frame)  # Save arena_frame with all flies
        cv2.imwrite(tracking_frame_filename, debug_frame)  # Save debug_frame

        frame_count += 1

        # Exit on pressing 'Esc'
        if cv2.waitKey(1) & 0xFF == 27:
            break

# Stop the camera pipeline and close all windows
Tis.Stop_pipeline()
cv2.destroyAllWindows()
