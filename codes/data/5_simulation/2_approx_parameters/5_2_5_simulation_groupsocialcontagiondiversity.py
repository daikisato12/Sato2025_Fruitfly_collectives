#!/usr/bin/env python3
# created by Daiki Sato and updated on 2025/3/20
# python3 5_2_5_simulation_groupsocialcontagiondiversity.py 10 0.1

import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import pandas as pd
import math
import sys
import os

args = sys.argv

num_agents = 6                   # Number of flies
agent_size = 1                   # agent size (mm)
arena_diameter = 30              # Arena diameter (mm)
perimeter_dist = 1               # Perimeter region distance to walls (mm)
decel_rate_peri = 0.5;           # agent_velocities reduction factor when located in the perimeter
average_speed = 3.8              # average speed (mm/sec)
power = 2.0                      # power law parameter
speed_influence = 0.6            # influence of speed at t-1
mu = 0.0                         # Rotation agent_velocities Gaussian distribution mean (rad/sec)
sigma = (180 / 360) * 2 * np.pi  # Rotation agent_velocities Gaussian distribution standard deviation (rad/sec)
dt = 0.02                        # Simulation-step time increment (sec)
T = 600                          # Simulation time (sec)

T_starting_stim = 300            # Starting time of stimulus (sec)
T_interval_stim = 15             # Interval of stimulus (sec)
T_freezing_stim = 0.46           # Time to start freezing after stimulus (sec)
arena_freezing_factor = 0.05     # Stimulus strength (1+arena_freezing_factor) when a fly is at y = 15
log_param_a = [0.005, 0.006, 0.007, 0.008, 0.009, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1]  
log_param_b = 1.02 - 2.45 * np.array(log_param_a)

freezing_rate = 0.9              # freezing rate

social_influence = float(args[3]) # Coefficient of aligning speed with neighboring flies

def minDistAngle(pos, agent_head_directions):
    x, y = pos
    dist_wall = arena_diameter / 2 - (x**2 + y**2)**(1/2)
    angle_wall = (agent_head_directions - math.atan2(y, x) + np.pi) % (2 * np.pi) - np.pi
    return dist_wall, angle_wall

def generate_trajectory(num_steps, param1_a, param1_b, param2_a, param2_b):
    # store arrays
    agent_positions = [np.zeros((num_agents, 2)) for _ in range(num_steps)]
    agent_velocities = [np.zeros((num_agents, 2)) for _ in range(num_steps)]
    agent_head_directions = [np.zeros(num_agents) for _ in range(num_steps)] # head direction
    agent_speeds = [np.random.pareto(power, num_agents) * average_speed for _ in range(num_steps)]
    motion_cues = [np.zeros(num_agents) for _ in range(num_steps)]
    random_turn = [np.random.normal(mu, sigma, num_agents) * dt for _ in range(num_steps)]
    # initial values
    agent_head_directions[0] = np.random.rand(num_agents) * 2 * np.pi
    agent_positions[0] = np.random.rand(num_agents, 2) * arena_diameter - arena_diameter / 2 
    # iteration of trajectory
    for t in range(num_steps):
        for i in range(num_agents):
            if t == 0:
                while agent_positions[t][i][0] ** 2 + agent_positions[t][i][1] ** 2 > (arena_diameter / 2) ** 2:
                    agent_positions[t][i] = np.random.rand(2) * arena_diameter - arena_diameter / 2

            # avoid wall
            turn_angle = random_turn[t][i]
            dist_wall, angle_wall = minDistAngle(agent_positions[t][i], agent_head_directions[t][i])
            if (dist_wall < perimeter_dist) & (abs(angle_wall) < np.pi / 2):
                agent_speeds[t][i] *= decel_rate_peri # deceleration
                sigma_wall = (10 / 360) * 2 * np.pi
                random_turn_tmp = np.random.normal(mu, sigma_wall) * dt
                turn_angle = random_turn_tmp + np.sign(angle_wall) * (np.pi / 2 - abs(angle_wall))

            # current speed is based on that of previous time step
            if agent_speeds[t][i] != 0:
                agent_speeds[t][i] = (1 - speed_influence) * agent_speeds[t][i] + speed_influence * agent_speeds[t-1][i]

            # freezing and recovery
            if round(t * dt - T_starting_stim, 2) >= T_freezing_stim:
                if round((t * dt - T_starting_stim) % T_interval_stim, 2) == T_freezing_stim:
                    arena_factor = (1 + arena_freezing_factor * (agent_positions[t][i][1] + arena_diameter / 2) / arena_diameter)
                    if np.random.uniform(0, 1) < freezing_rate * arena_factor:
                        agent_speeds[t][i] = 0
                else:
                    t_recent_stim = int((int((t * dt - T_starting_stim) / T_interval_stim) * T_interval_stim + T_starting_stim) / dt + T_freezing_stim / dt) 
                    if agent_speeds[t_recent_stim][i] == 0:
                        if i < 3:
                            recovery_factor = param1_a * math.log((t * dt - T_starting_stim - T_freezing_stim) % T_interval_stim) + param1_b
                        else:
                            recovery_factor = param2_a * math.log((t * dt - T_starting_stim - T_freezing_stim) % T_interval_stim) + param2_b
                        if np.random.uniform(0, 1) > recovery_factor:
                            agent_speeds[t][i] = 0

            # adjust speed with other flies
            motion_cues_diff = motion_cues[t-1][i] - motion_cues[t-2][i]
            agent_speeds[t][i] = (1 - social_influence) * agent_speeds[t][i] + social_influence * motion_cues_diff

            if agent_speeds[t][i] < 0:
                agent_speeds[t][i] = 0

            # check if agent walks through wall
            agent_velocities[t][i] = np.multiply([math.cos(agent_head_directions[t][i]), math.sin(agent_head_directions[t][i])], agent_speeds[t][i])
            agent_positions_tmp = agent_positions[t][i] + agent_velocities[t][i] * dt
            while agent_positions_tmp[0] ** 2 + agent_positions_tmp[1] ** 2 > (arena_diameter / 2) ** 2:
                agent_speeds[t][i] *= decel_rate_peri # deceleration
                agent_velocities[t][i] = np.multiply([math.cos(agent_head_directions[t][i]), math.sin(agent_head_directions[t][i])], agent_speeds[t][i])
                agent_positions_tmp = agent_positions[t][i] + agent_velocities[t][i] * dt

            distances = np.linalg.norm(agent_positions[t] - agent_positions[t][i], axis=1)
            distances[i] = np.inf
            a_li = np.divide(agent_size, np.multiply(2, distances)).tolist()
            visual_angles = np.multiply(list(map(math.atan, a_li)),2)
            motion_cues[t][i] = np.sum(np.multiply(visual_angles, agent_speeds[t]))

            agent_velocities[t][i] = np.multiply([math.cos(agent_head_directions[t][i]), math.sin(agent_head_directions[t][i])], agent_speeds[t][i])
            if t < num_steps-1:
                agent_head_directions[t+1][i] = (agent_head_directions[t][i] + turn_angle) % (2 * np.pi) # turn, 
                agent_positions[t+1][i] = agent_positions[t][i] + agent_velocities[t][i] * dt
            print(strain1_index, strain2_index, iter, i, t, (t * dt - T_starting_stim), (t * dt - T_starting_stim) % T_interval_stim, agent_speeds[t][i], motion_cues_diff)

    return agent_positions, agent_velocities, agent_speeds, agent_head_directions

def init():
    return p1,p2,p3,p4,p5 

def update(i):
    p1.set_data((agent_positionss[0][i+4,0], agent_positionss[0][i+4,1]))
    p1.set_alpha(1)
    p2.set_data((agent_positionss[0][i+3,0], agent_positionss[0][i+3,1]))
    p2.set_alpha(0.8)
    p3.set_data((agent_positionss[0][i+2,0], agent_positionss[0][i+2,1]))
    p3.set_alpha(0.6)
    p4.set_data((agent_positionss[0][i+1,0], agent_positionss[0][i+1,1]))
    p4.set_alpha(0.4)
    p5.set_data((agent_positionss[0][i,0], agent_positionss[0][i,1]))
    p5.set_alpha(0.2)
    return p1,p2,p3,p4,p5 


num_steps = round(T / dt)
num_index = len(log_param_a) #* len(log_param_b) 
f = lambda x: x ** 3
repn = int(args[1])
outdir = f"../data/5_simulation/2_approx_parameters/raw/groupsocialcontagiondiversity"
os.makedirs(outdir, exist_ok=True)
for strain1_index in range(num_index):
    param1_a = log_param_a[strain1_index] #int(strain1_index/len(log_param_a))
    param1_b = log_param_b[strain1_index] #int(strain1_index%len(log_param_b))
    for strain2_index in range(strain1_index, num_index):
        param2_a = log_param_a[strain2_index] #int(strain2_index/len(log_param_a))
        param2_b = log_param_b[strain2_index] #int(strain2_index%len(log_param_b))
        for iter in range(repn):
            agent_positions_res = generate_trajectory(num_steps, param1_a, param1_b, param2_a, param2_b)
            nd = [np.insert(agent_positions_res[0][i3], 0, [range(num_agents)], axis = 1) for i3 in range(len(agent_positions_res[0]))]
            df = pd.DataFrame([list(l) for l in nd]).stack().apply(pd.Series).reset_index(1, drop=True)
            df.reset_index(inplace=True)
            df.columns = ['frame', 'id', 'pos_x', 'pos_y']
            df['strain1'] = f'strain{strain1_index}'
            df['strain2'] = f'strain{strain2_index}'
            df['rep'] = iter
            df.to_csv(f'{outdir}/df_pos_{strain1_index}_{strain2_index}_{iter}.tsv', sep='\t', index=False)
