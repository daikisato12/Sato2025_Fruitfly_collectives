#!/usr/bin/env python3
# created by Daiki Sato and updated on 2023/12/8
# python3 5_1_2_simulation_group_socialcontagion.py 200

import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import pandas as pd
import math
import csv
import sys

args = sys.argv

num_agents = 6                   # Number of flies
social_influence =  0.1          # Coefficient of aligning speed with neighboring flies
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
T_freezing_stim = 0.46           # Time to start freezing after stimulus (sec) #0.4 #0.5 #0.46
arena_freezing_factor = 0.05     # Stimulus strength (1+arena_freezing_factor) when a fly is at y = 15
log_param_a = 0.02               # parameter a #0.048 #0.3 #0.04219259 #0.1 #0.2 #0.15 #0.05 #0.05 #0.04 #0.05 #0.048 #0.02 #0.02
log_param_b = 0.92               # parameter b #0.873 #0.2 #0.7744190  #0.5 #0.3 #0.35 #0.80 #0.70 #0.70 #0.85 #0.87  #0.85 #0.90

freezing_rate = 0.9              # freezing rate #0.25 #0.4 #0.6 #0.7 #0.8 #0.75 #0.9 #1.0 #0.9

def minDistAngle(pos, agent_head_directions):
    x, y = pos
    dist_wall = arena_diameter / 2 - (x**2 + y**2)**(1/2)
    angle_wall = (agent_head_directions - math.atan2(y, x) + np.pi) % (2 * np.pi) - np.pi
    return dist_wall, angle_wall

def generate_trajectory(num_steps):
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
                        recovery_factor = log_param_a * math.log((t * dt - T_starting_stim - T_freezing_stim) % T_interval_stim) + log_param_b
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
            print(iter, i, t, (t * dt - T_starting_stim), (t * dt - T_starting_stim) % T_interval_stim, agent_speeds[t][i], motion_cues_diff)

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

f = lambda x: x ** 3

for iter in range(int(args[1])):
    agent_positions_res = generate_trajectory(num_steps)

    nd = [np.insert(agent_positions_res[0][i3], 0, [range(num_agents)], axis = 1) for i3 in range(len(agent_positions_res[0]))]
    df = pd.DataFrame([list(l) for l in nd]).stack().apply(pd.Series).reset_index(1, drop=True)
    df.reset_index(inplace=True)
    df.columns = ['frame', 'id', 'pos_x', 'pos_y']
    df.to_csv(f'../data/group_socialcontagion/df_pos_{iter}.tsv', sep='\t', index=False)
