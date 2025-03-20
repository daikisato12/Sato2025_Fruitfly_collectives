#!/bin/bash

for t in 0.0 0.5 1.0 1.5 2.0 2.5 3.0 3.5 4.0 4.5; do
    for a in 0.01 0.02 0.03 0.04 0.05 0.1 0.2 0.3 0.4 0.5; do
        for social_factor in 0.0 0.01 0.02 0.03 0.04 0.05 0.1 0.2 0.3 0.4 0.5; do
            echo "$t\t$a\t$social_factor" >> "./grid_simulation_group_parameters.txt"
        done
    done
done