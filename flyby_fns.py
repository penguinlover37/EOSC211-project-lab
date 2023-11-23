import numpy as np
import math
import matplotlib.pyplot as plt

def grav_acc(s_x, s_y, planet_mass):
    '''
    takes in x,y components of position vector from the planet and mass of planet
    finds acceleration components

    input: s_x, s_y, planet_mass
    output: a_x (m/s), a_y (m/s)
    '''
    s = math.sqrt(s_x**2 + s_y**2)
    #formula for a is GMp/s^2
    big_G = 6.67e-11
    a = (big_G*planet_mass)/(s**2)

    sin_beta = -(s_x)/s
    cos_beta = -(s_y)/s

    a_x = a * sin_beta
    a_y = a * cos_beta

    return(a_x, a_y)

def checkinit(s_x0, s_y0, v_x0, v_y0, planet_radius):
    '''
    Checks starting position, s, determined by inputs s_x0 and s_y0 against planet_radius to prevent impossible values.
    Checks v_y0 to ensure y-axis velocity is positive. 
    '''
    s = math.sqrt(s_x0**2 + s_y0**2)
    if np.abs(s) <= planet_radius:
        raise ValueError('Please enter coordinates that are above the planetary surface i.e greater than the radius.')
    elif v_y0 < 0:
        raise ValueError('Please enter a positive velocity, in the direction towards the planet.')
    else:
        pass
    
def sc_vel_pos_change(a_x, a_y, v_x, v_y, time_step): 
    '''
    Takes in instantaneous acceleration and velocity x&y components at a certain timestep
    Returns the delta(change) of x-y position vectors and velocity vectors to track change over time.
    input: a_x, a_y, v_x, v_y, time_step
    output: ds_x, ds_y, dv_x, dv_y
    '''
    dv_x = a_x * time_step
    dv_y = a_y * time_step
    
    ds_x = v_x * time_step + 0.5*a_x*(time_step)**2
    ds_y = v_y * time_step + 0.5*a_y*(time_step)**2
    
    return ds_x, ds_y, dv_x, dv_y

def get_traj(s_x0, s_y0, v_x0, v_y0, time_step, total_time, planet_mass, planet_radius):
    '''
    defines four arrays with the size of every time interval determined by total_time and time_step: time, acc, vel, pos
    matches the instantaneous acceleration, velocity, and position determined by initial position and velocity inputs with each time interval
    returns data organized into x & y components, resulting in tracking of vector arrays over time.
    
    inputs: s_x0, s_y0, v_x0, v_y0, time_step, total_time, planet_mass, planet_radius
    output: arrays time, acc, vel, pos
    '''
    total_steps = int(total_time/time_step) + 1
    
    time = np.linspace(0, total_time, total_steps)
    
    acc = np.ones((time.size, 2))*np.nan
    vel = np.ones((time.size, 2))*np.nan
    pos = np.ones((time.size, 2))*np.nan
    
    checkinit(s_x0, s_y0, v_x0, v_y0, planet_radius)
        
    pos[0, 0] = s_x0
    pos[0, 1] = s_y0
    vel[0, 0] = v_x0
    vel[0, 1] = v_y0
    acc[0, 0], acc[0, 1] = grav_acc(pos[0, 0], pos[0, 1], planet_mass)

    for i in range(1, len(time)):

        
        ds_x, ds_y, dv_x, dv_y = sc_vel_pos_change(acc[i-1, 0], acc[i-1, 1], vel[i-1, 0], vel[i-1, 1], time_step)
        
        vel[i, 0] = vel[i-1, 0] + dv_x
        vel[i, 1] = vel[i-1, 1] + dv_y
        pos[i, 0] = pos[i-1, 0] + ds_x
        pos[i, 1] = pos[i-1, 1] + ds_y
        
        acc[i, 0], acc[i, 1] = grav_acc(pos[i, 0], pos[i, 1], planet_mass)

    
    return time, acc, vel, pos