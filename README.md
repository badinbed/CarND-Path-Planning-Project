[//]: # (Image References)


# Path Planning Project

## Introduction
The goal of this project is to safely navigate around a virtual highway with other traffic that is driving +-10 MPH of the 50 MPH speed limit. Sensor fusion data is provided by the simulator as well as a list of sparse waypoints around the highway. The car should try to go as close as possible to the 50 MPH speed limit, which means passing slower traffic when possible, note that other cars will try to change lanes too. The car should avoid hitting other cars at all cost as well as driving inside of the marked road lanes at all times, unless going from one lane to another. The car should be able to make one complete loop around the 6946m highway. Since the car is trying to go 50 MPH, it should take a little over 5 minutes to complete 1 loop. Also the car should not experience total acceleration over 10 m/s^2 and jerk that is greater than 10 m/s^3. 

## Algorithm
The basic planner consist of 4 parts. At first the data from the simulator is parsed. The information from sensor fusion is converted into vehicles organized into lanes. This allows to determine lane speeds, free spaces on neighboring lanes and to find close vehicles to avoid collisions.

The next step evaluates cost functions which will determine lane changes. This step is only evaluated once every second because it can be expensive and is only involved in safety measures indirectly. There are two cost functions. 

The first one evaluates the efficiency of each lane. This is determined by the lowest speed of another vehicle traveling ahead of the ego vehicle. I minor constant cost is added for other lanes than the host lane to avoid constant lane switches for minor gains. Afterwards a cost is added proportional to the distance to the optional lane regarding the speed value. This allows the vehicle to make intermediate lane switches to slower lanes than the host lane to reach a fast lane further away.

The second cost function evaluates the safety for a potential lane switch. It is binary and either return 0 or 1. The cost for a lane switch to the host lane is always 0. The cost for a lane change to a distant lane which would require a lane jump is always 1. The cost for a lane change to a neighboring lane is 0 of there is enough freespace at the projected location, 1 otherwise.

The two cost functions are weighted . The weight for the safety cost is higher than the efficiency cost. If another lane has a lower cost than the host lane, a lane change will be initialized.

In a third step a cubic spline is created that will later help generating the trajectory. The spline consist of 5 points. The first two points are the end points of the previous trajectory. This allows for a smooth transition. The last three points are generated from the waypoints using a distance of 30 meters for each point.

In a last step the trajectory is extended with new points. The vehicle's velocity is determined by teh target velocity and potential close vehicles to avoid collisions. The target lane for lane switches is considered as well.

## Rubic Points

### The car is able to drive at least 4.32 miles without incident.

This requirement basically states, that the next three conditions will be held for a full round of the map. The actual challenge though is to drive without incident more than 4.32 miles. The frenet coordinates reset to 0 at the end of the loop. This affects the sensor fusion when determining whether another vehicle is too close. 
This problem is solved by the Vehicle class and its function ```	
double Vehicle::getRelativeS(double rs, secs_look_ahead) const
```b.
It returns the relative delta S to the given reference value rs while also projecting the vehicle's current S-coori´dinate to the time of rs. It safely wraps around the reset of the waypoints.

### The car drives according to the speed limit.

The car's velocity is bound by the variable target_val which is set to 49.5 mph. When the trajectory is extended by new waypoints it is made sure, that the euclidean distance traveled from one point to another does not exceed 0.442 meter. The car travels to the next point of the trajectory queue with 50 fps. That means 0.442 meters per frame is roughly equals to 21 m/s which is equals to 49.5 mph.

### Max Acceleration and Jerk are not Exceeded.

Restricting acceleration works similar to restricting speed. For each point of the trajectory it is guaranteed during its generation, that the speed will neither increase nor decrease by more than 0.14 meter per frame which is close to 7 m/s and thus well within the requirements.

### Car does not have collisions.

For every point added to the trajectory the velocity is reduced when it gets closer than 30 meters to a vehicle. 30 meters is enough to get break from 50 mph to 0 with a negative acceleration of 10 m/s^2.

### The car stays in its lane, except for the time between changing lanes.

A cubic spline is used to interpolate the highway in-between the waypoints. Using frenet coordinates this approach allows for smooth trajectories that do not leave the lane. 

### The car is able to change lanes

Whenever a lane is detected with a lower cost, the car will try to reach it. The variable ref_lane will be set to the new lane. In the following iterations the spline curve interpolation will implicitly generate waypoints that smoothly switch to the new lane. The cost function for safety makes sure that the car
1. will not jump over a lane
2. will only change lane when there is free space available at the projected position of the trajectory

## Conclusion
The path planner is ale to drive quite efficiently and without collisions across the map. There are still some ways to improve the behavior. The car keeps breaking and accelerating when following a vehicle waiting for another lane to open up some free space for a lane switch. to solve this the existing trajectory would have to be recalculated in such that the vehicle reaches the blocking vehicle in a safe distance having the same velocity. Another approach would be to look further ahead than the horizon for generating the trajectory. Another potential issue which can happen rather rarely is a crash during lane changes. The projection of the free space in time could be improved.

