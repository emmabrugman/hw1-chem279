# HW1
For modularity and reusability, I separated my code into several different files for each task. 

## 1.1 files
- lj_energy.cpp
- lj_energy.h
- main.cpp

## 1.2 files
- lj_force.cpp
- lj_force.h
- lj_force_plot.py
- main.cpp

### Explanation of the plot:
This plot shows the log of the error versus the log of the step size for both forward difference and central difference approximations of the force. Ideally the plot should display a slope for forward and central differences, however, my plot displays a slope of 0 which suggests that there is something wrong with my calculations or that the error is not decreasing with smaller step sizes. 

## 1.3 file
- steepest_descent.cpp
- steepest_descent.h
- main.cpp

## main.cpp
The main file includes calculations for 1.1, 1.2, 1.3. To view the computations, you will need to run this file.