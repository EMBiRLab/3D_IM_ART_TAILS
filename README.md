# 3D_IM_ART_TAILS
 
## Description
This repository implements a generalized, optimization-based approach to quantitatively assess the performance of different inertial appendages in biomechanical and/or robotic systems. The approach leverages physics-based simulation and trajectory optimization techniques. Using this approach, one can study the impact of various inertial appendages on inertial maneuvering, gain functional insights into different inertial appendages, and inspire robotic appendage design. For more information, please refer to our accompanying paper:

## Dependencies
This repository has been verified on MATLAB R2022.

The repository depends on [Pinocchio](https://github.com/stack-of-tasks/pinocchio). To ensure the code functions as intended, please install and configure Pinocchio on your system.


## Repository Structure
The repository includes the following directories:

- **`DRD`**: Contains self-written MATLAB code for calculating rigid body dynamics based on Roy Featherstone's algorithm.
- **`Selfcollision_con_gen`**: Contains MATLAB code for generating constraints to prevent self-collisions between different parts of a rigid-body robotic model. 
- **`Target_trajectories`**: Contains target trajectories used for trajectory optimization.
- **`Trajectory_optimization`**: Contains MATLAB for performing trajecotry optimization.
- **`URDF`**: Contains URDF files for robot models with different inertial appendages. 
- **`Varied_length_derivative`**: Contains MATLAB code for calculating dynamics derivatives when the lengths of individual parts of the inertial appendages vary.

## Support and Contact
For support, questions, or feedback, please contact the authors directly.

## License
Licensed under the Creative Common Attribution 4.0 (CC-BY-4.0) License. For more information, see the [LICENSE](./LICENSE) file.

## Authors
Xun Fu (xunfu@umich.edu)
Bohao Zhang (jimzhang@umich.edu)
Ceri J. Weber (cweber@ucsd.edu)
Kimberly L. Cooper (kcooper@ucsd.edu)
Ram Vasudevan (ramv@umich.edu)
Talia Y. Moore (taliaym@umich.edu)

