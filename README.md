# Ababu
This is implementation of Robot arm Dynamics/Kinematics.
compute serial link robot joint torque in real-time.
Trajectory planing.
Jacobian of n-DOF robot arm.
It is only starter code for those struggling to figure out dyanmics of robot arm modeling as well as FKI
I have implemented impedance control  as a demo useage case. 
But becarfule the physical-parameters of your robot are diffrent from the robot i have used. change DH-parameters and Inertia property to yours.
This class can be used in real time but i belive there is still much sliming work to avoid watchdog for RTX and other RTOS technologies 
