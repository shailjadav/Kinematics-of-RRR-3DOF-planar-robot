% ME 639: Introduction to robotics
% Assignment 3: Problem 1
%               28 August 2018
% a. Formulate the forward kinematics transformation to calculate the end-effector position
%    and orientation.
% b. Program to solve for any two cases, assume L 1 = 4m, L 2 = 3m and L 3 = 2m.
%
%
% Author: Shail Jadav 18310039
%% Initialization
clear 
close all
clc
%%
%syms theta1 theta2 theta3 link1 link2 link3   

theta1=120; %input theta in degrees
theta2=70; %input theta in degrees
theta3=95; %input theta in degrees

link1=4; %Input the link length
link2=3; %Input the link length
link3=2; %Input the link length

theta1=theta1*(pi/180); %degree to radian conversion
theta2=theta2*(pi/180); %degree to radian conversion
theta3=theta3*(pi/180); %degree to radian conversion


% Homogeneus transformation matrix
H01 = [cos(theta1) -sin(theta1) 0 link1*cos(theta1);sin(theta1) cos(theta1) 0 link1*sin(theta1);0 0 1 0;0 0 0 1]; %Frame 0 to 1 tranformation
H12 = [cos(theta2) -sin(theta2) 0 link2*cos(theta2);sin(theta2) cos(theta2) 0 link2*sin(theta2);0 0 1 0;0 0 0 1]; %Frame 1 to 2 tranformation
H23 = [cos(theta3) -sin(theta3) 0 link3*cos(theta3);sin(theta3) cos(theta3) 0 link3*sin(theta3);0 0 1 0;0 0 0 1]; %Frame 2 to 3 tranformation

H03=H01*H12*H23;  %Frame 0 to 3 tranformation
H02=H01*H12;      %Frame 0 to 2 tranformation     
H01=H01;          %Frame 0 to 1 tranformation



%% Display The results
clc

O=[0,0];                   %Joint 1 position
P1=[H01(1,4) H01(2,4)];    %Joint 2 position
P2=[H02(1,4) H02(2,4)];    %Joint 3 position
P3=[H03(1,4) H03(2,4)];    %End effector position

Orn= atan2(H03(2,1),H03(1,1));  %Orientation of end effector
Orn=(Orn)*(180/pi);

fprintf('The end effector position x = %f y = %f' ,H03(1,4),H03(2,4))
fprintf('  and orientation= %f',Orn)

%% Plot results
figure
plot(P1(1),P1(2),'ok','LineWidth',5)
hold on
plot(P2(1),P2(2),'ok','LineWidth',5)
plot(P3(1),P3(2),'mX','LineWidth',10)
plot(0,0,'ok','LineWidth',10)
xlim([-10 10])
ylim([-10 10])
grid minor
plot([0 P1(1)], [0 P1(2)],'r','LineWidth',5)
plot([P1(1) P2(1)], [P1(2) P2(2)],'b','LineWidth',5)
plot([P2(1) P3(1)], [P2(2) P3(2)],'g','LineWidth',5)
title('3 DOF planar serial manipulator')
xlabel('X axis (m)')
ylabel('Y axis (m)')
set(gca,'FontSize',18)

