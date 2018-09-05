% ME 639: Introduction to robotics
% Assignment 3: Problem 1
%               28 August 2018
%
% Author: Shail Jadav 18310039
% a. Formulate the forward kinematics transformation to calculate the end-effector position
%    and orientation.
% b. Program to solve for any two cases, assume L 1 = 4m, L 2 = 3m and L 3 = 2m.
%% Initialization
clear 
close all
clc
 
%% Input cordinates and Phi
X =-2;
Y = 3;
phi = 90;

link1 =4; %link1 length
link2 =3; %link2 length
link3 =2; %link3 length


%% Inverse kinematics calculations
phi = deg2rad(phi);

nx = X - link3*cos(phi); %Calculations for nx and ny
ny = Y - link3*sin(phi);

delta = nx^2 + ny^2;     %calculations for theta2
c2 = ( delta -link1^2 -link2^2)/(2*link1*link2);
s2 = sqrt(1-c2^2);
if(imag(s2) ~= 0)
    error('The given point cannot be reached please try diffrent value') 
end
theta_2 = atan2(s2, c2);

s1 = ((link1+link2*c2)*ny - link2*s2*nx)/delta; %calculation for theta1
c1 = ((link1+link2*c2)*nx + link2*s2*ny)/delta;
theta_1 = atan2(s1,c1);


theta_3 = phi-theta_1-theta_2;  %calculation for theta3

%% Dispaly the results
fprintf('The theta1 = %f  theta2 = %f  theta3 = %f ' ,theta_1*(180/pi),theta_2*(180/pi),theta_3*(180/pi))


%% Validation by forward kinematics
theta1=theta_1;
theta2=theta_2;
theta3=theta_3;

link1=4;
link2=3;
link3=2;

%% Homogeneus transformation matrix
H01 = [cos(theta1) -sin(theta1) 0 link1*cos(theta1);sin(theta1) cos(theta1) 0 link1*sin(theta1);0 0 1 0;0 0 0 1]; %Frame 0 to 1 tranformation
H12 = [cos(theta2) -sin(theta2) 0 link2*cos(theta2);sin(theta2) cos(theta2) 0 link2*sin(theta2);0 0 1 0;0 0 0 1]; %Frame 1 to 2 tranformation
H23 = [cos(theta3) -sin(theta3) 0 link3*cos(theta3);sin(theta3) cos(theta3) 0 link3*sin(theta3);0 0 1 0;0 0 0 1]; %Frame 2 to 3 tranformation

H03=H01*H12*H23;  %Frame 0 to 3 tranformation
H02=H01*H12;      %Frame 0 to 2 tranformation     
H01=H01;          %Frame 0 to 1 tranformation



%% Display The results


O=[0,0];                   %Joint 1 position
P1=[H01(1,4) H01(2,4)];    %Joint 2 position
P2=[H02(1,4) H02(2,4)];    %Joint 3 position
P3=[H03(1,4) H03(2,4)];    %End effector position

Orn= atan2(H03(2,1),H03(1,1));  %Orientation of end effector
Orn=(Orn)*(180/pi);
%% Plot results of Forward kinematics
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


