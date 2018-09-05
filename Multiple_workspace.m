% ME 639: Introduction to robotics
% Assignment 3: Problem 3
%               28 August 2018
%
% Author: Shail Jadav 18310039
%
% Workspace is defined as the set of points reachable by the end-effector of a
% manipulator.
% a. Plot the workspace of the 3R manipulator when ? = 0 0 , L 1 = 4m, L 2 = 3m and L 3 = 2m.
% b. Plot the workspace of the 3R manipulator when ? = 45 0 , L 1 = 4m, L 2 = 3m and L 3 = 2m.
% c. Is it possible to calculate the workspace of the 3R manipulator when there is no
% constraint on ?? If so then plot it. Discuss your results in comparison with part (a and b).

%% Initialization
clear
close all
clc
%% Workspace calculation using Phi=45
M=1;
for i=-10:0.2:10
    for j=-10:0.2:10
        X=i; % X cordinate value
        Y=j; % Y cordinate value
        phi = 45; % Phi value
        
        link1 =4; %link dimentions
        link2 =3; %link dimentions
        link3 =2; %link dimentions
        
        phi = deg2rad(phi);
        
        nx = X - link3*cos(phi); %Calculations for nx and ny
        ny = Y - link3*sin(phi);
        
        delta = nx^2 + ny^2;
        c2 = ( delta -link1^2 -link2^2)/(2*link1*link2);
        s2 = sqrt(1-c2^2);
        
        if(imag(s2)==0)
            theta_2 = atan2(s2, c2); %calci=ultions for valid theta 2
            
            s1 = ((link1+link2*c2)*ny - link2*s2*nx)/delta;
            c1 = ((link1+link2*c2)*nx + link2*s2*ny)/delta;
            theta_1 = atan2(s1,c1);
            theta_3 = phi-theta_1-theta_2;
            
            theta1=theta_1*(180/pi);
            theta2=theta_2*(180/pi);
            theta3=theta_3*(180/pi);
            WS45(M,:)=num2cell([i j theta1 theta2 theta3]);
            M=M+1;
        end
    end
end
M=1;
for i=-10:0.2:10
    for j=-10:0.2:10
        X=i; % X cordinate value
        Y=j; % Y cordinate value
        phi = 60; % Phi value
        
        link1 =4; %link dimentions
        link2 =3; %link dimentions
        link3 =2; %link dimentions
        
        phi = deg2rad(phi);
        
        nx = X - link3*cos(phi); %Calculations for nx and ny
        ny = Y - link3*sin(phi);
        
        delta = nx^2 + ny^2;
        c2 = ( delta -link1^2 -link2^2)/(2*link1*link2);
        s2 = sqrt(1-c2^2);
        
        if(imag(s2)==0)
            theta_2 = atan2(s2, c2); %calci=ultions for valid theta 2
            
            s1 = ((link1+link2*c2)*ny - link2*s2*nx)/delta;
            c1 = ((link1+link2*c2)*nx + link2*s2*ny)/delta;
            theta_1 = atan2(s1,c1);
            theta_3 = phi-theta_1-theta_2;
            
            theta1=theta_1*(180/pi);
            theta2=theta_2*(180/pi);
            theta3=theta_3*(180/pi);
            WS60(M,:)=num2cell([i j theta1 theta2 theta3]);
            M=M+1;
        end
    end
end
M=1;
for i=-10:0.2:10
    for j=-10:0.2:10
        X=i; % X cordinate value
        Y=j; % Y cordinate value
        phi = 90; % Phi value
        
        link1 =4; %link dimentions
        link2 =3; %link dimentions
        link3 =2; %link dimentions
        
        phi = deg2rad(phi);
        
        nx = X - link3*cos(phi); %Calculations for nx and ny
        ny = Y - link3*sin(phi);
        
        delta = nx^2 + ny^2;
        c2 = ( delta -link1^2 -link2^2)/(2*link1*link2);
        s2 = sqrt(1-c2^2);
        
        if(imag(s2)==0)
            theta_2 = atan2(s2, c2); %calci=ultions for valid theta 2
            
            s1 = ((link1+link2*c2)*ny - link2*s2*nx)/delta;
            c1 = ((link1+link2*c2)*nx + link2*s2*ny)/delta;
            theta_1 = atan2(s1,c1);
            theta_3 = phi-theta_1-theta_2;
            
            theta1=theta_1*(180/pi);
            theta2=theta_2*(180/pi);
            theta3=theta_3*(180/pi);
            WS90(M,:)=num2cell([i j theta1 theta2 theta3]);
            M=M+1;
        end
    end
end
M=1;
for i=-10:0.2:10
    for j=-10:0.2:10
        X=i; % X cordinate value
        Y=j; % Y cordinate value
        phi = 180; % Phi value
        
        link1 =4; %link dimentions
        link2 =3; %link dimentions
        link3 =2; %link dimentions
        
        phi = deg2rad(phi);
        
        nx = X - link3*cos(phi); %Calculations for nx and ny
        ny = Y - link3*sin(phi);
        
        delta = nx^2 + ny^2;
        c2 = ( delta -link1^2 -link2^2)/(2*link1*link2);
        s2 = sqrt(1-c2^2);
        
        if(imag(s2)==0)
            theta_2 = atan2(s2, c2); %calci=ultions for valid theta 2
            
            s1 = ((link1+link2*c2)*ny - link2*s2*nx)/delta;
            c1 = ((link1+link2*c2)*nx + link2*s2*ny)/delta;
            theta_1 = atan2(s1,c1);
            theta_3 = phi-theta_1-theta_2;
            
            theta1=theta_1*(180/pi);
            theta2=theta_2*(180/pi);
            theta3=theta_3*(180/pi);
            WS180(M,:)=num2cell([i j theta1 theta2 theta3]);
            M=M+1;
        end
    end
end
M=1;
for i=-10:0.2:10
    for j=-10:0.2:10
        X=i;
        Y=j;
        phi = 0;
        
        link1 =4;
        link2 =3;
        link3 =2;
        
        phi = deg2rad(phi);
        nx = X - link3*cos(phi);
        ny = Y - link3*sin(phi);
        
        delta = nx^2 + ny^2;
        c2 = ( delta -link1^2 -link2^2)/(2*link1*link2);
        s2 = sqrt(1-c2^2);
        if(imag(s2)==0)
            theta_2 = atan2(s2, c2);
            
            s1 = ((link1+link2*c2)*ny - link2*s2*nx)/delta;
            c1 = ((link1+link2*c2)*nx + link2*s2*ny)/delta;
            theta_1 = atan2(s1,c1);
            theta_3 = phi-theta_1-theta_2;
            
            theta1=theta_1*(180/pi);
            theta2=theta_2*(180/pi);
            theta3=theta_3*(180/pi);
            WS0(M,:)=num2cell([i j theta1 theta2 theta3]);
            M=M+1;
        end
    end
end

Xdata45=cell2mat(WS45(:,1));
Ydata45=cell2mat(WS45(:,2));
figure
for i=1:1:length(Xdata45)
    plot(Xdata45(i,1),Ydata45(i,1),'*r')
    hold on
end
set(gca);
title('Workspace with Phi=45')
xlabel('X axis (m)')
ylabel('Y axis (m)')
grid minor
set(gca,'FontSize',18);
xlim([-10 10])
ylim([-10 10])

Xdata60=cell2mat(WS60(:,1));
Ydata60=cell2mat(WS60(:,2));
figure
for i=1:1:length(Xdata60)
    plot(Xdata60(i,1),Ydata60(i,1),'*r')
    hold on
end
set(gca);
title('Workspace with Phi=60')
xlabel('X axis (m)')
ylabel('Y axis (m)')
grid minor
set(gca,'FontSize',18);
xlim([-10 10])
ylim([-10 10])

Xdata90=cell2mat(WS90(:,1));
Ydata90=cell2mat(WS90(:,2));
figure
for i=1:1:length(Xdata90)
    plot(Xdata90(i,1),Ydata90(i,1),'*r')
    hold on
end
set(gca);
title('Workspace with Phi=90')
xlabel('X axis (m)')
ylabel('Y axis (m)')
grid minor
set(gca,'FontSize',18);
xlim([-10 10])
ylim([-10 10])

Xdata0=cell2mat(WS0(:,1));
Ydata0=cell2mat(WS0(:,2));
figure
for i=1:1:length(Xdata0)
    plot(Xdata0(i,1),Ydata0(i,1),'*g')
    hold on
end
set(gca);
title('Workspace  with Phi=0')
xlabel('X axis (m)')
ylabel('Y axis (m)')
grid minor
set(gca,'FontSize',18);
xlim([-10 10])
ylim([-10 10])

Xdata180=cell2mat(WS180(:,1));
Ydata180=cell2mat(WS180(:,2));
figure
for i=1:1:length(Xdata180)
    plot(Xdata180(i,1),Ydata180(i,1),'*r')
    hold on
end
set(gca);
title('Workspace with Phi=180')
xlabel('X axis (m)')
ylabel('Y axis (m)')
grid minor
set(gca,'FontSize',18);
xlim([-10 10])
ylim([-10 10])

figure
for i=1:1:length(Xdata0)
    plot3(Xdata0(i,1),Ydata0(i,1),0,'*r')
    hold on
end
for i=1:1:length(Xdata45)
    plot3(Xdata45(i,1),Ydata45(i,1),45,'*g')
    hold on
end
for i=1:1:length(Xdata60)
    plot3(Xdata60(i,1),Ydata60(i,1),60,'*b')
    hold on
end
for i=1:1:length(Xdata90)
    plot3(Xdata90(i,1),Ydata90(i,1),90,'*k')
    hold on
end
for i=1:1:length(Xdata180)
    plot3(Xdata180(i,1),Ydata180(i,1),180,'*m')
    hold on
end
set(gca);
title('Workspace')
xlabel('X axis (m)')
ylabel('Y axis (m)')
zlabel('Phi (degree)')
grid minor
set(gca,'FontSize',18);
xlim([-10 10])
ylim([-10 10])

