clc
clear
close all

global u1 u2 u3

u1 = 0; u2 =0; u3 = 0;

x0 = [0; pi/2; 0; 0; 0; 0];

t_span = 0:0.01:1;
options = odeset('maxstep',0.001);
[t,x] = ode45(@nonlinear_dynamic,t_span,x0,options);

theta1 = x(:,1); theta2 = x(:,2); theta3 = x(:,3); thetad4 = x(:,4); thetad5 = x(:,5); thetad6 = x(:,6);

figure
hold on
plot(t,theta1,'r-','linewidth',2)
plot(t,theta2,'b-.','linewidth',2)
plot(t,theta3,'g--','linewidth',2)
legend('Theta1','Theta2','Theta3')
set(gca,'fontsize',12);
xlabel('Time (s)','fontsize',12);
ylabel('Joint Angles (rad)','fontsize',12);