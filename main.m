close all
clear all
clc
addpath('./casadi-windows-matlabR2016a-v3.5.5')
import casadi.*

%% system setting
A=[-1 1;
  	0 1];
B=[1,3]';
% define the state and input variable
x=SX.sym('x',2);
u=SX.sym('u',1);
% define the dynamics
f=A*x+B*u;
dyn = Function('dyn', {x, u}, {f}, {'X','U'}, {'f'});

%% setup the cost function and the constraint function
features=[x(1)^2, x(2)^2, u^2]';
weights=[0.1,0.3,0.6]';
% cost function
cost= Function('cost',{x,u},{weights'*features}, {'X','U'}, {'c'});
% constraint function
Q1 = eye(2);
R1 = eye(1);
constraint_exp = x'*Q1*x + u'*R1*u;
constraint = Function('constraint', {x,u}, {constraint_exp}, {'X','U'}, {'cons'});
d = 70;

%% use the oc solver the solve the optimal control probme
x0=[0.01,-0.01]';
T=50;
sol=OCsolver_IntegralConstraint(x0,T,dyn,cost,constraint,d);
traj_x=sol.x;
traj_x(:,end)=[];
traj_u=sol.u;
traj_t=0:T;

%% do the plot 
figure(1)
subplot(2,1,1)
plot(traj_t,traj_x(1,:),'LineWidth',3)
xlim([0,51])
hold on
plot(traj_t,traj_x(2,:),'LineWidth',3)
grid on
ylabel('$x$','interpreter','latex')
legend('$x_1$', '$x_2$','interpreter','latex')
subplot(2,1,2)
plot(traj_t,traj_u,'LineWidth',3)
xlim([0,51])
grid on
ylabel('$u$','interpreter','latex')
xlabel('time')
legend('$u$','interpreter','latex')
saveas(gcf,'trajectory.png');
clc;