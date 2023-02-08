
eps1 = 0.01;
eps2=0.01;

%% Lagrangian
f = @(x,y) x;
fx = @(x,y) 1;
g = @(x,y) y;
s = @(x,y) sqrt(2*eps1)*cos(y)+sqrt(2*eps2)*sin(y);%Brownian motion of \dot x
gy = @(x,y) 1;
ga = @(x,y) -sqrt(2*eps1)*sin(y)/x+sqrt(2*eps2)*cos(y)/x;%Brownian motion of \dot y

L = @(x,y,u,v) (u-f(x,y))^2/(s(x,y))^2 + fx(x,y) + (v-g(x,y))^2/(ga(x,y))^2 + gy(x,y);

%% Compute the action functional
T = 3;
% uu = load('subdatau_opt.mat');
% vv = load('subdatav_opt.mat');
% u = uu.x1n;   
% v = vv.x3n;
% Act=ActionValue(u2(:,1),u2(:,3),L,T);
% [m n] = size(u);
 for i = 1 :l
     Act(i) = ActionValue(P1{i},P2{i},L,T);
 end
 ind = find(Act==min(min(Act)));
 
 plot(P1{ind}.*cos(P2{ind}),P1{ind}.*sin(P2{ind}),'b','LineWidth',1.3)
hold on;
w=0:0.01:2*pi;
plot(cos(w),sin(w),'g','LineWidth',1.3)