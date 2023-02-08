clear,clc
eta1 = 0;
eta2 = 0.4;
tspan = [0,5];
D1 = 0.01;
D2 = 0.01;
out1=-0.5:0.01:0.5;
out2=-0.5:0.01:0.5;
n1=length(out1);
n2=length(out2);
X1=zeros(4,n1*n2);
for i=1:n1
    out1(i)=-0.5+0.01*(i-1);
     for j=1:n2
         out2(j)=-0.5+(j-1)*0.01;
y0 = [0 out1(i) 0 out2(j)];
% out1 = -6.104436; out2 = -1.0454503;
% y0 = [-1 out1 0 out2];

%dy1 = @(t,y)[y(2);(y(1)-y(1)^3-y(1)*y(3)^2-eta1)*(1-3*y(1)^2-y(3)^2)+2*y(1)*y(3)*((1+y(1)^2)*y(3)+eta2)-4*y(1);y(4);-(y(1)-y(1)^3-y(1)*y(3)^2-eta1)*y(1)*y(3)+((1+y(1)^2)*y(3)+eta2)*(1+y(1)^2)-2*y(3)];
dy = @(t, y)[y(2); y(1)+y(1)*y(3)^2*(sin(y(3))^2/D1+cos(y(3))^2/D2)/(cos(y(3))^2/D1+sin(y(3))^2/D2);y(4);y(3)+sin(y(3))*cos(y(3))*(1/D2-1/D1)*(1-y(3)^2)/(sin(y(3))^2/D1+cos(y(3))^2/D2)];
%dy = @(t,y)[y(2);(y(1)-y(1)^3-y(1)*y(3)^2-eta1)*(1-3*y(1)^2-y(3)^2)+2*y(1)*y(3)*((1+y(1)^2)*y(3)+eta2)+3-y(1);y(4);-(y(1)-y(1)^3-y(1)*y(3)^2-eta1)*y(1)*y(3)+((1+y(1)^2)*y(3)+eta2)*(1+y(1)^2)-y(3)];

[x,u2] = ode45(dy,tspan,y0);
%plot(x,u2(:,1),'r');
%u2(:,1);   solution vector r 
%u2(:,3);   solution vector theta
X1(1,i*j)=out1(i); % initial velocity dot(r)
X1(2,i*j)=out2(j);  % initial velocity dot(theta)
%X1(3,i*j)=u2(end,1)*cos(u2(end,3)); %transform to E-coordinate-x
%X1(4,i*j)=u2(end,1)*sin(u2(end,3)); %transform to E-coordinate-y
X1(3,i*j)=u2(end,1);% the final point r at time T
X1(4,i*j)=u2(end,3);% the final point theta at time T
%plot(u2(:,1),u2(:,3),'g')
%plot(u2(:,1).*cos(u2(:,3)),u2(:,1).*sin(u2(:,3)))
     end 
end