% c ompute pdf of fpe corr to Levy noise with absorbing BC 
% in a bounded, symmetric domain (-r,r);
% central differencing for 2nd-order deriv 
% and using one-sided differencing near the boundaries;
% 3rd-order WENO for 1st-order deriv
% 3rd-order TVD RK in time


tic;

% function [U]=tgao7_Dwell_true(Alp,T)
clear ;
clc;





alp0=0.5;alp1 =0.5; dalp = 1.98/248; PP = (alp1 - alp0)/dalp; 
alpha = alp0 : dalp : alp1;

for m = 1:PP+1
T =5;
eps =0.25;
 a1 =0*pi;  b1=pi;
% %hh=pi/8; DD=(b-a)/hh;
%  %p2=zeros(16,1);
%  for i=1:16
%   a1(i)=pi*(i-1)/8;
%   b1(i)=pi*i/8;

eps = eps^(alpha(m))*(2/(b1-a1))^alpha(m);
d =0.25;
d = d*4/(b1-a1)^2;
h=1/5;
delta=(b1-a1)*h;

%k_f = 6;   K_d = 10;   k_d = 1;  R_b = 0.4;
%f = @(x)       (k_f*x.^2./(x.^2+K_d) - k_d*x + R_b)*2/(b1-a1);
f = @(x) 0;
rr=1;             %right
rl=1;             %left
fpmax=abs(f(b1));

h=1/5;
J=rr/h;             %right
L=rl/h;             %left
dt = 0.5*h^2;
% if alpha(m) <=1
% dt = 0.5*h;
% else 
% dt = 0.5*h^alpha(m); 
% end
dtdx = dt/h;
%T=2;

Jt=L+J-1;     %Jt=total number of unknowns  
x=-(rl+rr):h:(rl+rr);
C=-zeta(alpha(m)-1)*h^(2-alpha(m));          % correction term u''(x)
cons=alpha(m)*gamma((1+alpha(m))/2)/(2^(1-alpha(m))*sqrt(pi)*gamma(1-alpha(m)/2));
% coeff for the diffusion and the correction term
a=zeros(Jt,1);
c=zeros(Jt,1); 

Chh=( d/2 + cons*C*eps )/h^2;
c1=eps*cons/alpha(m);
c2=eps*cons*h;

b=zeros(Jt,1);
%nonintegral part

% coefficient of U_j
b(2:Jt-1) = (-2*Chh - c1*(1./(x(J+3:2*J+L-1)+rl).^alpha(m)))';
% one-sided diff near boundaries
b(1)  = 2*Chh - c1*(1/(rl+x(J+2))^alpha(m)); 
b(Jt) = 2*Chh - c1*(1/(rl+x(2*J+L))^alpha(m)); % one-sided diff
a= Chh*ones(Jt,1);  % coefficient of U_(j-1)
c= Chh*ones(Jt,1);  % coefficient of U_(j+1) 
c(1)  = -5*Chh; % one-sided diff
a(Jt) = -5*Chh; % one-sided diff
vp2 = zeros(Jt,1); vp2(3) = 4*Chh;  % one-sided diff
vp3 = zeros(Jt,1); vp3(4) =  -Chh; % one-sided diff 
vm2 = zeros(Jt,1); vm2(Jt-2) = 4*Chh;  % one-sided diff
vm3 = zeros(Jt,1); vm3(Jt-3) =  -Chh; % one-sided diff 

%only brownian motion
% b(1)  = -2*Chh - c1*(1/(rl+x(J+2))^alpha(m)+1/(rr-x(J+2))^alpha(m)); 
% b(Jt) = -2*Chh - c1*(1/(rl+x(2*J+L))^alpha(m)+1/(rr-x(2*J+L))^alpha(m)); % one-sided diff
% a = Chh*ones(Jt,1);  % coefficient of U_(j-1)
% c = Chh*ones(Jt,1);  % coefficient of U_(j+1) 
% vp2 = zeros(Jt,1); % one-sided diff
% vp3 = zeros(Jt,1); % one-sided diff 
% vm2 = zeros(Jt,1); % one-sided diff
% vm3 = zeros(Jt,1); % one-sided diff 

% integral part
for j=-L+1:J-1
   b(j+L)= b(j+L) - c2*( sum(1./abs(x(J+2-j:L+J)).^(1+alpha(m))) ...
                       + sum(1./abs(x(L+J+2:2*J+L-j)).^(1+alpha(m))) ...
         + .5/abs(x(J+1-j))^(1+alpha(m)) + .5/abs(x(2*J+L+1-j))^(1+alpha(m)) );  
end 
A=spdiags([vm3 vm2 [a(2:end); 0] b ...
           [0; c(1:end-1)] vp2 vp3],-3:3,Jt,Jt);

% coefficient of u_(j+k) 
B=zeros(size(A));
for j=-L+1:J-1
  B(L+j,:)=[1./abs(x(J+2-j:L+J)).^(1+alpha(m))  0  1./abs(x(L+J+2:2*J+L-j)).^(1+alpha(m))];
end


%A = diag(ones(Jt,1))+ dt*(A+c2*B); % the iterative matrix for each time-step
rx =2;
X=-rl:h:rr;
XX = a1 : (b1-a1)/2*h :b1;
%UU=sqrt(pi)*(r^2-X.^2).^(alpha/2)/(2^alpha*gamma(1+alpha/2)*gamma(1/2+alpha/2)); 
UU=sqrt(40/pi)*exp(-(XX - rx).^2 *40); %gaussian
% UU = 0.01./(pi*(0.01^2+X.^2));
%UU = zeros(2*J+1,1)';
%lo = (x0+1)/h+1;
%UU(lo) = 1/(2*h);    UU(lo+1) = 1/(2*h);  UU(lo-1) = 1/(2*h);
%UU=ones(size(X))./(rr+rl);
U=UU(2:end-1)';
Un=U; U1=U; U2=U;
%U=UU(2:end-1);

nft=round(T/dt);
%figure
%hold on
nu = length(U);
data = zeros(nu+4,1);


for nt=1:nft-1

    U1 = U + dt*(A+c2*B)*U;
    % global Lax-Friedrichs(LF) flux splitting
    data(3:nu+2) = (f(XX(2:end-1)) + fpmax)'.*U/2;
    fx1 = derWENOr2_minus(data,h);
    data(3:nu+2) = (f(XX(2:end-1)) - fpmax)'.*U/2;
    fx2 = derWENOr2_plus(data,h);
    U1 = U1 - dtdx*(fx1+fx2);
    
    U2 = 0.75*U + U1/4 + (dt/4)*(A+c2*B)*U1;
    data(3:nu+2) = (f(XX(2:end-1)) + fpmax)'.*U1/2;
    fx1 = derWENOr2_minus(data,h);
    data(3:nu+2) = (f(XX(2:end-1)) - fpmax)'.*U1/2;
    fx2 = derWENOr2_plus(data,h);
    U2 = U2 - (dtdx/4)*(fx1+fx2);
    
    Un = U/3 + 2*U2/3 + (2*dt/3)*(A+c2*B)*U2;
    data(3:nu+2) = (f(XX(2:end-1)) + fpmax)'.*U2/2;
    fx1 = derWENOr2_minus(data,h);
    data(3:nu+2) = (f(XX(2:end-1)) - fpmax)'.*U2/2;
    fx2 = derWENOr2_plus(data,h);
    Un = Un - (2*dtdx/3)*(fx1+fx2);
    
%     if mod(nt,100)==0
%        plot(X,[0; U; 0])    
%     end
         %pause
            
     U=Un;
     %t=nt*dt;
     
     
%      [max1, loc1]= max(U);
%      mp = loc1;
%      MP(nt) = XX(mp);


     % F(nt) = getframe;
      XU(:,nt) = U;
         
end
%save('plotdiffT_uniform','P1','P2','P3','P4')
%title(['final time t= ', num2str(t)])
% Ux = U; 
t = dt : dt : (nft-1)*dt;
p2=sum(XU(2:end,:)*delta);




  ppp=sum(XU(end,:)*delta);

% ppp

%  figure

 %% plot(XX(2:end-1),pp2,'b-')

%最大可能相图
% hold on;
% plot(t, MP,'k');
% xlabel('t');  ylabel('MP');
% [y1,nn]=min(MP),t(nn)
% % [y2,mm]=max(MP),t(mm)
% max(MP);
%figure  
%%%mesh(t,XX(2:end-1),XU);

%entropy
% ab(m) = max(MP);
% plot(rx,ab)

% [y1,nn] = max(aa);
% SS(m) = alpha(nn);
% SA(m) = alpha(y1);
 end
%plot(alpha, p2,'m-.')
 plot(t,p2,'k-')
 hold on
% figure
% aa = -1/h*sum(XU.*log(XU));
% plot(alpha, SS,'b.');

% save fpkalp01  alpha  SS   h  SA


%  save  fpkg05  rx  ab  XX  t MP  U
% m=find(U==a)
% figure
% plot(XX(2:end-1)',U,'g--')
% [y,l]=max(U),XX(l)

 toc;
%  MM=sum(U)*h
%  v1=abs(U);
%  save tgao7_Dwell_true.dat v1 -ascii

%  T = 30; h=1/600; alpha = 0.5; dt = 0.5*h^alpha;  nft=round(T/dt); t = dt : dt : (nft-1)*dt;
% t1 = 20; n1 = round(t1/dt); t2 = 25 ; n2 = round(t2/dt); plot(XX(2:end-1)',XU(:,n1),'r-'); 
% hold on; plot(XX(2:end-1)',XU(:,n2),'g-'); hold on;  plot(XX(2:end-1)',XU(:,end),'b-');
