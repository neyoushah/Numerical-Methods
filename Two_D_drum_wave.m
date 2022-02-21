function Two_D_drum_wave

clc; clear all; close all;
% This function solves the 2-D WAVE EQUATION U_tt = c^2(U_rr+(1/r) U_r+ (1/r^2) U_theta theta ) 
% Using second-order finite difference approximation for the partial derivatives
% The method is explicit---------------------------------------------------
ni=10;%%%%% If we increase ni and nj you should change dt in order to comply with Courant-Friedrich Stability Condition
nj=10;%--------------------------------------------------------------------
dx=1/ni;%------------------------------------------------------------------
dy=1/nj;%------------------------------------------------------------------
x = 0.1:dx:1;%%%%% ----x-------- will be  the radial coordinate------------
y = 0.1:dy:1;%%%%% ---y-----will be transformed to theta coordinate  ------
c =1;  % Wave velocity-----------------------------------------------------
sigma = 1/sqrt(2); gamma = 1/sqrt(2); %Courant-Friedrich Stability Condition
dt = sigma*(dx/c); %gamma=dt/(2*0.01*dx);
dt=0.001;%-----------------------------------------------------------------
t = 0:dt:1;  
u = zeros(length(x),length(y),length(t));
%--------------------------------------------------------------------------
for j=1:length(y)-1%-------------------------------------Initial conditions
%u(1,j,1)=-0.02;
u(1,j,1)=-0.5;
end
%__________________________________________________________________________
%p = 2; q = 1;--------------Another initial condition
%u(:,:,1) = transpose(0.02*sin(p.*pi.*x))*sin(q.*pi.*y); %u(x,y,0) = sin(p*pi*x)*sin(q*pi*y)
%__________________________________________________________________________
%__________________________________________________________________________
for i=2:length(x)-2 
    for j=2:length(y)-2
        u(i,j,2)=u(i,j,1)+(1/2)*((c^2*dt^2/dx^2)*(u(i+1,j,1)-2*u(i,j,1)+u(i-1,j,1))...
                +((c^2*dt^2)/((x(i))^2*dy^2))*(u(i,j+1,1)-2*u(i,j,1)+u(i,j-1,1))...
            +(c^2*dt^2/(2*x(i)*dx))*(u(i+1,j,1)-u(i-1,j,1))); 
    end
end
%__________________________________________________________________________
if i==length(x)-1 
    if j==length(y)-1
        u(i,j,2)=u(i,j,1)+(1/2)*((c^2*dt^2/dx^2)*(0-2*u(i,j,1)+u(i-1,j,1))...
                +((c^2*dt^2)/((x(i))^2*dy^2))*(0-2*u(i,j,1)+u(i,j-1,1))...
            +(c^2*dt^2/(2*x(i)*dx))*(0-u(i-1,j,1))); 
    end
end
%__________________________________________________________________________
for n=2:length(t)-1
    for i=2:length(x)-2
        for j=2:length(y)-2
            u(i,j,n+1)=2*u(i,j,n)-u(i,j,n-1)+(c^2*dt^2/dx^2)*(u(i+1,j,n)-2*u(i,j,n)+u(i-1,j,n))...
                +((c^2*dt^2)/((x(i))^2*dy^2))*(u(i,j+1,n)-2*u(i,j,n)+u(i,j-1,n))...
            +(c^2*dt^2/(2*x(i)*dx))*(u(i+1,j,n)-u(i-1,j,n)); 
        end
    end
%__________________________________________________________________________    
    if i==length(x)-1
        if j==length(y)-1
            u(i,j,n+1)=2*u(i,j,n)-u(i,j,n-1)+(c^2*dt^2/dx^2)*(0-2*u(i,j,n)+u(i-1,j,n))...
                +((c^2*dt^2)/((x(i))^2*dy^2))*(0-2*u(i,j,n)+u(i,j-1,n))...
            +(c^2*dt^2/(2*x(i)*dx))*(0-u(i-1,j,n)); 
        end
    end
%__________________________________________________________________________
end
%__________________________________________________________________________
% 
[T,R] = meshgrid(linspace(0.1,2*pi,10),linspace(0.1,1,10));
X = R.*cos(T);
Y = R.*sin(T);
figure
for j=1:length(t)
       colormap(cool);
       p1 = surf(X,Y,u(:,:,j));  view(15, 35)    
       title(sprintf('Two-D wave equation polar coordinates  at t = %1.2f ',t(j)),'Fontsize',11, 'interpreter', 'latex');
       xlabel('r','Fontsize',11); ylabel('\theta','Fontsize',11);
       zlabel(sprintf('u(r,\theta,t = %1.2f)',t(j)),'Fontsize',11);
       axis ([-1 1 0 2*pi -0.0025 0.0025]);
       pause(0.001);
       hold on; 
       if(j~=length(t))
       delete(p1);
       end
end
