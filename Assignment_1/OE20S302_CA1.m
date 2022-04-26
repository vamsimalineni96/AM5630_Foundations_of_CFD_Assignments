% Plot temperature vs Length graph at : 
% t = 0 ; 0.5 ; 1 ; 2 ; 5 ; 10 s
% dt= 0.1 ; 0.01 ; 0.001  

clc;
clear;

l=1;                                                % Length of the wire
t=10;                                               % Final time
a=input('Enter the value of thermal diffusivity:'); % Thermal Diffusivity 
dt=input('Enter the value of time spacing(dt):');   % Time spacing
dx=input('Enter the value of grid spacing(dx):');   % Grid spacing


timesteps=t/dt;     % Number of time steps       
nx=(l/dx)+1;        % Number of space steps for ftcs scheme
g=a*dt/(dx*dx);     % Constant in FDM form

% Grid Generation and Initial Conditions :

for ss=1:nx
    x(ss)=(ss-1)*dx;
    T(ss,1)=0;
end

% Boundary Conditions :

for k=2:timesteps+1
    T(1,k)=1;
    T(nx,k)=0;
    time(k)=(k-1)*dt;
end

% Scheme Implementation :
for ts=1:timesteps
    for ss=2:nx-1
        T(ss,ts+1)=T(ss,ts)+g*(T(ss-1,ts)+T(ss+1,ts)-2*T(ss,ts));
    end
end

% Plotting the Graphs :
if(dt==0.1)
    figure(1)
    plot(x,T(:,1),'ro',x,T(:,5),x,T(:,10),x,T(:,20),x,T(:,50),x,T(:,100))
    legend("t=0s","t=0.5s","t=1s","t=2s","t=5s","t=10s")
    xlabel("Distance(m)")
    ylabel("Temperature (^0C)")
    title("Temperature vs Distance for dt=0.1")
end

if(dt==0.01)
    figure(2)
    plot(x,T(:,1),'ro',x,T(:,50),x,T(:,100),x,T(:,200),x,T(:,500),x,T(:,1000))
    legend("t=0s","t=0.5s","t=1s","t=2s","t=5s","t=10s")
    xlabel("Distance(m)")
    ylabel("Temperature (^0C)")
    title("Temperature vs Distance for dt=0.01")
end

if(dt==0.001)
    figure(3)
    plot(x,T(:,1),'ro',x,T(:,500),x,T(:,1000),x,T(:,2000),x,T(:,5000),x,T(:,10000))
    legend("t=0s","t=0.5s","t=1s","t=2s","t=5s","t=10s")
    xlabel("Distance(m)")
    ylabel("Temperature (^0C)")
    title("Temperature vs Distance for dt=0.001")
end

