
% Computer Assignment 3 - AM5630
% Vamsi Sai Krishna Malineni - OE20S302
% Lid Driven Cavity Problem - Stream Function-Vorticity Approach
% Unsteady Navier Stokes Equation is solved for incompressible fluid with
% constant properties

clear all;
clc;
tic;
% Problem Definition
a = 1;                  % Side Length
u_lid = 1;             % Lid Velocity
rho = 1;                % Density
mu = 0.01;              % Dynamic Viscosity.
nu=mu/rho;               % Kinematic Visocity
dt = 0.001;             % Time step
max_iter = 100000;          % Maximum Number of Iterations
max_error = 1e-7;            % Error criteria for convergence

% Ghia .et al, Centerline velocity data
v_ghia = [0.00000 0.09233 0.10091 0.1089 0.12317 0.16077 0.17507 0.17527 0.05454 -0.24533 -0.22445 -0.16914 -0.10313 -0.08864 -0.07391 -0.05906 0.00000];
x_ghia = [1 9 10 11 13 21 30 31 65 104 111 117 122 123 124 125 129]/129; 

u_ghia = [0.00000 -0.03717 -0.004192 -0.04775 -0.06434 -0.10150 -0.15662 -0.21090 -0.20581 -0.13641 0.00332 0.23151 0.68717 0.737722 0.78871 0.84123 1.00000];
y_ghia = [1 8 9 10 14 23 37 59 65 80 95 110 123 124 125 126 129]/129;

% Defining mesh parameters
nx = 51;  
ny = nx;  

dx=(a/(nx-1));
dy=(a/(ny-1));

x = 0:dx:a; 
y = 0:dx:a;  

im = 1:nx-2; 
i = 2:nx-1; 
ip = 3:nx;  
jm = 1:ny-2; 
j = 2:ny-1; 
jp = 3:ny;

p = zeros(nx,ny);    % Pressure 
pn = zeros(nx,ny);

% Assesing the stability criteria :
st=(nu*dt/(dx^2))+(nu*dt/(dy^2));

if st <= 0.5
    disp("CFL Criteria is met: FTCS is stable ")    
    
    % Allocating matrix sizes for various parameters.
    psi = zeros(nx,ny);         % Stream Function
    
    omega_n1 = zeros(nx,ny);    % Vorticity at time step n+1
    omega_n = zeros(nx,ny);     % Vorticity at time step n
    
    u = zeros(nx,ny);         % X - direction velocity
    v = zeros(nx,ny);         % Y - direction velocity
    
    rhs = zeros(nx,ny);         % RHS of pressure poisson equation
    
    % Algorithm Starts

    for iter = 1:max_iter     
        % Boundary Conditions for Vorticity.
        omega_n1(1:nx,ny) = -2*psi(1:nx,ny-1)/(dx^2) - u_lid*2/dx;        % Top
        omega_n1(1:nx,1)  = -2*psi(1:nx,2)   /(dx^2);                     % Bottom
        omega_n1(1,1:ny)  = -2*psi(2,1:ny)   /(dx^2);                     % Left
        omega_n1(nx,1:ny) = -2*psi(nx-1,1:ny)/(dx^2);                     % Right
    
    % Partially solving the Vorticity Transport Equation
        omega_n = omega_n1;
        omega_n1(i,j) = omega_n(i,j) + (-1*(psi(i,jp)-psi(i,jm))/(2*dx) .* (omega_n(ip,j)-omega_n(im,j))/(2*dx)+...
                        (psi(ip,j)-psi(im,j))/(2*dx) .* (omega_n(i,jp)-omega_n(i,jm))/(2*dx)+...
                        (mu/rho)*(omega_n(ip,j)+omega_n(im,j)-4*omega_n(i,j)+omega_n(i,jp)+omega_n(i,jm))/(dx^2))*dt;
            
    % Partially solving the Vorticity Elliptic Equation to obtain Stream Function(phi)
        psi(i,j) = (omega_n1(i,j)*dx^2 + psi(ip,j) + psi(i,jp) + psi(i,jm) + psi(im,j))/4;
    
    % Convergence Criteria
        if iter > 10
            error = max(max(omega_n1 - omega_n));
            if error < max_error
                disp("Iterations to Steady State")
                disp(iter)
                break;
            end
        end
end
 
% Calculating velocity fields from Stream Function
    u(2:nx-1,ny) = u_lid; 
    u(i,j) = (psi(i,jp)-psi(i,jm))/(2*dx);  
    v(i,j) = (-psi(ip,j)+psi(im,j))/(2*dx);
    

    
    const=2*rho/(4*dx*dx);
% Calculating Pressure field throughout the cavity
  for k=1:max_iter
   
for i=2:(nx-1)
    for j=2:(ny-1)
     rhs(i,j)=const*(((u(i+1,j)-u(i-1,j))*(v(i,j+1)-v(i,j-1)))-((v(i+1,j)-v(i-1,j))*(u(i,j+1)-u(i,j-1))));
   end
end

    for i=2:nx-1
        for j=2:ny-1
            p(i,j)=(.25*(pn(i+1,j)+pn(i-1,j) + pn(i,j+1)+pn(i,j-1)))-(0.25*(rhs(i,j)*dx^2));
        end
    end

%    Boundary conditions
     for i=1:nx
    p(i,1)=p(i,2);
    p(i,ny)=p(i,ny-1);
    end
    for j=1:ny
    p(nx,j)=p(nx-1,j);
    p(1,j)=p(2,j);
    end
    
    pn=p;
  end
  
% Algorithm ends

% Plotting 
    % Center-line X - direction Velocity
    figure(1);  
    plot(y,u(round(ny/2),:));
    hold on 
    scatter(y_ghia,u_ghia,'filled')
    title('U-Velocity Distribution (@ y = 0.5)');
    xlabel('y');  
    ylabel('u');  
    axis('square');  
    legend('Obtained','Ghia')
    xlim([0 a]);

    % Center-line Y - Direction velocity
    figure(2);  
    plot(x,v(:,round(ny/2)));
    hold on 
    scatter(x_ghia,v_ghia,'filled')
    title('V-Velocity Distribution(@ x = 0.5)');
    xlabel('x');  
    ylabel('v');  
    axis('square');
    legend('Obtained','Ghia')
    xlim([0 a]);

    % Mesh for stream function and contour plots
    grid on  
    N = 1000;  
    xst = max(x)*rand(N,1);  
    yst = max(y)*rand(N,1);
    [X,Y] = meshgrid(x,y);

    % Stream Function
    figure(3);  
    g=streamline(X,Y,u',v',xst,yst,[0.1, 200]);
    title('Stream Lines');  
    xlabel('x');  
    ylabel('y')
    axis('equal',[0 a 0 a]);    
    set(g,'color','k')

    % Vorticity
    figure(4); 
    contourf(x,y,omega_n1',23,'LineColor','none');
    title('Vorticity Contours');  
    xlabel('x');  
    ylabel('y')  
    axis('equal',[0 a 0 a]);  
    colormap(hot);  
    colorbar('eastoutside');   
    
    % U- Velocity Contours
    figure(5); 
    contourf(x,y,u',23,'LineColor','none');
    title('U-Velocity Contours');  
    xlabel('x');  
    ylabel('y')  
    axis('equal',[0 a 0 a]);  
    colormap(hot);  
    colorbar('eastoutside');   
    
    % V- Velocity Contours
    figure(6); 
    contourf(x,y,v',23,'LineColor','none');
    title('V-Velocity Contours');  
    xlabel('x');  
    ylabel('y')  
    axis('equal',[0 a 0 a]);  
    colormap(hot);  
    colorbar('eastoutside');   
    
    % Pressure
    figure(7);  
    contourf(x,y,pn',50,'LineColor','none');
    title('Pressure Distribution');  
    xlabel('x');  
    ylabel('y')  
    axis('equal',[0 a 0 a]);  
    colormap(jet);  
    colorbar('eastoutside');   
else 
    disp("FTCS is unstable: CFL Criteria is not met")
end

toc