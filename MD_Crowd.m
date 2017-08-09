function y_out = MD_Crowd()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% July 2017, Orit Peleg, opeleg@seas.harvard.edu
% Edited August 2017, Ethan Hobbs and Rebecca Wayne, ehobbs@carthage.edu
% and Rebecca_Wayne@student.uml.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
dbstop if error
global stage_pos N 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize parameters 
% N LJ particles in a box, fixed temperature

N = 300; % # particles,
T = 20; %temperature
dt = 0.0005; %integration time step
steps = 10001; %time steps
cutoff = 2.5;
wallcutoff = 2.5;
L = 35;
min_sep = 1.122;
stage_pos = [0, L/2];
print_interval = 50;
target1 = [15,15];
target2 = [0,12.5];
z = 0;
smart_x = [0, -15];
smart_F = zeros(1, 2);
mag_F = 0;
Forces_output = zeros(steps,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize position velocities and forces  
[x,v] =useful_initial_configuration(T,L-cutoff,min_sep,smart_x); %init. coordinates x,
% F_particles = forces(N,x,cutoff); %velocities v and forces F
% F_wall= forces_wall (N,x,L,wallcutoff);
% F_stage = forces_stage (N,x,L,L);
% F = F_particles + F_wall + F_stage;

ax = gca;
ax.NextPlot = 'replaceChildren';
clear F
F(floor(steps/print_interval)) = struct('cdata',[],'colormap',[]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% main loop 
for step_i=1:steps, %molecular dynamics loop
    %[x,v,F]=velocity_verlet(N,x,v,F,dt,L,cutoff,wallcutoff); %propagate trajectory
    %v = temperature_control(N,v,T); %velocity rescaling
    [x, smart_F_particles] = steepest_descent (N,x,smart_x,dt,L,cutoff,wallcutoff, T); 
    [smart_x, z] = smart_steepest_descent (N, x, smart_x, dt, L, cutoff, wallcutoff, target1, target2, z, T);
    if mod(step_i-1,print_interval)==0
        mytitle = ['step=',num2str(step_i), ' N=',num2str(N), ' L=',num2str(L)] ;
        visualize_particles (N,x,smart_x,L,mytitle);
        F(floor(step_i/print_interval)+1) = getframe;
    end
    Forces_output(step_i) = norm(smart_F_particles);
end

%Saving the graphs as MP4
movie_name = ['movie_','T_',num2str(T),'_N_',num2str(N),'_steps_',num2str(steps),'_target1_[',num2str(target1),']_target2_[',num2str(target2),']']
v = VideoWriter([movie_name,'.mp4'],'MPEG-4'); 
open(v);
writeVideo(v,F);
close(v);

%Saving final plot of Force vs. Time
graph_name = ['Force_vs_Time_','T_',num2str(T),'_N_',num2str(N),'_steps_',num2str(steps),'_target1_[',num2str(target1),']_target2_[',num2str(target2),']'] 
figure(2)
plot(Forces_output, 'b-')
saveas(gcf,graph_name,'png')

%Outputting the Force on the Smart Agent every timestep
csv_name = ['Force_','T_',num2str(T),'_N_',num2str(N),'_steps_',num2str(steps),'_target1_[',num2str(target1),']_target2_[',num2str(target2),']'] 
csvwrite(csv_name, Forces_output)

y_out = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Position update for normal agents
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x, smart_F_particles ]= steepest_descent (N,x,smart_x,dt,L,cutoff,wallcutoff,T) %velocity Verlet integrator
F_particles = forces(N,x,cutoff);
[smart_F_particles, interact_F] = smart_forces(N, x, smart_x, cutoff); 
F_wall= forces_wall (N,x,L,wallcutoff);
F_stage = forces_stage (N,x,L,L);
noise = T * normrnd(0,1,size(x));
F = F_particles + F_wall + F_stage + interact_F + noise;
x = x + dt*F;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Position update for the Smart Agent
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [smart_x, z] = smart_steepest_descent (N,x,smart_x,dt,L,cutoff,wallcutoff, target1, target2, z, T) %velocity Verlet integrator
[smart_F_particles, interact_F] = smart_forces(N,x,smart_x,cutoff);
smart_F_wall= smart_forces_wall (smart_x,L,wallcutoff);
[smart_dir, z] = smart_direction(smart_x, target1, target2, z);
smart_noise = T * normrnd(0,1,1);
smart_F = smart_F_particles + smart_F_wall + smart_dir(1,:) + smart_noise;
smart_x = smart_x + dt*smart_F;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculating the direction and force towards the desired targets
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [smart_dir, z] = smart_direction(smart_x, target1, target2, z) %the smart agents smart direction
a = (target1 - smart_x)/norm(target1 - smart_x);
b = (target2 - smart_x)/norm(target2 - smart_x);
c = [0,0];

if (z == 1)
    smart_dir = 40 * [b; c];
elseif (z == 2)
    smart_dir = c;
else
    smart_dir = 40 * [a; b; c];
end

if abs(smart_x - target1) <= .005
    smart_dir(1,:) = [];
    z = z + 1;
end

if abs(smart_x - target2) <= .01
    if smart_dir(1,:) ~= c;
        smart_dir(1,:) = [];
        z = z + 1;
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Unused but possible implementation of Velocity Verlet integration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x,v,F]= velocity_verlet (N,x,v,F,dt,L,cutoff,wallcutoff) %velocity Verlet integrator
x = x + dt*v + dt^2/2*F; %new x with boundary conditions
v = v + dt/2*F;
F_particles = forces(N,x,cutoff);
F_wall= forces_wall (N,x,L,wallcutoff);
F_stage = forces_stage (N,x,L,L);
F = F_particles + F_wall + F_stage;
v = v + dt/2*F;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Unused Temperature control for Velocity Verlet interation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function v = temperature_control (N,v,T) %rescaling velocities according ?wanted? temperature
T_measured=sum(v(:).^2)/(2*N); v=v*sqrt(T/T_measured);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Caclulates the interactions between all naive agents
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ip,pair,connector]= all_interactions (N,x,cutoff) %obtain interacting pairs
ip=0; connector=zeros(1,2); pair=zeros(1,2);
for i=1:N-1,
    for j=i+1:N,
        distance = (x(j,:)-x(i,:));
        if norm(distance) < cutoff, %only interacting pairs (cutoff+shell for neighbor lists)
            ip = ip + 1; %interaction pair counter
            pair(ip,:) = [i j]; %particle numbers (i,j) belonging to pair (ip)
            connector(ip,:) = distance;
        end; %connecting vector x j-x i for pair (i,j)
    end;
end; %end both ?for? loops
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculates the interaction pairs for the smart agent
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ip,pair,connector]= smart_interactions (N,x,smart_x,cutoff) %obtain interacting pairs
ip=0; connector=zeros(1,2); pair=zeros(1);
for i=1:N,
    distance = (x(i,:)-smart_x);
    if norm(distance) < cutoff, %only interacting pairs (cutoff+shell for neighbor lists)
        ip = ip + 1; %interaction pair counter
        pair(ip) = i; %particle numbers (i,j) belonging to pair (ip)
        connector(ip,:) = distance;
    end; %connecting vector x j-x i for pair (i,j)
end; %end both ?for? loops
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculates all forces from interactions between the normal agents
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function F= forces (N,x,cutoff) %clear forces F, then calculate them using ..
F=zeros(N,2); 
[no,pair,connector]=all_interactions(N,x,cutoff); %interacting pairs

for i=1:no, FORCE=force_LJ(connector(i,:));
    F(pair(i,1),:)=F(pair(i,1),:)-FORCE; 
    F(pair(i,2),:)=F(pair(i,2),:)+FORCE; %actio=reactio;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculates all forces from interactions between the Smart Agent and the
%Noraml Agents
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [smart_F, F]= smart_forces (N,x,smart_x,cutoff) %clear forces F, then calculate them using ..
smart_F = zeros(1,2);
F = zeros(N,2);

[no,pair,connector]=smart_interactions(N,x,smart_x,cutoff); %interacting pairs

for i=1:no, FORCE=force_LJ(connector(i,:));
    smart_F = smart_F - FORCE; 
    F(pair(i),:)=F(pair(i),:) + FORCE; %actio=reactio;
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Creates the interaction pairs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ip,pair,connector]= all_interactions_stage (N,x,cutoff) %obtain interacting pairs
global stage_pos
ip=0; connector=zeros(1,2); pair=zeros(1,2);
for i=1:N,
    distance1 = x(i,:)-stage_pos;
    distance = distance1;
    %(true) connecting vector x j-x i
    if norm(distance) < cutoff, %only interacting pairs (cutoff+shell for neighbor lists)
        ip = ip + 1; %interaction pair counter
        pair(ip,:) = [i NaN]; %particle numbers (i,j) belonging to pair (ip)
        connector(ip,:) = distance;
    end; %connecting vector x j-x i for pair (i,j)
    
end; %end both for loops
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculates the pull towards the stage for the normal particles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function F_stage= forces_stage (N,x,L,cutoff) %clear forces F, then calculate them using ..
[no,pair,connector]=all_interactions_stage(N,x,cutoff); %interacting pairs
F=zeros(N,2);
for i=1:no, FORCE=force_spring_stage(connector(i,:));
    F(pair(i,1),:)=F(pair(i,1),:)-FORCE;
end
F_stage = F;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Finds the Interaction pairs with the Walls
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ip,pair,connector]= all_interactions_wall (N,x,L,cutoff) %obtain interacting pairs
ip=0; connector=zeros(1,2); pair=zeros(1,2);
for i=1:N,
    distance1 = x(i,:)-[L/2,x(i,2)];
    distance2 = x(i,:)-[-L/2,x(i,2)];
    distance3 = x(i,:)-[x(i,1),L/2];
    distance4 = x(i,:)-[x(i,1),-L/2];
    distance_array = {distance1,distance2,distance3,distance4};
    
    for distance_i = 1:length(distance_array)
        distance = distance_array{distance_i};
        %(true) connecting vector x j-x i
        if norm(distance) < cutoff, %only interacting pairs (cutoff+shell for neighbor lists)
            ip = ip + 1; %interaction pair counter
            pair(ip,:) = [i NaN]; %particle numbers (i,j) belonging to pair (ip)
            connector(ip,:) = distance;
        end; %connecting vector x j-x i for pair (i,j)
    end
    
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate the Smart Agents Interaction with the Wall
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ip,pair,connector]= smart_interactions_wall (smart_x,L,cutoff) %obtain interacting pairs
ip=0; connector=zeros(1,2); pair=zeros(1,2);

distance1 = smart_x-[L/2,smart_x(:,2)];
distance2 = smart_x-[-L/2,smart_x(:,2)];
distance3 = smart_x-[smart_x(:,1),L/2];
distance4 = smart_x-[smart_x(:,1),-L/2];
distance_array = {distance1,distance2,distance3,distance4};
    
    for distance_i = 1:length(distance_array)
        distance = distance_array{distance_i};
        %(true) connecting vector x j-x i
        if norm(distance) < cutoff, %only interacting pairs (cutoff+shell for neighbor lists)
            ip = ip + 1; %interaction pair counter
            pair(ip) = [1]; %particle numbers (i,j) belonging to pair (ip)
            connector(ip,:) = distance;
        end; %connecting vector x j-x i for pair (i,j)
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate the forces on the Normal Agents from the Wall
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function F_wall= forces_wall (N,x,L,cutoff) %clear forces F, then calculate them using ..
[no,pair,connector]=all_interactions_wall(N,x,L,cutoff); %interacting pairs
F=zeros(N,2);
for i=1:no, FORCE=force_LJ_wall(connector(i,:));
    F(pair(i,1),:)=F(pair(i,1),:)-FORCE;
end
F_wall = F;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculates the forces on the Smart Agents from the Wall
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function F_wall= smart_forces_wall (smart_x,L,cutoff) %clear forces F, then calculate them using ..
[no,pair,connector]=smart_interactions_wall(smart_x,L,cutoff); %interacting pairs
F=zeros(1,2);
for i=1:no, FORCE=force_LJ_wall(connector(i,:));
    F(pair(i),:)=F(pair(i),:)-FORCE;
end
F_wall = F;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Lennard-Jones Force calculation - Particle to particle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function force_LJ= force_LJ (r_vector); r=norm(r_vector); %two-body force
force_LJ = 24*(2*r.^(-14)-r^(-8)) * r_vector; %here: Lennard?Jones
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Lennard-Jones Force calculation - Particle to wall
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function force_LJ_wall= force_LJ_wall (r_vector); r=norm(r_vector); %two-body force
if r>0
    force_LJ_wall = 24*(-r^(-8)) * r_vector; %here: Lennard?Jones
else
    force_LJ_wall = 0;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Spring force to the Stage - Pull towards the stage
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function force_spring_stage =  force_spring_stage (r_vector); r=norm(r_vector); %two-body force
if r>0
    force_spring_stage = 0.5*(r_vector./r); 
else
    force_spring_stage = 0;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Ploting mechanism for each iteration 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function visualize_particles (N,x,smart_x,L,mytitle) %visualize N spheres, radius r, at positions x
global stage_pos
scatter3(stage_pos(1),stage_pos(2),0,500,'filled','blue','square');  hold on;
scatter3(smart_x(1), smart_x(2), 0, 40, 'filled', 'cyan'); hold on;
for i=1:N, scatter3(x(i,1),x(i,2),0,40,'filled','black'); view(2); hold on; end; 
xlim([-L,L]/2); ylim([-L,L]/2);
axis square; hold off; 
title(mytitle); pause(0.001); %customize axes and title, pause
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Initializes the Particles in an already Crowd formation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x,v]= useful_initial_configuration (T,L,min_sep,smart_x)
global N

tentative_N = N*10;

particles_per_edge = ceil((tentative_N^(1/2)));
initial_L = particles_per_edge*min_sep;

v_indiv = pi*((min_sep/2)^2);

curr_R = sqrt(2*N*v_indiv/pi);

[X,Y] = meshgrid(0:(particles_per_edge-1),0:-1:-(particles_per_edge-1));
even_rows = find(mod(1:particles_per_edge-1,2)==0);
X(even_rows,:)=X(even_rows,:)+(min_sep/2);
i=1; counter = 1;
while counter<tentative_N
    curr_x = min_sep*X(counter) - (initial_L/2) + (min_sep/2);
    curr_y = L/2 -(min_sep/2) + (sqrt(3)/2)*min_sep*Y(counter);
    if ((curr_x^2 + (curr_y-(L/2))^2) < curr_R^2)
        x(i,1) =curr_x;
        x(i,2) =curr_y;
        i = i+1;
    end
    counter = counter+1;
end
N=i-1;
visualize_particles (N,x,smart_x,L,'');

v=temperature_control(N,rand(N,2)-0.5,T); %uniformely distributed random v?s
vcm=sum(v)/N; for k=1:2, v(:,k)=v(:,k)-vcm(k); end; %ensure center of mass vcm=0
end


