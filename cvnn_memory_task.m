%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% linear cv-NN for designing computations    %%%
%%%                                            %%%
%%% Roberto Budzinski, Alex Busch, Lyle Muller %%%
%%% May 2024                                   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Memory task

%% setup

clearvars; clc;

addpath( './analysis' );
addpath( './graphs' );
addpath( './helper_functions' );

%% building network and parameters

% parameters
dt = 0.001;         %timestep
T1 = 1.0; t1 = 0:dt:T1;
T2 = 3.0; t2 = 0:dt:T2;
T3 = 3.0; t3 = 0:dt:T3;
T4 = 1.0; t4 = 0:dt:T4;
N = 321;            %number of nodes
epsilon = 45;       %coupling strength
phi = 1.55;         %phase-lag
f_mu = 10;          %(Hz) natural frequency
omega = f_mu*2*pi;  %(rad) natural frequency

% adjacency matrix (distance-depedent graph)
alpha = 1.0; %power-law exponent
a = distance_dependent_graph( N, alpha);

% matrix K
K = (epsilon) .* exp(-1i*phi) .* a;

% eigensystem
[v,d] = circulant_eigensystem( K ); %analytical eigensystem

%% Memory task

%number of items in the task
N_i = 8;

%first item
n_i_1 = 2; %item number

%second item
n_i_2 = 6; %item number

%number of nodes in each synchronized cluster
N_m = floor(N/N_i);

%first memory target
cluster_idx = 1+(n_i_1 - 1)*N_m:N_m+(n_i_1 - 1)*N_m; phase_angle = 0.0;
theta_target = 2*pi*( rand(N,1) - 0.5 );
theta_target(cluster_idx) = ones(size(theta_target(cluster_idx)))*phase_angle;

time_1 = T2;

x_target_1 = (rand(N,1)*0.5+2) .* exp(1i*theta_target);

%second memory target
cluster_idx = 1+(n_i_2 - 1)*N_m:N_m+(n_i_2 - 1)*N_m; phase_angle = 0.0;
theta_target = 2*pi*( rand(N,1) - 0.5 );
theta_target(cluster_idx) = ones(size(theta_target(cluster_idx)))*phase_angle;

time_2 = T3;

x_target_2 = (rand(N,1)*0.5+2) .* exp(1i*theta_target);

% initial conditions
x_0_1 = exp(1i * 2*pi*( rand(N,1) - 0.5 ));  %random initial conditions for async behavior

x_0_2 = exp(-1i*omega*time_1) .* expm(-time_1 *K) * x_target_1; %initial conditons for the first memory

x_0_3 = exp(-1i*omega*time_2) .* expm(-time_2 *K) * x_target_2; %initial conditons for the second memory

x_0_4 = exp(1i * 2*pi*( rand(N,1) - 0.5 ));  %random initial conditions for async behavior

t = 0:dt:T1+T2+T3+T4;

x = zeros( length(t), N );

jj_time = 1;
% cv-NN dynamics before first item
for jj = 1:length(t1)
    x(jj_time,:) = exp( 1i * omega * t1(jj) ) .* ( v * diag( exp( diag(d) * t1(jj) ) ) * v' ) * x_0_1 ;
    jj_time = jj_time +1;
end

jj_time = jj_time -1;
% cv-NN dynamics first item
for jj = 1:length(t2)
    x(jj_time,:) = exp( 1i * omega * t2(jj) ) .* ( v * diag( exp( diag(d) * t2(jj) ) ) * v' ) * x_0_2 ;
    jj_time = jj_time +1;
end

jj_time = jj_time -1;
% cv-NN dynamics second item
for jj = 1:length(t3)
    x(jj_time,:) = exp( 1i * omega * t3(jj) ) .* ( v * diag( exp( diag(d) * t3(jj) ) ) * v' ) * x_0_3 ;
    jj_time = jj_time +1;
end

jj_time = jj_time -1;
% cv-NN dynamics after second item
for jj = 1:length(t4)
    x(jj_time,:) = exp( 1i * omega * t4(jj) ) .* ( v * diag( exp( diag(d) * t4(jj) ) ) * v' ) * x_0_4 ;
    jj_time = jj_time +1;
end

%% plot

%fig - spatiotemporal phase CVNN 
fg1 = figure;
imagesc( t, 1:N, angle(x).' );
xlabel( 'time (s)' ); ylabel( 'nodes' ); %title( 'CVNN - phase' );
set( gca, 'fontname', 'arial', 'fontsize', 22, 'linewidth', 2 );
colormap bone
set(gcf,'position',[224         227        1121         524]);
