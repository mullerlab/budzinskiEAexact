%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% linear cv-NN for designing computations    %%%
%%%                                            %%%
%%% Roberto Budzinski, Alex Busch, Lyle Muller %%%
%%% May 2024                                   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% XOR gate

%% setup

clearvars; clc;

addpath( './analysis' );
addpath( './graphs' );
addpath( './helper_functions' );

%% building network and parameters

% parameters
dt = 0.001;         %timestep
T = 3.0;            %total time
t = 0:dt:T;         %array with time
N = 201;            %number of nodes
epsilon = 50;       %coupling strength
phi = 1.56;         %phase-lag
f_mu = 10;          %(Hz) natural frequency
omega = f_mu*2*pi;  %(rad) natural frequency

% adjacency matrix (distance-depedent graph)
alpha = 1.0; %power-law exponent
a = distance_dependent_graph( N, alpha );

% matrix K
K = (epsilon) .* exp(-1i*phi) .* a;

% eigensystem
[v,d] = circulant_eigensystem( K ); %analytical eigensystem

%% initial state and input to the network

N_m = floor(N/4);
cluster_idx = [51:150];        %position of synchronized cluster

%defining input x
phase_angle = -1.5;
theta_target = 2*pi*( rand(N,1) - 0.5 );
theta_target(cluster_idx) = ones(size(theta_target(cluster_idx)))*phase_angle;
time = T;
x_target = (rand(N,1)*2+1.5) .* exp(1i*theta_target);
%input for x
x_0_x = exp(-1i*omega*time) .* expm(-time *K) * x_target; %initial conditons for the target

%defining input y
phase_angle = 1.5;
theta_target = 2*pi*( rand(N,1) - 0.5 );
theta_target(cluster_idx) = ones(size(theta_target(cluster_idx)))*phase_angle;
time = T;
x_target = (rand(N,1)*2+1.5) .* exp(1i*theta_target);
%input for y
x_0_y = exp(-1i*omega*time) .* expm(-time *K) * x_target; %initial conditons for the target


%% realizing XOR gate for x=0; y=0

x_0 = exp(1i * 2*pi*( rand(N,1) - 0.5 ));

x = zeros( length(t), N );

% cv-NN dynamics
for jj = 1:length(t)
    x(jj,:) = exp( 1i * omega * t(jj) ) .* ( v * diag( exp( diag(d) * t(jj) ) ) * v' ) * x_0 ;
end

%fig - spatiotemporal phase cv-NN
fg1 = figure;
imagesc( t, 1:N, angle(x).' );
xlabel( 'time (s)' ); ylabel( 'nodes' ); %title( 'CVNN - phase' );
set( gca, 'fontname', 'arial', 'fontsize', 16, 'linewidth', 2 );
colormap bone
set(gcf,'position',[1   755   783   192])

%fig - final state cv-NN
fg2 = figure;
scatter(1:N,angle(x(end,:)),'o','filled','black')
xlabel('nodes');ylabel('phase');
xlim([0 N+1]);ylim([-4 4]);
set( gca, 'fontname', 'arial', 'fontsize', 16, 'linewidth', 2 );
set(gcf,'position',[785   763   405   184])

%% realizing XOR gate for x=1; y=0

x_0 = x_0_x;

x = zeros( length(t), N );

% cv-NN dynamics
for jj = 1:length(t)
    x(jj,:) = exp( 1i * omega * t(jj) ) .* ( v * diag( exp( diag(d) * t(jj) ) ) * v' ) * x_0 ;
end

%fig - spatiotemporal phase cv-NN
fg3 = figure;
imagesc( t, 1:N, angle(x).' );
xlabel( 'time (s)' ); ylabel( 'nodes' ); %title( 'CVNN - phase' );
set( gca, 'fontname', 'arial', 'fontsize', 16, 'linewidth', 2 );
colormap bone
set(gcf,'position',[1   491   783   192])

%fig - final state cv-NN
fg4 = figure;
scatter(1:N,angle(x(end,:)),'o','filled','black')
xlabel('nodes');ylabel('phase');
xlim([0 N+1]);ylim([-4 4]);
set( gca, 'fontname', 'arial', 'fontsize', 16, 'linewidth', 2 );
set(gcf,'position',[785   499   405   184])

%% realizing XOR gate for x=0; y=1

x_0 = x_0_y;

x = zeros( length(t), N );

% cv-NN dynamics
for jj = 1:length(t)
    x(jj,:) = exp( 1i * omega * t(jj) ) .* ( v * diag( exp( diag(d) * t(jj) ) ) * v' ) * x_0 ;
end

%fig - spatiotemporal phase cv-NN
fg5 = figure;
imagesc( t, 1:N, angle(x).' );
xlabel( 'time (s)' ); ylabel( 'nodes' ); %title( 'CVNN - phase' );
set( gca, 'fontname', 'arial', 'fontsize', 16, 'linewidth', 2 );
colormap bone
set(gcf,'position',[1   245   783   192])

%fig - final state cv-NN
fg6 = figure;
scatter(1:N,angle(x(end,:)),'o','filled','black')
xlabel('nodes');ylabel('phase');
xlim([0 N+1]);ylim([-4 4]);
set( gca, 'fontname', 'arial', 'fontsize', 16, 'linewidth', 2 );
set(gcf,'position',[785   253   405   184])

%% realizing XOR gate for x=1; y=1

x_0 = x_0_x + x_0_y;

x = zeros( length(t), N );

% cv-NN dynamics
for jj = 1:length(t)
    x(jj,:) = exp( 1i * omega * t(jj) ) .* ( v * diag( exp( diag(d) * t(jj) ) ) * v' ) * x_0 ;
end

%fig - spatiotemporal phase cv-NN
fg7 = figure;
imagesc( t, 1:N, angle(x).' );
xlabel( 'time (s)' ); ylabel( 'nodes' ); %title( 'CVNN - phase' );
set( gca, 'fontname', 'arial', 'fontsize', 16, 'linewidth', 2 );
colormap bone
set(gcf,'position',[1    1   783   192])

%fig - final state cv-NN
fg8 = figure;
scatter(1:N,angle(x(end,:)),'o','filled','black')
xlabel('nodes');ylabel('phase');
xlim([0 N+1]);ylim([-4 4]);
set( gca, 'fontname', 'arial', 'fontsize', 16, 'linewidth', 2 );
set(gcf,'position',[785   1   405   184])
