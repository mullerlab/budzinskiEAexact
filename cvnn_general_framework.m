%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% linear cv-NN for designing computations    %%%
%%%                                            %%%
%%% Roberto Budzinski, Alex Busch, Lyle Muller %%%
%%% May 2024                                   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% setup

clearvars; clc;

addpath( './analysis' );
addpath( './graphs' );
addpath( './helper_functions' );

%% building network and parameters

% parameters
dt = 0.001;         %timestep
T = 6.0;            %total time
t = 0:dt:T;         %array with time
N = 200;            %number of nodes
epsilon = 50;       %coupling strength
phi = 1.55;         %phase-lag
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

%designing target state - chimera
cluster_idx = [51:100];        %position of synchronized cluster
phase_angle = 0.0;             %phase of synchronized cluster

rng(0);
theta_target = 2*pi*( rand(N,1) - 0.5 );
theta_target(cluster_idx) = ones(size(theta_target(cluster_idx)))*phase_angle;

time = T; %time when the target state should appear

% desiging target state in the complex plane
x_target = (rand(N,1)*0.5+2) .* exp(1i*theta_target);

% analytical solution for the initial state 
x_0 = exp(-1i*omega*time) .* expm(-time *K) * x_target; %initial conditons for the target


%% dynamics of cv-NN

x = zeros( length(t), N );

% closed-form solution for cv-NN dynamics
for jj = 1:length(t)
    x(jj,:) = exp( 1i * omega * t(jj) ) .* ( v * diag( exp( diag(d) * t(jj) ) ) * v' ) * x_0 ;
end

%% plot

%fig - spatiotemporal phase CVNN 
fg1 = figure;
imagesc( t, 1:N, angle(x).' );
xlabel( 'time (s)' ); ylabel( 'nodes' );
colormap bone
cb = colorbar();
cb.Label.String = 'phase (rad)';
title('phase dynamics cv-NN');
set( gca, 'fontname', 'arial', 'fontsize', 16, 'linewidth', 2 );
set(gcf,'position',[2   671   823   276]);

%fig - spatiotemporal amplitude CVNN 
fg2 = figure;
imagesc( t, 1:N, abs(x).' );
xlabel( 'time (s)' ); ylabel( 'nodes' );
cb = colorbar();
cb.Label.String = 'amplitude';
colormap sky
title('amplitude dynamics cv-NN');
set( gca, 'fontname', 'arial', 'fontsize', 16, 'linewidth', 2 );
set(gcf,'position',[2   317   823   276]);

%target pattern (complex plane)
fg3 = figure;
hold on
scatter(real(x(end,:)), imag(x(end,:)), 150, 'filled');
scatter(real(x_target), imag(x_target), 30, 'filled');
xlabel('real');ylabel('imag');
legend({'simulation', 'target'});
title('final state');
set( gca, 'fontname', 'arial', 'fontsize', 16, 'linewidth', 2 );
set(gcf,'position',[1     0   402   261]);

%target pattern (phase)
fg4 = figure;
hold on
scatter(1:N, angle(x(end,:)), 150, 'filled');
scatter(1:N, angle(x_target), 30, 'filled');
xlabel('nodes');ylabel('phase (rad)');
xlim([0 N+1]); ylim([-3.5 3.5]);
legend({'simulation', 'target'});
title('final state (phase)');
set( gca, 'fontname', 'arial', 'fontsize', 16, 'linewidth', 2 );
set(gcf,'position',[404     0   402   261]);

%initial state (complex plane)
fg5 = figure;
scatter(real(x_0), imag(x_0), 100, 'filled');
xlabel('real');ylabel('imag');
title('initial state')
set( gca, 'fontname', 'arial', 'fontsize', 16, 'linewidth', 2 );
set(gcf,'position',[806     0   402   261]);

%initial state (phase)
fg6 = figure;
scatter(1:N, angle(x_0), 100, 'filled');
xlabel('nodes');ylabel('phase (rad)');
xlim([0 N+1]); ylim([-3.5 3.5]);
title('initial state (phase)');
set( gca, 'fontname', 'arial', 'fontsize', 16, 'linewidth', 2 );
set(gcf,'position',[1208     0   402   261]);

%decoder
N_output = 4;             %number of decoders
N_s = floor(N/N_output);
threshold = 0.7;          %used in the Heaviside function

fg7 = figure;
set(gcf,'position',[827   337   853   607]);
for ii = 1:N_output
    rr = order_parameter(angle(x(:,1+(ii-1)*N_s:N_s+(ii-1)*N_s)),N_s);  %order parameter of each patch in the cv-NN
    oo = nan(size(rr));
    oo( rr > threshold) = 1; oo( rr<=threshold) = 0; 
    subplot(2,2,ii);
    hold on
    h2 = plot(t,rr, 'linewidth', 2);
    h1 = plot(t,oo, 'linewidth', 3, 'color', [0.5 0.5 0.5 0.2]);
    xlim([0 T]); ylim([-0.05 1.05]);
    xlabel('time (s)');
    box off
    set( gca, 'fontname', 'arial', 'fontsize', 14, 'linewidth', 2 );  
    uistack( h1, 'bottom' );
    str1 = sprintf('$R_{%d}[x(t)]$',ii);
    str2 = sprintf('$o_{%d}(t)$',ii);
    legend({str2, str1}, 'interpreter', 'latex', 'location', 'northwest');
end

