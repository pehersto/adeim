%%%%
% Simple toy problem to demonstrate the adaptive DEIM (aDEIM) introduced in
% Peherstorfer, B. & Willcox, K. Online Adaptive Model Reduction for Nonlinear Systems via Low-Rank Updates.
% SIAM Journal on Scientific Computing, 37(4):A2123-A2150, SIAM, 2015.
%%%%
% https://github.com/pehersto/adeim
%%%%

clear variables;
close all;

N = 64; % number of grid points in spatial domain
Nt = 10000; % number of time steps
Tstart = 1; % start time
Tend = 3; % end time
dim = 2; % dimension of (a)DEIM space
w = 25; % window size
ms = 5; % number of additional(!) sampling points

% construct static DEIM using SVD
t = linspace(Tstart, Tend, Nt)';
x = linspace(0, 2*pi, N)';
[lS, ~, ~] = svd(nonlin(x, t), 0);
% static DEIM
Ustatic = lS(:, 1:dim);
Pstatic = deim(Ustatic);
% aDEIM
U = Ustatic;
P = Pstatic;

for i=w:Nt
    
    % generate sampling points
    Pu = unique([P; randi(N, ms, 1)]);
    
    % evaluate nonlinear function at sampling points for time steps within
    % the window given by w
    F = nonlin(x(Pu), t(i - w + 1:i)); % window size is w
    % adapt DEIM
    [U, P] = aDEIM(U, P, Pu, F);
    % approximate solution with adaptive DEIM
    solapprox = U*(U(P, :)\nonlin(x(P), t(i)));
    
    % generate a plot every 1000th time step that compares aDEIM with DEIM
    if(mod(i, 1000) == 0)
        figure('visible', 'off');
        plot(nonlin(x, t(i)), '-xr');
        hold on;
        plot(Ustatic*(Ustatic(Pstatic, :)\nonlin(x(Pstatic), t(i))), '-ms');
        hold on;
        plot(solapprox, '-b'); 
        xlabel('spatial domain');
        ylabel('function output');
        legend('nonlinear function (ground truth)', 'DEIM', 'aDEIM');
        saveas(gcf, ['output_t', num2str(i), '.pdf']);
    end
end



