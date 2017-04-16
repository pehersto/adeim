function [ y ] = nonlin( x, t )
%NONLIN Toy nonlinear function
% In
%   x       ...     spatial coordinate
%   t       ...     time
% Out
%   y       ...     function value
%%% Used to demonstrate adaptive DEIM (aDEIM) introduced in
% Peherstorfer, B. & Willcox, K. Online Adaptive Model Reduction for Nonlinear Systems via Low-Rank Updates.
% SIAM Journal on Scientific Computing, 37(4):A2123-A2150, SIAM, 2015.
%%%%
% https://github.com/pehersto/adeim
%%%%

y = zeros(length(x), length(t));
for i=1:length(t)
    y(:, i) = (sin(t(i)*x) + sin(t(i)*pi/2*x) + sin(t(i)*pi*x))*t(i);
end

end

