function [ P ] = updateP( Unew, Uold, P )
%UPDATEP Updates DEIM points following Algorithm 2 in
% Peherstorfer, B. & Willcox, K. Online Adaptive Model Reduction for Nonlinear Systems via Low-Rank Updates.
% SIAM Journal on Scientific Computing, 37(4):A2123-A2150, SIAM, 2015.
% In
%   Unew    ... new DEIM basis
%   Uold    ... DEIM basis before update
%   Z       ... old DEIM points
% Out
%   Z       ... new DEIM points
%%%%
% https://github.com/pehersto/adeim
%%%%

% find basis vector which changed most w.r.t. dot product
[~, minI] = min(abs(diag(Unew'*Uold))./(sqrt(diag(Unew'*Unew)).*sqrt(diag(Uold'*Uold))));

% remove maxI-th vector from Unew
Utmp = Unew;
Utmp(:, minI) = [];
% remove maxI-th point from P
Ptmp = P;
Ptmp(minI) = [];

% find new point following standard DEIM
c = Utmp(Ptmp, :)\Unew(Ptmp, minI);
r = abs(Unew(:, minI) - Utmp*c); % residual
[~, I] = sort(r, 'descend');
curI = 1;
while(~isempty(find(Ptmp == I(curI), 1))) % check if point is already interp
    curI = curI + 1;
end
P(minI) = I(curI); % add first point that is not already an interp point

end

