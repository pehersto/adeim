function [ Unew, Pnew ] = aDEIM( Uold, Pold, S, F )
%ADEIM Computes update to DEIM space following Algorithm 1 and 2 in
% Peherstorfer, B. & Willcox, K. Online Adaptive Model Reduction for Nonlinear Systems via Low-Rank Updates.
% SIAM Journal on Scientific Computing, 37(4):A2123-A2150, SIAM, 2015.
% In
%   Uold        ...         old DEIM basis matrix
%   Pold        ...         old DEIM interpolation points matrix
%   S           ...         sampling points
%   F           ...         function values in window at sampling points
% Out
%   Unew        ...         adapted DEIM basis matrix
%   Pnew        ...         adapted DEIM interpolation points matrix
%%%%
% https://github.com/pehersto/adeim
%%%%

% relative tolerance for rank truncation
rankTol = 1e-14;

% compute DEIM coefficients matrix and residual at sampling points
C = Uold(S, :)\F; % DEIM coefficients
res = Uold(S, :)*C - F; % residual at sampling points

% reveal rank of C matrix following Lemma 3.4
[Q, R, E] = qr(C);
I = mean(abs(R), 2)/norm(R, 'fro')^2 > rankTol;
RTilde = R(I, :);
RR = RTilde*RTilde';

% now use RTilde instead of C to solve eigenvalue problem (Lemma 3.5)
% compute update vectors alpha and beta
[rMat, eMat] = eig(RTilde*(res*E)'*(res*E)*RTilde', RR);
[~, maxI] = max(diag(eMat));
normBSquare = rMat(:, maxI)'*RR*rMat(:, maxI);
beta = Q(:, 1:size(RTilde, 1))*rMat(:, maxI);
alpha = -1/normBSquare*res*C'*beta;

% apply update to basis
Unew = Uold;
Unew(S, :) = Unew(S, :) + alpha*beta';

% update DEIM interpolation points
Pnew = updateP(Unew, Uold, Pold);

end

