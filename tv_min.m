% adapted from the Matlab code by Benjamin Tremoulheac, available at
% https://www.mathworks.com/matlabcentral/fileexchange/36278-split-bregman-
% method-for-total-variation-denoising?s_tid=srchtitle
%
% the scheme is based on 
% T. Goldstein and S. Osher, “The split bregman method for l1-regularized
% problems,” SIAM Journal on Imaging Sciences, vol. 2, no. 2, pp.
% 323–343, 2009. [Online]. Available: https://doi.org/10.1137/080725891

function u = tv_min(u0, mu)
    n = length(u0);
    uin = u0(:);
    up = uin;
    u = up;
    N = length(up);
    % forward and backward difference matrices, 
    % and the laplacian
    [Df, Db, Lapl] = ders(n);
    % d = [dx;dy], same for b
    d = zeros(2*N,1); b = zeros(2*N,1);
    err = 1; lam = 1; tol = 1e-3; k = 1;
    while err > tol
        up = u;
        % optimizing over u
        [u,~] = cgs(speye(N)-Lapl, uin+lam*Db*(d-b),1e-5,10);
        Dfub = Df*u+b;
        % optimizing over dx, dy via shrinkage
        d = max(abs(Dfub)-mu/lam,0).*sign(Dfub);
        b = Dfub-d;
        err = norm(u-up)/norm(u);
        k = k+1;
    end
    u = reshape(u,n,n);
    fprintf('Mu = %.2f | Number of iterations: %i.\n', mu, k);
end

function [Df, Db, Lapl] = ders(N)
    e = ones(N,1);
    D = spdiags([-e e], [0 1], N, N+1);
    D(:,1) = [];
    Df = [kron(speye(N),D); kron(D,speye(N))];
    Db = Df';
    Lapl = -Db*Df;
end