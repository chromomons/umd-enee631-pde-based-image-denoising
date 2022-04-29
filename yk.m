% implementation by Mansur Shakipov
%
% based on
%
% Y.-L. You and M. Kaveh, “Fourth-order partial differential equations
% for noise removal,” IEEE Transactions on Image Processing,
% vol. 9, no. 10, pp. 1723–1730, Oct. 2000. [Online]. Available:
% http://ieeexplore.ieee.org/document/869184/
function In = yk(I0,K)
    tol = 1e-3; dt = 1/50;
    err = 1; k = 1;
    Io = I0; In = zeros(size(Io));
    % mask
    lpl = [0 1 0; 1 -4 1; 0 1 0];
    while err > tol
        LIo = imfilter(Io,lpl,'symmetric','conv');
        g = gfun(LIo,K).*LIo;
        Lg = imfilter(g,lpl,'symmetric','conv');
        In = Io - dt*Lg;
        err = norm(In(:)-Io(:))/norm(In(:));
        Io = In;
        k = k+1;
    end
    fprintf('Kappa = %.3f | Number of iterations: %i.\n', K, k);
end

function res = gfun(x,K)
    res = K^2./(K^2+x.^2);
end