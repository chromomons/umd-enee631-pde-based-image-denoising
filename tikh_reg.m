% implementation by Mansur Shakipov
function In = tikh_reg(I,lambda)
    In = I; Io = In;
    mask = [0 1 0; 1 -4 1; 0 1 0];
    dt = .1; err = 1; k = 1;
    maxit = 5e3; tol = 1e-3;
    while err > tol && k < maxit
        In = Io + dt*(imfilter(Io,mask,'symmetric','conv') - lambda*(Io-I));
        err = norm(In(:)-Io(:))/norm(In(:));
        Io = In; k = k+1;
    end
    fprintf('Lambda = %.3f | Number of iterations: %i.\n', lambda, k);
end