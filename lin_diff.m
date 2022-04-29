% implementation by Mansur Shakipov
function I = lin_diff(I,sigma_f)
    [Nt,lambda] = lin_dif_params(I,sigma_f/length(I));
    mask = [0 1 0; 1 -4 1; 0 1 0];
    for t=2:Nt
        I = I + lambda*imfilter(I,mask,'symmetric','conv');
    end
end

function [Nt,lambda] = lin_dif_params(I,sigma_p)
    Tf = 1/2*sigma_p^2;
    [Nx,Ny] = size(I); hx = 1/(Nx-1); hy = 1/(Ny-1);
    Nt = ceil(4*Tf/(hx*hy));
    lambda = ht/(hx*hy);
end