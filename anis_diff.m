% implementation by Mansur Shakipov
% 
% adapted from the Matlab code by Daniel Lopes, available at
% https://www.mathworks.com/matlabcentral/fileexchange/14995-anisotropic-
% diffusion-perona-malik

% the scheme is based on
% P. Perona and J. Malik, “Scale-space and edge detection using
% anisotropic diffusion,” IEEE Transactions on Pattern Analysis and
% Machine Intelligence, vol. 12, no. 7, pp. 629–639, Jul. 1990. [Online].
% Available: http://ieeexplore.ieee.org/document/56205/

% some of the other kernels are due 

% Charbonnier, P., Blanc-Féraud, L., Aubert, G., Barlaud, M.: Deterministic 
% edge-preserving regularization in computed imaging. IEEE Trans. Image 
% Process. 6(2), 298–311 (1997)

% Zhichang Guo, Jiebao Sun, Dazhi Zhang and Boying Wu, Adaptive Perona–
% Malik Model Based on the Variable Exponent for Image Denoising, IEEE
% Transactions on Image Processing

% Joachim Weicket, Coherence Enhancing diffusion Filtering, International 
% journal of Computer Vision, 1999, 31(2-3): 111-127

function In = anis_diff(I0, K, tp, logs)
    tol = 1e-3; dt = 0.1;
    err = 1; k = 1;
    Io = I0; In = zeros(size(Io));
    % masks 
    N = 4; masks = get_masks();
    grads = zeros([size(Io) N]); cs = zeros(size(grads));
    while err > tol
        sm = 0.0;
        for i=1:N
            grads(:,:,i) = imfilter(Io,masks(:,:,i),'symmetric','conv');
            cs(:,:,i) = pm(grads(:,:,i),K,tp);
            sm = sm+cs(:,:,i).*grads(:,:,i);
        end
        In = Io + dt*sm;
        err = norm(In(:)-Io(:))/norm(In(:));
        Io = In;
        k = k+1;
    end
    if logs
        fprintf('Kappa = %.3f | Number of iterations: %i.\n', K, k);
    end
end

function res = pm(x,K,tp)
    if strcmp(tp,'ch')
        res = 1./((1+(x/K).^2).^0.5);
    elseif strcmp(tp,'pm2')
        res = exp(-(x/K).^2);
    elseif strcmp(tp,'zh')
        alp = 2-2./(1+(x/K).^2);
        res = 1./(1+abs((x/K)).^alp);
    elseif strcmp(tp,'wi')
        res = 1-exp(-3.31488*(abs(K./x)).^8);
    else 
        res = 1./(1+(x/K).^2);
    end
end

function msk = get_masks()
    % N S E W
    msk = zeros(3,3,4);
    msk(:,:,1) = [0 1 0; 0 -1 0; 0 0 0];
    msk(:,:,2) = [0 0 0; 0 -1 0; 0 1 0];
    msk(:,:,3) = [0 0 0; 0 -1 1; 0 0 0];
    msk(:,:,4) = [0 0 0; 1 -1 0; 0 0 0];
end