function [ z ] = fitness_sat( lambda,m,dx,TH,RHO,skx,sky,xt,dn )
% Objective function for parameters generating the self trapped beams
% in a saturated media.
% The ansatz that is used is a Gaussian beam vortex.
% The function can be changed to analyze other ansatz functions or
% different values of saturation.

%   lambda, and m:  parameters from the PDE
%   dx: the differential of the space variable
%   TH, RHO: matrix of the space in cilindrical coordinates
%   skx,sky: fftshift vectors of the space frequency
%   xt: optimization variables
%   dn: the differential of the optimized variables
s = 0.05;   % saturation
A = xt(1);
b = xt(2);

dna=dn(1);
dnb=dn(2);

% avoiding the use of  a cycle so we have a 3D matrix to compute U(n) and U(n+dn)
dim=size(TH,1);
mat1=ones(dim,dim);
b1=b^2*mat1;     b2=(b+dnb)^2*mat1;
% cat is faster than repmat
bsq = cat(3,b1,b2);

% --------------------- Compute the derivative -------------------------
% Proposed function
u= (A)*exp(-RHO.^2./bsq).*exp(1i*m*TH).*RHO.^m;
u(:,:,3) = ((A+dna)/A)*u(:,:,1);

% Derivative of the function calculated by Fourier soectral method
ff = fft2(u);
dux = ifft2(skx.*ff);
duy = ifft2(sky.*ff);

% it is faster to calculate the conjugate than to calculate the absolute
% squared
sqdx = dux.*conj(dux);
sqdy = duy.*conj(duy);
squ = u.*conj(u);
logsqu = (log(s*squ + 1) - s*squ)./(s^2);

%------------------------- Objective function -----------------------------
% Substitute in the functional
G = (lambda*squ + sqdx + sqdy + logsqu);

% integrate the functional
f1 = sum(sum(G(:,:,1)))*dx^2;
f2 = sum(sum(G(:,:,2)))*dx^2;
f3 = sum(sum(G(:,:,3)))*dx^2;

% calculate the fitness
za2 = ((f3 - f1)./dna).^2;
zb2 = ((f2 - f1)./dnb).^2;
z = sqrt(za2 + zb2); 
end