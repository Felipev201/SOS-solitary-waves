function [ dx,TH,RHO,skx,sky ] = gridfft2( L,N,dimskxy,dimtr )
%gridfft2: Returns several parameters that are going to be needed to
% calculate the fitness function. However, unlike the fitness function, 
% these parameters are going to be calculated only once. Therefore this 
% separate function is added.

%   2*N: the number of points taken for each dimension of the grid
%   L: window length, (in range -L:L)
%   dimskxy: size of the third dimension of frequency tensor
%   dimskr: size of the third dimensiona of space tensor

dx=L/N;     % Delta x (space)
dk=pi/L;    % Delta k (frquency)
vec=-N:N-1;

x=vec*dx;   % space vector
k=vec*dk;   % frequency vector

% generating two dimensional grid
[X,Y]=meshgrid(x,x);
[kx,ky]=meshgrid(k,k);
% Change to cylindrical coordinates
[TH,RHO] = cart2pol(X,Y);

% compute the fourier shift of the frequency
nd=1;
skx = fftshift((1i*kx).^nd);
sky = fftshift((1i*ky).^nd);

% to avoid the use of loops we use a 3D array
TH=repmat(TH,[1,1,dimtr]);  RHO=repmat(RHO,[1,1,dimtr]);
skx=repmat(skx,[1,1,dimskxy]);  sky=repmat(sky,[1,1,dimskxy]);

end