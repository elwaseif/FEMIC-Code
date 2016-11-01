function[ WTW] = calcWTW(dx,dz,dy,wt,alx,aly,alz,als)
% [WTW] = calcWTW(MTX,wt)
% Calculate WTW - the model regularization matrix
% USE: grad, kron3

% Copyright (c) 2007 by the Society of Exploration Geophysicists.
% For more information, go to http://software.seg.org/2007/0001 .
% You must read and accept usage terms at:
% http://software.seg.org/disclaimer.txt before use.
% 
% Revision history:
% Original SEG version by Adam Pidlisecky and Eldad Haber
% Last update, July 2005

nx=length(dx);
ny=length(dy);
nz=length(dz);

[Gx,Gz,Gy] = grad(dx,dy,dz);


%%Create a weighted smallness term 
V = spdiags(mkvc(wt), 0, nx*ny*nz, nx*ny*nz); 

%Assemble the Anisotropic gradient operrator
Gs = [alx*Gx;alz*Gz;aly*Gy];

%Weights certain points more than others
Wt = spdiags(mkvc(wt),0,nx*ny*nz,nx*ny*nz);


%assemble the 3d weighting matrix
WTW = Wt' * ( Gs' * Gs + als * V) * Wt;

