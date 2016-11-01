function[Gx,Gy,Gz] = grad(dx,dy,dz)
% [G] = grad(dx,dy,dz)
%Creates the 3D finite volume gradient operator
%operator is set up to handle variable grid discretization
%dx,dy,dz are vectors containing the cell widths in the x y and z
%directions, respectively
%
% Copyright (c) 2007 by the Society of Exploration Geophysicists.
% For more information, go to http://software.seg.org/2007/0001 .
% You must read and accept usage terms at:
% http://software.seg.org/disclaimer.txt before use.
% 
% Revision history:
% Original SEG version by Adam Pidlisecky and Eldad Haber
% Last update, July 2006

dx = shiftdim(dx);
dy = shiftdim(dy);
dz = shiftdim(dz);

Nx = length(dx)-2;
Ny = length(dy)-2;
Nz = length(dz)-2;

%dx = [dx(1);dx;dx(end)];
%dy = [dy(1);dy;dy(end)];
%dz = [dz(1);dz;dz(end)];


% Number the phi grid 
np = (Nx+2)*(Ny+2)*(Nz+2);
GRDp = reshape(1:1:np,Nx+2,Ny+2,Nz+2);


% Number the Ax grid
nax = (Nx+1)*(Ny+2)*(Nz+2); 
GRDax = reshape(1:1:nax, (Nx+1),(Ny+2),(Nz+2));

% Number the Ay grid
nay = (Nx+2)*(Ny+1)*(Nz+2); 
GRDay = reshape(1:1:nay, (Nx+2),(Ny+1),(Nz+2));

% Number the Az grid
naz = (Nx+2)*(Ny+2)*(Nz+1); 
GRDaz = reshape(1:1:naz, (Nx+2),(Ny+2),(Nz+1));


%%%%   Generate d/dx  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lx = []; jx = []; kx = [];

% Generate grid
ex = ones(Nx+2,1);
ey = ones(Ny+2,1);
ez = ones(Nz+2,1);

Dx = kron3(dx,ey,ez);
% Entries (l,j,k)

lx = mkvc(GRDax);
jx = mkvc(GRDp(1:end-1,:,:));
kx = mkvc(-2./(Dx(1:end-1,:,:) + Dx(2:end,:,:)));

% Entries (l+1,j,k)

lx = [lx; lx];
jx = [jx;mkvc(GRDp(2:end,:,:))];
kx = [kx;-kx];


%%%%   Generate d/dy  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ly = []; jy = []; ky = [];

Dy = kron3(ex,dy,ez);

% Entries (l,j,k)

ly = mkvc(GRDay);
jy = mkvc(GRDp(:,1:end-1,:));
ky = mkvc(-2./(Dy(:,1:end-1,:) + Dy(:,2:end,:)));

% Entries (l,j+1,k)

ly = [ly; ly];
jy = [jy;mkvc(GRDp(:,2:end,:))];
ky = [ky;-ky];


%%%%   Generate d/dz  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lz = []; jz = []; kz = [];

Dz = kron3(ex,ey,dz);

% Entries (l,j,k)

lz = mkvc(GRDaz);
jz = mkvc(GRDp(:,:,1:end-1));
kz = mkvc(-2./(Dz(:,:,1:end-1) + Dz(:,:,2:end)));;

% Entries (l,j,k+1)

lz = [lz; lz];
jz = [jz;mkvc(GRDp(:,:,2:end))];
kz = [kz; -kz];

Gx = sparse(lx,jx,kx,nax,np);
Gy = sparse(ly,jy,ky,nay,np); 
Gz = sparse(lz,jz,kz,naz,np);


G = [Gx;Gy;Gz];

