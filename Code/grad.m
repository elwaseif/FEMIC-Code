function[Gx,Gz] = grad(dx,dz)
% [Gx, Gz] = grad(dx,dy,dz)
% code modified from Pidlisecky et al., 2007
%Reference: Pidlisecky, A., Haber, E., and Knight, R., 2007, RESINVM3D: 
% A 3D resistivity inversion package: GEOPHYSICS,VOL. 72, NO. 2
% MARCH-APRIL 2007; P. H1–H10,

dx = shiftdim(dx);
dz = shiftdim(dz);
Nx = length(dx)-2;
Nz = length(dz)-2;

% Number the phi grid 
np = (Nx+2)*(Nz+2);
GRDp = reshape(1:np,Nx+2,Nz+2);


% Number the Ax grid
nax = (Nx+1)*(Nz+2); 
GRDax = reshape(1:nax, (Nx+1),(Nz+2));


% Number the Az grid
naz = (Nx+2)*(Nz+1); 
GRDaz = reshape(1:naz, (Nx+2),(Nz+1));


%%%%   Generate d/dx  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lx = []; jx = []; kx = [];

ex = ones(Nx+2,1);

ez = ones(Nz+2,1);

Dx = kron(dx',ez)';


lx = mkvc(GRDax);
jx = mkvc(GRDp(1:end-1,:));
kx = mkvc(-2./(Dx(1:end-1,:) + Dx(2:end,:)));

lx = [lx; lx];
jx = [jx;mkvc(GRDp(2:end,:))];
kx = [kx;-kx];



%%%%   Generate d/dz  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lz = []; jz = []; kz = [];

Dz = kron(ex',dz)';

% Entries (l,j,k)

lz = mkvc(GRDaz);
jz = mkvc(GRDp(:,1:end-1));
kz = mkvc(-2./(Dz(:,1:end-1) + Dz(:,2:end)));;


lz = [lz; lz];
jz = [jz;mkvc(GRDp(:,2:end))];
kz = [kz; -kz];

Gx = sparse(lx,jx,kx,nax,np);
 
Gz = sparse(lz,jz,kz,naz,np);


G = [Gx;Gz];

