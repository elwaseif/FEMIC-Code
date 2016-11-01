function K = kron3(A,B,C)
% K = kron3(A,B,C)
%calculates a 3D kronecker product
%A, B, C are vectors of arbitrary length

% Copyright (c) 2007 by the Society of Exploration Geophysicists.
% For more information, go to http://software.seg.org/2007/0001 .
% You must read and accept usage terms at:
% http://software.seg.org/disclaimer.txt before use.
% 
% Revision history:
% Original SEG version Eldad Haber
% Last update, July 2003

ma = length(A);
mb = length(B);
mc = length(C);

A = reshape(A,ma,1,1);
B = reshape(B,1,mb,1);
C = reshape(C,1,1,mc);

[ma,na,ka] = size(A);
[mb,nb,kb] = size(B);
[mc,nc,kc] = size(C);


t = 0:(ma*mb*mc-1);
ia = fix(t/(mb*mc))+1;
ib = rem(t,mb)+1;
ic = rem(t,mc)+1;

t = 0:(na*nb*nc-1);
ja = fix(t/(nb*nc))+1;
jb = rem(t,nb)+1;
jc = rem(t,nc)+1;

t = 0:(ka*kb*kc-1);
ka = fix(t/(kb*kc))+1;
kb = rem(t,kb)+1;
kc = rem(t,kc)+1;

K = A(ia,ja,ka).*B(ib,jb,kb).*C(ic,jc,kc);


