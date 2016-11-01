function [k_corner,cornerval,nx,ny] = FEMIC_gcorner(xx,yy)
%keyboard;
x = log10(xx);
y = log10(yy);
minX = min(x);%(min(x)-(max(x)-min(x))/5);
maxX = max(x);%(max(x)+(max(x)-min(x))/5);
dx = (max(x)-min(x))/255;
xp = [minX:dx:maxX]';
yp = spline(x,y,xp);
grr = -2.*ones(length(xp),1);
for i=21:length(xp),
    grr(i) = sum( diff( yp(i-20:i),2 ) );
end
% while icc<=21,
%     [mxc,icc] = max(grr(ix:end));
%     ix=icc;
%     return
% end
%keyboard
[mxc,icc] = max(grr);
k_corner = icc;
cornerval = 10.^xp(k_corner);
nx = 10.^xp;
ny = 10.^yp;
