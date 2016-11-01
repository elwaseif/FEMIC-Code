function varargout = cubicSpline(coeff,smoothness); 
 % NATURAL CUBIC SPLINE ALGORITHM 3.4 
 % from Book: Numerical Analysis 6th edition by Richard L. Burden and 
 % J Douglas Faires; Published by Brooks/Cole Publishing Company
 %
 % To construct the cubic spline interpolant S for the function f,
 % defined at the numbers x(0) < x(1) < ... < x(n), satisfying
 % S''(x(0)) = S''(x(n)) = 0:
 %
 % INPUT:   A matrix with two columns and n rows (nx2) where column 1
 %          is the X data and column 2 is f(X) or Y values. [Xo,Yo]
 %                                                          [X1,Y1]
 %                                                          [:  : ]
 %                                                          [Xn,Yn]
 %
 % OUTPUT:  The new X and Y coordinates of the smoothed curve ready to plot.
 %
 % NOTE:    S(x) = A(J) + B(J)*( x - x(J) ) + C(J)*( x - x(J) )**2 +
 %          D(J) * ( x - x(J) )**3 for x(J) <= x < x(J + 1) 
 %
 % NOTE:    A smoothness factor of 32 works well. This spline technique 
 %          gives a more accurate curve than the built in Matlab Spline 
 %          function. In that it comes very close to the original data 
 %          points.
 
N = size(coeff,1)-1;
 
if N > 0 
    X = zeros(1,N+1);
    A = zeros(1,N+1);
    for I = 0:N 
        X(I+1) = coeff(I+1,1);
        A(I+1) = coeff(I+1,2);
    end
end

M = N - 1;
% STEP 1
H = zeros(1,M+1);
for I = 0:M
    H(I+1) = X(I+2) - X(I+1);
end

% STEP 2
% Use XA in place of ALPHA
XA = zeros(1,M+1);
for I = 1:M
    XA(I+1) = 3.0*(A(I+2)*H(I)-A(I+1)*(X(I+2)-X(I))+A(I)*H(I+1))/(H(I+1)*H(I));
end
    
% STEP 3
% STEPs 3, 4, 5 and part of 6 solve the tridiagonal system using 
% Crout reduction.
% use XL, XU, XZ in place of L, MU, Z resp.
XL = zeros(1,N+1);
XU = zeros(1,N+1);
XZ = zeros(1,N+1);
XL(1) = 1;
XU(1) = 0;
XZ(1) = 0;

% STEP 4
for I = 1:M
    XL(I+1) = 2*(X(I+2)-X(I))-H(I)*XU(I);
    XU(I+1) = H(I+1)/XL(I+1);
    XZ(I+1) = (XA(I+1)-H(I)*XZ(I))/XL(I+1);
end

% STEP 5
XL(N+1) = 1;
XZ(N+1) = 0;
B = zeros(1,N+1);
C = zeros(1,N+1);
D = zeros(1,N+1);
C(N+1) = XZ(N+1);

% STEP 6
for I = 0:M
    J = M-I;
    C(J+1) = XZ(J+1)-XU(J+1)*C(J+2);
    B(J+1) = (A(J+2)-A(J+1))/H(J+1) - H(J+1) * (C(J+2) + 2.0 * C(J+1)) / 3.0;
    D(J+1) = (C(J+2) - C(J+1)) / (3.0 * H(J+1));
end

% Use the above data A, B, C, and D, along with x to calculate the ySpline points.
xSpline = [];
xpoints = cell(1);
for i=1:length(X)-1
    xpoints{i} = X(i):(X(i+1)-X(i))/smoothness:X(i+1);
    xSpline = cat(2,xSpline,xpoints{i});
end

xSpline = unique(xSpline);
ySplineIdx = 1;
for i=1:length(xpoints)    
    for j=1:smoothness
        xCalc = xpoints{i}(j)-xpoints{i}(1);
        ySpline(ySplineIdx) = A(i) + B(i)*(xCalc) + C(i)*power(xCalc,2) + D(i)*power(xCalc,3);
        ySplineIdx=ySplineIdx+1;
    end
end
ySpline(end+1) = coeff(end,2); % add in the final y data point from the original data.

varargout{1} = xSpline;
varargout{2} = ySpline;