h = 1e-50 ;
p = 1+randn*1e-1 ;
xobs = randn(2,1) ;
x0 = [-1;2.5]+randn*1e-3 ;

% Calculate analaytical derivatives
[calC,DpcalC,DppcalC,astar] = mycalCost(p,xobs,x0) ;

[f,Dxf,Dpf] = myf(astar,p) ;
bstar = - Dxf \ Dpf * h ;
xstar = astar + i * bstar ;
% Calculate using complex tricks
[~,Dxci] = mycost(xstar,xobs) ;
[~,Dxfi,Dpfi] = myf(xstar,p+i*h) ;

CSDppcalC = imag(- Dxci * (Dxfi \ Dpfi)) / h 
DppcalC
max_error = max(abs((CSDppcalC(:) - DppcalC(:)) ./ DppcalC(:))) ;
fprintf('Max relative difference in ∂²C/∂p²: %g\n',max_error)
