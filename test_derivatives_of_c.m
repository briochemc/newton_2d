x = randn(2,1)
xobs= randn(2,1)

[c,Dxc,Dxxc] = mycost(x,xobs) ;

h = 1e-50 ;
xi1 = x + i*h*[1;0] ;
[cxi1,Dxcxi1,Dxxcxi1] = mycost(xi1,xobs) ;
xi2 = x + i*h*[0;1] ;
[cxi2,Dxcxi2,Dxxcxi2] = mycost(xi2,xobs) ;

CSDxc = imag([cxi1 cxi2])/h ;
max_error = max(abs((CSDxc(:) - Dxc(:)) ./ Dxc(:))) ;
fprintf('Max relative difference in ∂c/∂x: %g\n',max_error)

CSDxxc = zeros(1,2,2) ;
CSDxxc(:,:,1) = imag(Dxcxi1) / h ;
CSDxxc(:,:,2) = imag(Dxcxi2) / h ;
max_error = max(abs((CSDxxc(:) - Dxxc(:)) ./ Dxxc(:))) ;
fprintf('Max relative difference in ∂²c/∂x²: %g\n',max_error)

