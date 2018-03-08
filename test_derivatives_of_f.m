x = randn(2,1)
p = randn(1,1)

[f,Dxf,Dpf,Dxxf,Dppf,Dxpf] = myf(x,p) ;

h = 1e-50 ;
[fpi,Dxfpi,Dpfpi,Dxxfpi,Dppfpi,Dxpfpi] = myf(x,p+i*h) ;
xi1 = x + i*h*[1;0] ;
[fxi1,Dxfxi1,Dpfxi1,Dxxfxi1,Dppfxi1,Dxpfxi1] = myf(xi1,p) ;
xi2 = x + i*h*[0;1] ;
[fxi2,Dxfxi2,Dpfxi2,Dxxfxi2,Dppfxi2,Dxpfxi2] = myf(xi2,p) ;

CSDpf = imag(fpi)/h ;
max_error = max(abs((CSDpf(:) - Dpf(:)) ./ Dpf(:))) ;
fprintf('Max relative difference in ∂f/∂p: %g\n',max_error)

CSDxf = imag([fxi1 fxi2])/h ;
max_error = max(abs((CSDxf(:) - Dxf(:)) ./ Dxf(:))) ;
fprintf('Max relative difference in ∂f/∂x: %g\n',max_error)

CSDppf = imag(Dpfpi)/h ;
max_error = max(abs((CSDppf(:) - Dppf(:)) ./ Dppf(:))) ;
fprintf('Max relative difference in ∂²f/∂p²: %g\n',max_error)

CSDxxf = zeros(2,2,2) ;
CSDxxf(:,:,1) = imag(Dxfxi1) / h ;
CSDxxf(:,:,2) = imag(Dxfxi2) / h ;
max_error = max(abs((CSDxxf(:) - Dxxf(:)) ./ Dxxf(:))) ;
fprintf('Max relative difference in ∂²f/∂x²: %g\n',max_error)

