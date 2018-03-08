function [f,Dxf,Dpf,Dxxf,Dppf,Dxpf] = myf(x,p)

f = [        x(1)^2 - x(2)
     exp(x(1)-1) + exp(p*x(2)) - 1] ;

Dxf = [     2*x(1)       -1
       exp(x(1)-1)   p*exp(p*x(2))] ;

Dpf = [        0
       x(2) * exp(p*x(2))] ;

Dxxf = zeros(2,2,2) ;
Dxxf(:,:,1) = [     2        0
               exp(x(1)-1)   0] ;
Dxxf(:,:,2) = [0         0
               0   p^2*exp(p*x(2))] ;

Dxpf = [0 0
        0 (1+p*x(2))*exp(p*x(2))] ;
