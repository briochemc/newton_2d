function [f,Dxf,Dpf,Dxxf,Dppf,Dxpf] = myf(x,p)

myeps = 1e-7 ;
g = x(1)^2*x(2)^2*p^2*myeps ;
Dx1g = 2*x(1)*x(2)^2*p^2*myeps ;
Dx2g = 2*x(1)^2*x(2)*p^2*myeps ;
Dpg = 2*x(1)^2*x(2)^2*p*myeps ;
Dx1x1g = 2*x(2)^2*p^2*myeps ;
Dx2x2g = 2*x(1)^2*p^2*myeps ;
Dx1x2g = 4*x(1)*x(2)*p^2*myeps ;
Dpx1g = 4*x(1)*x(2)^2*p*myeps ;
Dpx2g = 4*x(1)^2*x(2)*p*myeps ;
Dppg = 2*x(1)^2*x(2)^2*myeps ;

f = [        x(1) + x(2)^2 + g
     exp(x(1)-1) + exp(p*x(2)) - 1 ] ;

Dxf = [        1+Dx1g            2*x(2)+Dx2g
             exp(x(1)-1)       p*exp(p*x(2))] ;

Dpf = [        Dpg
       x(2) * exp(p*x(2))] ;

Dppf = [        Dppg
       x(2)^2 * exp(p*x(2))] ;

Dxxf = zeros(2,2,2) ;
Dxxf(:,:,1) = [     Dx1x1g                   Dx1x2g
                    exp(x(1)-1)               0] ;
Dxxf(:,:,2) = [Dx1x2g         2+Dx2x2g
               0          p^2*exp(p*x(2))] ;

Dxpf = [Dpx1g           Dpx2g
        0 (1+p*x(2))*exp(p*x(2))] ;
