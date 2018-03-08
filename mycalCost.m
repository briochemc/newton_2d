function [calC,DpcalC,DppcalC,xstar] = mycalCost(p,xobs,x0)

% find xstar(p)
fun = @(x) myf(x,p) ;
xstar = nsold(fun,x0) ;

[c,Dxc,Dxxc] = mycost(xstar,xobs) ;
calC = c ;

[f,Dxf,Dpf,Dxxf,Dppf,Dxpf] = myf(xstar,p) ;

Dpxstar = - Dxf \ Dpf ;

DpcalC = Dxc * Dpxstar ;

Dppxstar = - Dxf \ (Dxxf(:,:,1) * Dpxstar * Dpxstar(1) + Dxxf(:,:,2) * Dpxstar * Dpxstar(2)...
                + 2 * Dxpf * Dpxstar...
                + Dppf) ;
DppcalC = Dpxstar.' * Dxxc * Dpxstar + Dxc * Dppxstar ;
