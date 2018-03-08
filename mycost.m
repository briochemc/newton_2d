function [c,Dxc,Dxxc] = mycost(x,xobs)

c = (x - xobs).' * (x - xobs) ;
Dxc = 2 * (x - xobs).' ;
Dxxc = 2 * eye(2) ;
