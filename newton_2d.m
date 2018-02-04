clear all
close all

% function
Fx_af = @(X,Y) Y - X ;
dFxdx_af = @(X,Y) - ones(size(X)) ;
dFxdy_af = @(X,Y) ones(size(X)) ;
Fy_af = @(X,Y) X.^2 + Y.^2 - 1 ;
dFydx_af = @(X,Y) 2*X ;
dFydy_af = @(X,Y) 2*Y ;

f = @(x) [Fx_af(x(1),x(2)); Fy_af(x(1),x(2))] ;
Df = @(x) [dFxdx_af(x(1),x(2)) dFxdy_af(x(1),x(2)); ...
           dFydx_af(x(1),x(2)) dFydy_af(x(1),x(2))] ;

% make 2d mesh
x = -5:0.05:6 ;
y = -5:0.05:5 ;
[X,Y] = meshgrid(y,x) ;

% Calculate functions
Fx = Fx_af(X,Y) ;
Fy = Fy_af(X,Y) ;
figure ;
surf(X,Y,Fx,'EdgeColor','none')
hold on
contour3(X,Y,Fx,-10:10,'k')
contour3(X,Y,Fx,-0.1:0.1:0.1,'r')

figure ;
surf(X,Y,Fy,'EdgeColor','none')
hold on
contour3(X,Y,Fy,-5:5,'k')
contour3(X,Y,Fy,-0.1:0.1:0.1,'r')

% norm of F
normF = sqrt(Fx.^2 + Fy.^2) ;

% Jacobian of F (four components)
dFxdx = dFxdx_af(X,Y) ;
dFydx = dFydx_af(X,Y) ;
dFxdy = dFxdy_af(X,Y) ;
dFydy = dFydy_af(X,Y) ;

% Newton Direction
%for ix = 1:numel(x) ;
%for iy = 1:numel(y) ;
%  local_Jacobian = [dFxdx(iy,ix) dFxdy(iy,ix) ; dFydx(iy,ix) dFydy(iy,ix)] ;
%  local_F = [Fx(iy,ix) ; Fy(iy,ix)] ;
%  local_Newton_dir = - local_Jacobian \ local_F ;
%  Newton_dirx(iy,ix) = local_Newton_dir(1) ;
%  Newton_diry(iy,ix) = local_Newton_dir(2) ;
%end
%end

figure ;
contourf(X,Y,normF,0:1:40)
hold on
%streamslice(X,Y,Newton_dirx,Newton_diry)

% chose X0 = (2,1) and use the jacobian there only
x0 = 3 ;
y0 = 1 ;
%Jacobian_at_X0 = [dFxdx_af(y0,x0) dFxdy_af(y0,x0) ; dFydx_af(y0,x0) dFydy_af(y0,x0)] ;
%Jacobian_at_X0 = [dFxdx_af(x0,y0) dFxdy_af(x0,y0) ; dFydx_af(x0,y0) dFydy_af(x0,y0)] ;
Jacobian_at_X0 = Df([x0;y0]) ;
Jacobian_at_X0_f = linfactor(Jacobian_at_X0) ;
plot(x0,y0,'xk')

% Chord Direction
for ix = 1:numel(x) ;
for iy = 1:numel(y) ;
  local_F = f([x(ix),y(iy)]) ;
  local_Chord_dir = - linfactor(Jacobian_at_X0_f,local_F) ;
  local_Jacobian = [dFxdx_af(x(ix),y(iy)) dFxdy_af(x(ix),y(iy)) ; dFydx_af(x(ix),y(iy)) dFydy_af(x(ix),y(iy))] ;
  local_Newton_dir = - local_Jacobian \ local_F ;
  Chord_dirx(ix,iy) = local_Chord_dir(1) ;
  Chord_diry(ix,iy) = local_Chord_dir(2) ;
  Newton_dirx(ix,iy) = local_Newton_dir(1) ;
  Newton_diry(ix,iy) = local_Newton_dir(2) ;
end
end
h = streamslice(X,Y,Chord_dirx,Chord_diry) ;
set( h, 'Color', 'r' )
h2 = streamslice(X,Y,Newton_dirx,Newton_diry) ;
set( h2, 'Color', 'k' )

x_n = [x0;y0] ;

for ii = 1:29
  x_n(:,ii+1) = x_n(:,ii) - Df(x_n(:,ii)) \ f(x_n(:,ii)) ;
end

plot(x_n(1,:),x_n(2,:),'o-g')
