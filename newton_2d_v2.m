clear all
close all

f = @(chi) [chi(1) + chi(2) ; chi(1)^2 + chi(2)^2 - 1] ;
Df = @(chi) [1 1 ; 2*chi(1) 2*chi(2)] ;

% make 2d mesh
x = (-5:0.05:6)+0.001 ;
y = -5:0.05:5 ;
[X,Y] = meshgrid(x,y) ;

% Calculate functions at each point fo the grid
for ix = 1:numel(x)
for iy = 1:numel(y)
  chi = [x(ix);y(iy)] ;
  fchi = f(chi) ;
  Fx(ix,iy) = fchi(1) ;
  Fy(ix,iy) = fchi(2) ;
  Dfchi = Df(chi) ;
  DFxx(ix,iy) = Dfchi(1,1) ;
  DFxy(ix,iy) = Dfchi(1,2) ;
  DFyx(ix,iy) = Dfchi(2,1) ;
  DFyy(ix,iy) = Dfchi(2,2) ;
end
end

figure ;
surf(X,Y,Fx','EdgeColor','none')
hold on
contour3(X,Y,Fx',-10:10,'k')
contour3(X,Y,Fx',-0.1:0.1:0.1,'r')

figure ;
surf(X,Y,Fy','EdgeColor','none')
hold on
contour3(X,Y,Fy',-5:5,'k')
contour3(X,Y,Fy',-0.1:0.1:0.1,'r')

% norm of F
normF = sqrt(Fx.^2 + Fy.^2) ;

figure ;
contourf(X,Y,normF',0:1:40)
hold on
%streamslice(X,Y,Newton_dirx,Newton_diry)

% chose X0 = (2,1) and use the jacobian there only
chi0 = [3+0.001;1] ;
%Dfchi0 = [dFxdx_af(y0,x0) dFxdy_af(y0,x0) ; dFydx_af(y0,x0) dFydy_af(y0,x0)] ;
%Dfchi0 = [dFxdx_af(x0,y0) dFxdy_af(x0,y0) ; dFydx_af(x0,y0) dFydy_af(x0,y0)] ;
plot(chi0(1),chi0(2),'or','markersize',10,'linewidth',2)

% Chord Direction
for ix = 1:numel(x) ;
for iy = 1:numel(y) ;
  chi = [x(ix);y(iy)] ;
  local_Chord_dir = - Df(chi0) \ f(chi) ;
  local_Newton_dir = - Df(chi) \ f(chi) ;
  Chord_dirx(iy,ix) = local_Chord_dir(1) ;
  Chord_diry(iy,ix) = local_Chord_dir(2) ;
  Newton_dirx(iy,ix) = local_Newton_dir(1) ;
  Newton_diry(iy,ix) = local_Newton_dir(2) ;
end
end
h = streamslice(X,Y,Chord_dirx,Chord_diry) ;
set( h, 'Color', 'r' )
h2 = streamslice(X,Y,Newton_dirx,Newton_diry) ;
set( h2, 'Color', 'k' )
xlim_save = xlim ;
ylim_save = ylim ;

chi_start = [1.1;3+i*1e-10] ;
chi_Chord_n = chi_start ;
chi_Newton_n = chi_Chord_n ;
for ii = 1:49
  chi_Chord = chi_Chord_n(:,ii) ;
  chi_Newton = chi_Newton_n(:,ii) ;
  chi_Chord_n(:,ii+1) = chi_Chord - Df(chi0) \ f(chi_Chord) ;
  chi_Newton_n(:,ii+1) = chi_Newton - Df(chi_Newton) \ f(chi_Newton) ;
end
[sol, it_hist, ierr, x_hist] = nsold(f,Df,chi_start) ;
plot(chi_Chord_n(1,:),chi_Chord_n(2,:),'o-r','linewidth',2)
plot(chi_Newton_n(1,:),chi_Newton_n(2,:),'o-k','linewidth',2)
plot(x_hist(1,:),x_hist(2,:),'o-y','linewidth',2)
xlim(xlim_save) ;
ylim(ylim_save) ;


figure
for ii = 1:50
  fchi_Chord_n(:,ii) = f(chi_Chord_n(:,ii)) ;
  fchi_Newton_n(:,ii) = f(chi_Newton_n(:,ii)) ;
end
for ii = 1:size(x_hist,2) 
  fchi_Kelley_n(:,ii) = f(x_hist(:,ii)) ;
end
clear h
h(1) = semilogy(sqrt(sum(real(fchi_Chord_n).^2))) ;
hold on
h(2) = semilogy(sqrt(sum(real(fchi_Newton_n).^2))) ;
h(3) = semilogy(sqrt(sum(real(fchi_Kelley_n).^2))) ;
iend = size(fchi_Kelley_n,2) ;
h(4) = semilogy(iend,sqrt(sum(real(fchi_Kelley_n(:,iend)).^2)),'xk') ;
legend(h,{'Chord','Newton','Kelley (Default nsold)','Armijo fail'})

figure ;
h(1) = semilogy(sqrt(sum(imag(fchi_Chord_n).^2))) ;
hold on
h(2) = semilogy(sqrt(sum(imag(fchi_Newton_n).^2))) ;
h(3) = semilogy(sqrt(sum(imag(fchi_Kelley_n).^2))) ;
iend = size(fchi_Kelley_n,2) ;
h(4) = semilogy(iend,sqrt(sum(imag(fchi_Kelley_n(:,iend)).^2)),'xk') ;
legend(h,{'Chord','Newton','Kelley (Default nsold)','Armijo fail'})


