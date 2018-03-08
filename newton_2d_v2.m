clear all
close all

ieps = 1e-10 ;
p = 1 + i * ieps ;


cost = @(chi) sum((chi - [1;1]).^2) ;

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
surf(X,Y,real(Fx)','EdgeColor','none')
hold on
contour3(X,Y,real(Fx)',-10:10,'k')
contour3(X,Y,real(Fx)',-0.1:0.1:0.1,'r')

figure ;
surf(X,Y,real(Fy)','EdgeColor','none')
hold on
contour3(X,Y,real(Fy)',-5:5,'k')
contour3(X,Y,real(Fy)',-0.1:0.1:0.1,'r')

% norm of F
normF = sqrt(real(Fx).^2 + real(Fy).^2) ;

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


chi_start = [1.1;3] ;
chi_Chord_n = chi_start ;
chi_Newton_n = chi_Chord_n ;
for ii = 1:49
  chi_Chord = chi_Chord_n(:,ii) ;
  chi_Newton = chi_Newton_n(:,ii) ;
  chi_Chord_n(:,ii+1) = chi_Chord - Df(chi0) \ f(chi_Chord) ;
  chi_Newton_n(:,ii+1) = chi_Newton - Df(chi_Newton) \ f(chi_Newton) ;
end
[sol, it_hist, ierr, chi_Kelley_n] = nsold(f,Df,chi_start) ;
plot(chi_Chord_n(1,:),chi_Chord_n(2,:),'o-r','linewidth',2)
plot(chi_Newton_n(1,:),chi_Newton_n(2,:),'o-k','linewidth',2)
plot(chi_Kelley_n(1,:),chi_Kelley_n(2,:),'o-y','linewidth',2)
xlim(xlim_save) ;
ylim(ylim_save) ;


figure
for ii = 1:50
  fchi_Chord_n(:,ii) = f(chi_Chord_n(:,ii)) ;
  fchi_Newton_n(:,ii) = f(chi_Newton_n(:,ii)) ;
  cchi_Chord_n(ii) = cost(chi_Chord_n(:,ii)) ;
  cchi_Newton_n(ii) = cost(chi_Newton_n(:,ii)) ;
end
for ii = 1:size(chi_Kelley_n,2)
  fchi_Kelley_n(:,ii) = f(chi_Kelley_n(:,ii)) ;
  cchi_Kelley_n(ii) = cost(chi_Kelley_n(:,ii)) ;
end
clear h
h(1) = semilogy(sqrt(sum(real(fchi_Chord_n).^2))) ;
hold on
h(2) = semilogy(sqrt(sum(real(fchi_Newton_n).^2))) ;
h(3) = semilogy(sqrt(sum(real(fchi_Kelley_n).^2))) ;
iend = size(fchi_Kelley_n,2) ;
h(4) = semilogy(iend,sqrt(sum(real(fchi_Kelley_n(:,iend)).^2)),'xk') ;
legend(h,{'Chord','Newton','Kelley (Default nsold)','Armijo fail'})
title('|f(x_n,p)|, real part')

figure ;
h(1) = semilogy(sqrt(sum(real(fchi_Chord_n).^2))/sqrt(sum(real(fchi_Chord_n(:,1)).^2))) ;
hold on
h(2) = semilogy(sqrt(sum(real(fchi_Newton_n).^2))/sqrt(sum(real(fchi_Newton_n(:,1)).^2))) ;
h(3) = semilogy(sqrt(sum(real(fchi_Kelley_n).^2))/sqrt(sum(real(fchi_Kelley_n(:,1)).^2))) ;
ax = gca;
ax.ColorOrderIndex = 1;
h(4) = semilogy(sqrt(sum(imag(fchi_Chord_n).^2))/sqrt(sum(imag(fchi_Kelley_n(:,1)).^2)),'--') ;
h(5) = semilogy(sqrt(sum(imag(fchi_Newton_n).^2))/sqrt(sum(imag(fchi_Kelley_n(:,1)).^2)),'--') ;
h(6) = semilogy(sqrt(sum(imag(fchi_Kelley_n).^2))/sqrt(sum(imag(fchi_Kelley_n(:,1)).^2)),'--') ;
legend(h,{'Chord R','Newton R','Kelley (Default nsold) R','Chord I','Newton I','Kelley (Default nsold) I'})
grid on
title('relative error in |f(x_n,p)|, real and imaginary part')

figure ;
h(1) = semilogy(abs((real(cchi_Chord_n)-real(cchi_Chord_n(:,end)))/real(cchi_Chord_n(:,end))),'-o') ;
hold on
h(2) = semilogy(abs((real(cchi_Newton_n)-real(cchi_Newton_n(:,end)))/real(cchi_Newton_n(:,end))),'-o') ;
h(3) = semilogy(abs((real(cchi_Kelley_n)-real(cchi_Kelley_n(:,end)))/real(cchi_Kelley_n(:,end))),'-o') ;
ax = gca;
ax.ColorOrderIndex = 1;
h(4) = semilogy(abs((imag(cchi_Chord_n)-imag(cchi_Chord_n(:,end)))/imag(cchi_Chord_n(:,end))),'-x') ;
h(5) = semilogy(abs((imag(cchi_Newton_n)-imag(cchi_Newton_n(:,end)))/imag(cchi_Newton_n(:,end))),'-x') ;
h(6) = semilogy(abs((imag(cchi_Kelley_n)-imag(cchi_Kelley_n(:,end)))/imag(cchi_Kelley_n(:,end))),'-x') ;
legend(h,{'Chord R','Newton R','Kelley (Default nsold) R','Chord I','Newton I','Kelley (Default nsold) I'})
grid on
title('relative error in cost, real and imaginary part')


figure ;
h(1) = semilogy(abs(diff(real(cchi_Chord_n))/real(cchi_Chord_n(end))),'-o') ;
hold on
h(2) = semilogy(abs(diff(real(cchi_Newton_n))/real(cchi_Newton_n(end))),'-o') ;
h(3) = semilogy(abs(diff(real(cchi_Kelley_n))/real(cchi_Kelley_n(end))),'-o') ;
ax = gca;
ax.ColorOrderIndex = 1;
h(4) = semilogy(abs(diff(imag(cchi_Chord_n))/imag(cchi_Chord_n(end))),'-x') ;
h(5) = semilogy(abs(diff(imag(cchi_Newton_n))/imag(cchi_Newton_n(end))),'-x') ;
h(6) = semilogy(abs(diff(imag(cchi_Kelley_n))/imag(cchi_Kelley_n(end))),'-x') ;
legend(h,{'Chord R','Newton R','Kelley (Default nsold) R','Chord I','Newton I','Kelley (Default nsold) I'})
grid on
title('Cauchy relative diff of cost, real and imaginary part')

figure ;
h(1) = semilogy(sqrt(sum(imag(fchi_Chord_n / ieps).^2))) ;
hold on
h(2) = semilogy(sqrt(sum(imag(fchi_Newton_n / ieps).^2))) ;
h(3) = semilogy(sqrt(sum(imag(fchi_Kelley_n / ieps).^2))) ;
iend = size(fchi_Kelley_n,2) ;
h(4) = semilogy(iend,sqrt(sum(imag(fchi_Kelley_n(:,iend) / ieps).^2)),'xk') ;
legend(h,{'Chord','Newton','Kelley (Default nsold)','Armijo fail'})
title('|f(x_n,p)|, imaginary part')


figure
h(1) = semilogy(sqrt(sum(real(chi_Chord_n - chi_Chord_n(:,end)).^2))) ;
hold on
h(2) = semilogy(sqrt(sum(real(chi_Newton_n - chi_Newton_n(:,end)).^2))) ;
h(3) = semilogy(sqrt(sum(real(chi_Kelley_n - chi_Kelley_n(:,end)).^2))) ;
legend(h,{'Chord','Newton','Kelley (Default nsold)'})
title('|x_n-x_{50}|, real part')


figure ;
h(1) = semilogy(sqrt(sum(imag(chi_Chord_n - chi_Chord_n(:,end)).^2))) ;
hold on
h(2) = semilogy(sqrt(sum(imag(chi_Newton_n - chi_Newton_n(:,end)).^2))) ;
h(3) = semilogy(sqrt(sum(imag(chi_Kelley_n - chi_Kelley_n(:,end)).^2))) ;
legend(h,{'Chord','Newton','Kelley (Default nsold)'})
title('|x_n-x_{50}|, imag part')



