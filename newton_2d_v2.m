clear all
close all

h = 1e-50 ;
p = 1+randn*1e-1 ;
nit = 300 ;

xobs = randn(2,1) ;

% make 2d mesh
x1 = (-5:0.05:6)+0.001 ;
x2 = -5:0.05:5 ;
[X1,X2] = meshgrid(x1,x2) ;

% Calculate functions at each point fo the grid
for ix1 = 1:numel(x1)
for ix2 = 1:numel(x2)
  x = [x1(ix1);x2(ix2)] ;
  f = myf(x,p) ;
  Fx1(ix1,ix2) = f(1) ;
  Fx2(ix1,ix2) = f(2) ;
end
end

figure ;
surf(X1,X2,real(Fx1)','EdgeColor','none')
hold on
contour3(X1,X2,real(Fx1)',-10:10,'k')
contour3(X1,X2,real(Fx1)',-0.1:0.1:0.1,'r')

figure ;
surf(X1,X2,real(Fx2)','EdgeColor','none')
hold on
contour3(X1,X2,real(Fx2)',-5:5,'k')
contour3(X1,X2,real(Fx2)',-0.1:0.1:0.1,'r')

% norm of F
normF = sqrt(real(Fx1).^2 + real(Fx2).^2) ;

figure ;
contourf(X1,X2,normF',0:1:40)
hold on
%streamslice(X1,X2,Newton_dirx1,Newton_dirx2)

% chose X10 = (2,1) and use the jacobian there onlx2
x0 = [-1;2.5]+randn*1e-3 ;
[f0,Dxf0] = myf(x0,p) ;
%Dfx0 = [dFx1dx1_af(x20,x10) dFx1dx2_af(x20,x10) ; dFx2dx1_af(x20,x10) dFx2dx2_af(x20,x10)] ;
%Dfx0 = [dFx1dx1_af(x10,x20) dFx1dx2_af(x10,x20) ; dFx2dx1_af(x10,x20) dFx2dx2_af(x10,x20)] ;
plot(x0(1),x0(2),'or','markersize',10,'linewidth',2)

% Chord Direction
for ix1 = 1:numel(x1) ;
for ix2 = 1:numel(x2) ;
  x = [x1(ix1);x2(ix2)] ;
  [f,Dxf] = myf(x,p) ;
  local_Chord_dir = - Dxf0 \ f ;
  local_Newton_dir = - Dxf \ f ;
  Chord_dirx1(ix2,ix1) = local_Chord_dir(1) ;
  Chord_dirx2(ix2,ix1) = local_Chord_dir(2) ;
  Newton_dirx1(ix2,ix1) = local_Newton_dir(1) ;
  Newton_dirx2(ix2,ix1) = local_Newton_dir(2) ;
end
end
h1 = streamslice(X1,X2,Chord_dirx1,Chord_dirx2) ;
set( h1, 'Color', 'r' )
h2 = streamslice(X1,X2,Newton_dirx1,Newton_dirx2) ;
set( h2, 'Color', 'k' )
x1lim_save = xlim ;
x2lim_save = xlim ;


x_start = x0 ;
x_Chord_n = x_start ;
x_Newton_n = x_Chord_n ;
for ii = 1:nit-1
  x_Chord = x_Chord_n(:,ii) ;
  [f_Chord] = myf(x_Chord,p) ;
  x_Chord_n(:,ii+1) = x_Chord - Dxf0 \ f_Chord ;
  x_Newton = x_Newton_n(:,ii) ;
  [f_Newton,Dxf_Newton] = myf(x_Newton,p) ;
  x_Newton_n(:,ii+1) = x_Newton - Dxf_Newton \ f_Newton ;
end
fun = @(x) myf(x,p) ;
[sol, it_hist, ierr, x_Kelley_n] = nsold(fun,x_start) ;
plot(x_Newton_n(1,:),x_Newton_n(2,:),'o-k','linewidth',2)
plot(x_Kelley_n(1,:),x_Kelley_n(2,:),'o-m','linewidth',2)
plot(x_Chord_n(1,:),x_Chord_n(2,:),'o-g','linewidth',2)
contour3(X1,X2,real(Fx1)',-0.1:0.1:0.1,'r')
contour3(X1,X2,real(Fx2)',-0.1:0.1:0.1,'r')
xlim(x1lim_save) ;
xlim(x2lim_save) ;


figure
for ii = 1:nit
  fx_Chord_n(:,ii) = myf(x_Chord_n(:,ii),p) ;
  fx_Newton_n(:,ii) = myf(x_Newton_n(:,ii),p) ;
  cx_Chord_n(ii) = mycost(x_Chord_n(:,ii)) ;
  cx_Newton_n(ii) = mycost(x_Newton_n(:,ii)) ;
end
for ii = 1:size(x_Kelley_n,2)
  fx_Kelley_n(:,ii) = myf(x_Kelley_n(:,ii),p) ;
  cx_Kelley_n(ii) = mycost(x_Kelley_n(:,ii)) ;
end
hax(1) = semilogy(sqrt(sum(real(fx_Chord_n).^2))) ;
hold on
hax(2) = semilogy(sqrt(sum(real(fx_Newton_n).^2))) ;
hax(3) = semilogy(sqrt(sum(real(fx_Kelley_n).^2))) ;
iend = size(fx_Kelley_n,2) ;
hax(4) = semilogy(iend,sqrt(sum(real(fx_Kelley_n(:,iend)).^2)),'xk') ;
legend(hax,{'Chord','Newton','Kelley (Default nsold)','Armijo fail'})
title('|f(x1_n,p)|, real part')

figure ;
hax(1) = semilogy(sqrt(sum(real(fx_Chord_n).^2))/sqrt(sum(real(fx_Chord_n(:,1)).^2))) ;
hold on
hax(2) = semilogy(sqrt(sum(real(fx_Newton_n).^2))/sqrt(sum(real(fx_Newton_n(:,1)).^2))) ;
hax(3) = semilogy(sqrt(sum(real(fx_Kelley_n).^2))/sqrt(sum(real(fx_Kelley_n(:,1)).^2))) ;
ax1 = gca;
ax1.ColorOrderIndex = 1;
hax(4) = semilogy(sqrt(sum(imag(fx_Chord_n).^2))/sqrt(sum(imag(fx_Kelley_n(:,1)).^2)),'--') ;
hax(5) = semilogy(sqrt(sum(imag(fx_Newton_n).^2))/sqrt(sum(imag(fx_Kelley_n(:,1)).^2)),'--') ;
hax(6) = semilogy(sqrt(sum(imag(fx_Kelley_n).^2))/sqrt(sum(imag(fx_Kelley_n(:,1)).^2)),'--') ;
legend(hax,{'Chord R','Newton R','Kelley (Default nsold) R','Chord I','Newton I','Kelley (Default nsold) I'})
grid on
title('relative error in |f(x1_n,p)|, real and imaginarx2 part')

figure ;
hax(1) = semilogy(abs((real(cx_Chord_n)-real(cx_Chord_n(:,end)))/real(cx_Chord_n(:,end))),'-o') ;
hold on
hax(2) = semilogy(abs((real(cx_Newton_n)-real(cx_Newton_n(:,end)))/real(cx_Newton_n(:,end))),'-o') ;
hax(3) = semilogy(abs((real(cx_Kelley_n)-real(cx_Kelley_n(:,end)))/real(cx_Kelley_n(:,end))),'-o') ;
ax1 = gca;
ax1.ColorOrderIndex = 1;
hax(4) = semilogy(abs((imag(cx_Chord_n)-imag(cx_Chord_n(:,end)))/imag(cx_Chord_n(:,end))),'-x') ;
hax(5) = semilogy(abs((imag(cx_Newton_n)-imag(cx_Newton_n(:,end)))/imag(cx_Newton_n(:,end))),'-x') ;
hax(6) = semilogy(abs((imag(cx_Kelley_n)-imag(cx_Kelley_n(:,end)))/imag(cx_Kelley_n(:,end))),'-x') ;
legend(hax,{'Chord R','Newton R','Kelley (Default nsold) R','Chord I','Newton I','Kelley (Default nsold) I'})
grid on
title('relative error in cost, real and imaginarx2 part')


figure ;
hax(1) = semilogy(abs(diff(real(cx_Chord_n))/real(cx_Chord_n(end))),'-o') ;
hold on
hax(2) = semilogy(abs(diff(real(cx_Newton_n))/real(cx_Newton_n(end))),'-o') ;
hax(3) = semilogy(abs(diff(real(cx_Kelley_n))/real(cx_Kelley_n(end))),'-o') ;
ax1 = gca;
ax1.ColorOrderIndex = 1;
hax(4) = semilogy(abs(diff(imag(cx_Chord_n))/imag(cx_Chord_n(end))),'-x') ;
hax(5) = semilogy(abs(diff(imag(cx_Newton_n))/imag(cx_Newton_n(end))),'-x') ;
hax(6) = semilogy(abs(diff(imag(cx_Kelley_n))/imag(cx_Kelley_n(end))),'-x') ;
legend(hax,{'Chord R','Newton R','Kelley (Default nsold) R','Chord I','Newton I','Kelley (Default nsold) I'})
grid on
title('Cauchx2 relative diff of cost, real and imaginarx2 part')

figure ;
hax(1) = semilogy(sqrt(sum(imag(fx_Chord_n / h).^2))) ;
hold on
hax(2) = semilogy(sqrt(sum(imag(fx_Newton_n / h).^2))) ;
hax(3) = semilogy(sqrt(sum(imag(fx_Kelley_n / h).^2))) ;
iend = size(fx_Kelley_n,2) ;
hax(4) = semilogy(iend,sqrt(sum(imag(fx_Kelley_n(:,iend) / h).^2)),'xk') ;
legend(hax,{'Chord','Newton','Kelley (Default nsold)','Armijo fail'})
title('|f(x1_n,p)|, imaginarx2 part')


figure ;
hax(1) = semilogy(sqrt(sum(real(x_Chord_n - x_Chord_n(:,end)).^2))) ;
hold on
hax(2) = semilogy(sqrt(sum(real(x_Newton_n - x_Newton_n(:,end)).^2))) ;
hax(3) = semilogy(sqrt(sum(real(x_Kelley_n - x_Kelley_n(:,end)).^2))) ;
legend(hax,{'Chord','Newton','Kelley (Default nsold)'})
title('|x1_n-x1_{50}|, real part')


figure ;
hax(1) = semilogy(sqrt(sum(imag(x_Chord_n - x_Chord_n(:,end)).^2))) ;
hold on
hax(2) = semilogy(sqrt(sum(imag(x_Newton_n - x_Newton_n(:,end)).^2))) ;
hax(3) = semilogy(sqrt(sum(imag(x_Kelley_n - x_Kelley_n(:,end)).^2))) ;
legend(hax,{'Chord','Newton','Kelley (Default nsold)'})
title('|x1_n-x1_{50}|, imag part')

astar = x_Kelley_n(:,end) ;
[fstar,Dxfstar,Dpfstar,Dxxfstar,Dppfstar,Dxpfstar] = myf(astar,p)

bstar = - Dxfstar \ (Dpfstar * h) ;

xstar = astar + i * bstar ;



[~,Dxcstari] = mycost(xstar,p+i*h) ;

[~,Dxfstari] = myf(xstar,p+i*h) ;
CSDxxfstar = imag(Dxfstar) / h
Dxxfstar
max_error = max(abs((CSDxxf(:) - Dxxf(:)) ./ Dxxf(:))) ;
fprintf('Max relative difference in ∂²f/∂x²: %g\n',max_error)



