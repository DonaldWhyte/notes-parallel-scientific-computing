clear all;
close all;

load heat_solution.dat

  sol = heat_solution(:,3);

  dimensions = size(sol);
  nnodes     = dimensions(1);
  nx         = sqrt(nnodes);
  ny         = nx;

  SOL = reshape(sol,ny,nx);

  x = -2:(4.0/(nx-1)):2;
  y = -2:(4.0/(ny-1)):2;
  [X,Y] = meshgrid(x,y);

mesh(X,Y,SOL)
set(gca,'PlotBoxAspectRatio',[1,1,1])
xlabel('x')
ylabel('y')
