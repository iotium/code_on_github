function m = tank_wall_mass(V,D,rho,t)

L = V/(pi*D^2/4);
A = pi*D*L + 0.25*pi*D^2;
m = A*t*rho;