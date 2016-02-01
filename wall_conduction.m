% function that computes the rate of change of temperature within the tank wall
% uses a discretization of the transient conduction equation
% computes the spatial derivatives so only the time derivatives are returned

function Tdot = wall_conduction(T, q_in, q_out, constants)
% T comes in without ghost points
% q_in is heat flux from inner wall into the tank
% q_out is heat flux from ambient into wall

% this scheme solves for dT_dt in radial coordinates
% uses central differencing for the 1st and 2nd order derivatives

k = constants.k_w;
rho = constants.rho_w;
cv = constants.cv_w;
N_r = constants.N_rw;
D = constants.D_tank;
r_i = D/2;
t_w = constants.t_w;
r_o = r_i + t_w;

dr = (r_o - r_i)/(N_r - 1);
r = [r_i:dr:r_o];

alpha = k/(rho*cv); % thermal diffusivity

T = [0; T(:); 0]; % add ghost points
r = [0; r(:); 0]; 

% define ghost point values based on heat transfer
% (make the heat flux in/out of wall match k*dT_dr)
T(1) = T(3) - 2*dr/k*q_in;
T(end) = T(end-2) + 2*dr/k*q_out;

% planar conduction (note you have to add 2 to N_r to account for ghosts)
A = diag(-2/(dr^2)*ones(N_r+2,1)) ...
    + diag( ( 1/dr^2) * ones(N_r+1,1), +1 ) ...
    + diag( ( 1/dr^2) * ones(N_r+1,1), -1 );

% % make it cylindrical
A = A  ...
    + diag( (-1/(2*dr) ) ./r(1:end-1), +1 ) ...
    + diag( (1/(2*dr) ) ./r(2:end), -1 );

A = alpha*A;
    
Tdot = A*T;

% remove ghost points
Tdot = Tdot(2:end-1);