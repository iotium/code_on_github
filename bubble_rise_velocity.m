% bubble rise velocity
function u_rise = bubble_rise_velocity(r_q, T)

paramfitvals = ...
[   0.000016760525593  -0.002762373944479;
   0.000388390480938  -0.069244519862044;
   0.002600683548867  -0.485479569513648;
   0.003858997028662  -0.718475683744660;
   0.000002383141343  -0.000105490890616;
   0.000099360365247  -0.014875924816018;
   0.001289325930961  -0.234305313015507;
   0.005379762723785  -1.056540290396839];

paramfitvals = 1e5*paramfitvals';

params = paramfitvals(1,:)*T + paramfitvals(2,:);

logr = log(r_q);

logu = params(1)*logr.^3 + params(2)*logr.^2 + params(3)*logr + params(4);
logu = logu./(logr.^4 + params(5)*logr.^3 + params(6)*logr.^2 + params(7)*logr + params(8));
% logu = zeros(size(r));
% exponent = [9:-1:0];
% for i = 1:10
%     logu = logu + params(i)*logr.^exponent(i);
% end


u_rise = exp(logu);


% if length(r_q) > 1
%     [r, sort_ind] = sort(r_q(:),'descend');
% else
%     r = r_q;
% end
% 
% g = 9.81;
% 
% u_rise1 = sqrt( 2.14*sigma./( rho_l*(2*abs(r) ) + 0.505*g*( 2*abs(r) ) ) );
% 
% u_rise2 = 0.71*sqrt( g * 2 * r);
% 
% u_rise3 = 2*rho_l*g*r.^2/(9*sigma);
% 
% Mo = g*mu^4*(rho_l - rho_v)/(rho_l^2 * sigma^3);
% 
% for i = 1:length(r)
%     if i == 1
% %     IC = min([u_rise2(i) u_rise3(i)]);
% % IC = u_rise2(i);
%         IC = u_rise2(i)*( rho_l^2 * 4 * r(i)^2 / (sigma*mu) )^(1/3);
%     else
% %         IC = u_rise4(i-1);
%         IC = V(i-1);
%     end
% %     [Cd, err_val(i)] = fminbnd(@(Cd) error_fn(Cd, Mo, r_q(i), rho_l, rho_v, mu), 1e-2, 1e3);
% %     u_rise4(i) = sqrt( 2/3 * g *(1 - rho_v/rho_l)*r_q(i)/Cd);
% 
%     [V(i), err_val(i)] = fminbnd(@(V) error_fn(V, Mo, r(i), rho_l, sigma, mu), IC*1e-1, IC*1e1);
%     u_rise4(i) = V(i)/( rho_l^2 * 4 * r(i)^2 / (sigma*mu) )^(1/3);
%  
% end
% 
% if length(r_q) > 1
%     u_rise(sort_ind) = u_rise4;
%     u_rise = reshape(u_rise, size(r_q));
%     
% else
%     u_rise = u_rise4;
% end
% 
% function E = error_fn(V, Mo, r_bub, rho_l, sigma, mu)
% g = 9.81;
% 
% % V = u*( rho_l^2 * 4 * r_bub^2 / (sigma*mu) )^(1/3);
% 
% F = g * ( rho_l^5 * (2*r_bub)^8/( sigma * mu^4 ) )^(1/3);
% 
% E = abs( (V - F/12 *( (1 + 1.31e-5*Mo^(11/20)*F^(73/33))^(21/176) )/( (1 + 0.020*F^(10/11))^(10/11)) )/V );


% function E = error_fn(Cd, Mo, r_bub, rho_l, rho_v, mu)
% 
% g = 9.81;
% 
% u_rise = sqrt( 2/3 * g *(1 - rho_v/rho_l)*r_bub/Cd);
% 
% Re = rho_l*u_rise*2*r_bub/mu;
% 
% Y = (3/4*Cd*Re^2*Mo)^(8/9);
% 
% X = (1 + 1.31e-5*Mo^(11/20)*Y^(73/33) );
% 
% E = abs(Cd - 16/Re * ( 1 + 0.02*Y^(10/11) )^(10/11)/X^(21/176))/Cd;
