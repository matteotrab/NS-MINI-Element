%% compute the loading term at a given coordinate (x,y)
function z=carico1(x,y,k,nu)

%nu=1;

if k==1        % first component
 z = nu*(pi^2)*sin(pi*(y-0.5))*cos(pi*(y-0.5))*(1-4*(cos(pi*(x-0.5)))^2) ...
     - pi*cos(pi*(x-0.5));
elseif k==2    % second component
 z = nu*(-pi^2)*sin(pi*(x-0.5))*cos(pi*(x-0.5))*(1-4*(cos(pi*(y-0.5)))^2) ...
     + pi*cos(pi*(y-0.5));
end

% ALTERNATIVE VERSION
% u1 = -(cos(pi*(x-0.5))^2)*cos(pi*(y-0.5))*sin(pi*(y-0.5))/2;
% u2 = (cos(pi*(y-0.5))^2)*cos(pi*(x-0.5))*sin(pi*(x-0.5))/2;
% 
% x = pi*(x-0.5);
% y = pi*(y-0.5);
% 
% u1x=pi*cos(x)*sin(x)*cos(y)*sin(y);
% u1y=pi*(-1/2*(cos(x)^2)*(-sin(y)^2 + cos(y)^2));
% u2x=pi*1/2*(cos(y)^2)*(-sin(x)^2 + cos(x)^2);
% u2y=pi*(-cos(x)*sin(x)*cos(y)*sin(y));
% 
% 
% if k==1        % prima componente
%  z = z + u1*u1x + u2*u1y;
% elseif k==2    % seconda componente
%  z = z + u1*u2x + u2*u2y;
% end





