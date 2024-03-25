%% compute the gradient of the exact solution at a given coordinate (x,y)
function Z=uexactG2(x,y,k)

if k == 1 
    Z = [0.5*pi*sin(2*pi*x).*sin(2*pi*y); pi*cos(2*pi*y).*(sin(pi*x)).^2];
elseif k == 2
    Z = [-pi*cos(2*pi*x).*(sin(pi*y).^2); -0.5*pi*sin(2*pi*x).*sin(2*pi*y)];
end




% if k==1        % first component
%     Z = [0; -2*y];
% elseif k==2    % second component
%     Z = [0; 0];
% end