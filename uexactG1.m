%% compute the gradient of the exact solution at a given coordinate (x,y)
function Z=uexactG1(x,y,k)

x = pi*(x-0.5);
y = pi*(y-0.5);
if k==1        % first component
    Z = pi*[cos(x)*sin(x)*cos(y)*sin(y); -1/2*(cos(x)^2)*(-sin(y)^2 + cos(y)^2)];
elseif k==2    % second component
    Z = pi*[1/2*(cos(y)^2)*(-sin(x)^2 + cos(x)^2);-cos(x)*sin(x)*cos(y)*sin(y)];
end