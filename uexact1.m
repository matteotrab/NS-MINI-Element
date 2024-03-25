%% compute the exact velocity at a given coordinate (x,y)
function z=uexact1(x,y,k)

if k==1        % first component
    z = -(cos(pi*(x-0.5)).^2).*cos(pi*(y-0.5)).*sin(pi*(y-0.5))/2;
elseif k==2    % second component
    z = (cos(pi*(y-0.5)).^2).*cos(pi*(x-0.5)).*sin(pi*(x-0.5))/2;
end

