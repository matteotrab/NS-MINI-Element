function z=uexact2(x,y,k)

if k == 1
    z = (sin(pi*x)).^2.*sin(pi*y).*cos(pi*y);
elseif k == 2
    z = -(sin(pi*y)).^2.*sin(pi*x).*cos(pi*x);
end