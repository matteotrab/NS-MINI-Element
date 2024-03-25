%% basis functions on reference element of P1+bubble
function z=phih2(i,xh,yh)

switch i
    case 1
        z = 1 - xh - yh;
    case 2
        z = xh;
    case 3
        z = yh;
    case 4                 % bubble
        z = 27*xh*yh*(1-xh-yh);
end
