%% gradients of basis functions on reference element of P1+bubble
function [gx,gy]=gradhphih2(i,x,y)

switch i
    case 1
        gx = -1;
        gy = -1;
    case 2
        gx = 1;
        gy = 0;
    case 3
        gx = 0;
        gy = 1;
    case 4                 % bubble
        gx = 27*(y-2*x*y-y^2);
        gy = 27*(x-x^2-2*x*y);
end
