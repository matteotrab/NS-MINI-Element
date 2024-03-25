%% compute the exact pressure at a given coordinate (x,y)
function z=pexact1(x,y)

z = (sin(pi*(x-0.5)) - sin(pi*(y-0.5)));