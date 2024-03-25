%% compute the exact pressure at a given coordinate (x,y)
function z=pexact2(x,y)

z = 50*(x.*y.*(1-x-y+x.*y)-1/36);

