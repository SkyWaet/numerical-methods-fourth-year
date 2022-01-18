function [f, mu_bot, mu_top, mu_left, mu_right] = initConditions(u)
global lx;
global ly;
syms x;
syms y;
    
    f = -( diff((3*x+2)*diff(u,x),x) + diff(diff(u,y),y) );
    mu_bot = subs(u, y, 0);
    mu_top = subs(u, y, ly);
    mu_left = subs(u, x, 0);
    mu_right = subs(u, x, lx);
    
end