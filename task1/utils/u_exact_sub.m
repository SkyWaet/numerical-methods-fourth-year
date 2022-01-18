function u_ex = u_exact_sub(N, M)
global lx;
global ly;
global U_exact;
hx = lx/N;
hy = ly/M;
syms x;
syms y;

    u_ex = zeros(N+1, M+1);
    for j=0:M
        for i=0:N
            u_ex(j+1,i+1) = subs(subs(U_exact, x, i*hx), y, j*hy);
        end
    end
     
end