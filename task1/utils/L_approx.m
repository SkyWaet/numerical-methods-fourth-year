function L_h = L_approx(U, N, M)
global p;
global q;
global lx;
global ly;
syms x;
hx = lx/N;
hy = ly/M;

L_h = zeros(M-1, N-1);
    for j=1:M-1
        for i=1:N-1
            temp1 = subs(p,x,i*hx+hx/2)*(U(j+1,i+2)-U(j+1,i+1))/(hx^2);
            temp2 = subs(p,x,i*hx-hx/2)*(U(j+1,i+1)-U(j+1,i))/(hx^2);
            temp3 = q*(U(j+2,i+1)-2*U(j+1,i+1)+U(j,i+1))/(hy^2);
            L_h(j, i) = temp1-temp2+temp3;
        end
    end
end