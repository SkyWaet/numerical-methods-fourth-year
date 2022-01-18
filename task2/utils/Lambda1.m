function L1 = Lambda1(U, N, M)
global p;
global lx;
syms x;
hx = lx/N;

    L1 = zeros(M-1, N-1);
    for j=1:M-1
        for i=1:N-1
            temp1 = subs(p,x,i*hx+hx/2)*(U(j+1,i+2)-U(j+1,i+1))/(hx^2);
            temp2 = subs(p,x,i*hx-hx/2)*(U(j+1,i+1)-U(j+1,i))/(hx^2);
            L1(j, i) = temp1-temp2;
        end
    end

end