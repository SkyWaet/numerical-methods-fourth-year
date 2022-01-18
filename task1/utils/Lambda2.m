function L2 = Lambda2(U, N, M)
global q;
global ly;
hy = ly/M;

    L2 = zeros(M-1, N-1);
    for j=1:M-1
        for i=1:N-1
            L2(j, i) = q*(U(j+2,i+1)-2*U(j+1,i+1)+U(j,i+1))/(hy^2);
        end
    end

end