function [f_h, norm1, norm2, norm3, ro] = Norms14(flag, U_ex_sub, N, M, U_0, ksi)
global eps;
global k_max;

ro = (1-ksi)/(1+ksi);

    disp('1. Мера аппроксимации точного решения ||F-Au*||:');
    f_h = f_approx(N, M);
    if (flag==1)||(flag==2)||(flag==3)||(flag==4)
        L_h_ex = L_approx(U_ex_sub, N, M);
        norm1 = norm_C(L_h_ex+f_h)
    end
    if (flag==5)
        L_h_ex = Lambda1(U_ex_sub, N, M) + Lambda2(U_ex_sub, N, M);
        norm1 = norm_C(L_h_ex+f_h)
    end
    
    disp('2. Мера аппроксимации нулевого приближения ||F-Au0||:');
    if (flag==1)||(flag==2)||(flag==3)||(flag==4)
        L_h_u0 = L_approx(U_0, N, M);    
        norm2 = norm_C(L_h_u0+f_h)
    end
    if (flag==5)
        L_h_u0 = Lambda1(U_0, N, M)+Lambda2(U_0, N, M);
        norm2 = norm_C(L_h_u0+f_h)
    end
    
    disp('3. Оценка количества итераций');
    if flag == 1
        norm3 = ceil(log(1/eps)/(2*ksi))
    end
    if flag == 2
        norm3 = ceil(log(1/eps)/(4*ksi))
    end
    if flag == 3
        norm3 = ceil(log(2/eps)/(2*sqrt(ksi)))
    end
    if flag == 4
        disp('Считаем отдельно.');
        norm3 = 0;
    end
    if flag == 5
        norm3 = ceil(N/(2*pi)*log(1/eps))
    end
    
    disp('4. Спектральный радиус матрицы перехода:');
    if (flag == 1)||(flag == 3)||(flag == 4)
        ro 
    end
    if flag == 2
        ro = ro^2
    end
    if flag == 5
        disp('-');
        ro = 0;
    end
   k_max = min(50,norm3);
end

function f_h = f_approx(N, M)
global lx;
global ly;
global f;
syms x;
syms y;
hx = lx/N;
hy = ly/M;

    f_h = zeros(N-1, M-1);
    for j=1:M-1
        for i=1:N-1
            f_h(j, i) = subs(subs(f, x, i*hx), y, j*hy);
        end
    end
end