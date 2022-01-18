function U = AlternateTriangularMethod(N, M, par)
global lx;
global ly;
global c1;
global c2;
global d12;
global p;
global q;
global mu_bot;
global mu_top;
global mu_left;
global mu_right;
global eps;
global pk;
global k_max
syms x;
syms y;

n = 2^par;
hx = lx/N;
hy = ly/M;

disp('Поперемено-треугольный метод с оптимальными чебышевскими параметрами:');

%вычисление упорядоченных параметров
    teta_prev = [1];
    m=1;
    while 2*m<=n
        for i=1:m
            teta(2*i-1) = teta_prev(i);
            teta(2*i) = 4*m - teta(2*i-1);
        end
        m = m*2;
        teta_prev = teta;
    end

%вычисление оптимального параметра
    delta = c1*4/(hx^2)*(sin(pi*hx/(2*lx)))^2 + d12*4/(hy^2)*(sin(pi*hy/(2*ly)))^2;
    Delta = c2*4/(hx^2) + d12*4/(hy^2);
    eta = delta/Delta;
    omega = 2/sqrt(delta*Delta);
    gamma1 = delta/(2+2*sqrt(eta));
    gamma2 = delta/(4*sqrt(eta));
    ksi = gamma1/gamma2;
    for k=1:n
        sigma = cos(pi*teta(k)/(2*n));
        tau(k) = 2/(gamma1+gamma2+(gamma2-gamma1)*sigma);
    end
    k1 = omega/(hx^2);
    k2 = omega/(hy^2);
    
%первое приближение
    k = 0;
    U_0 = zeros(M+1,N+1);
    for i=0:N
        U_0(0 +1, i+1) = subs(subs(mu_bot, y, 0), x, i*hx);
        U_0(M+1, i+1) = subs(subs(mu_top, y, ly), x, i*hx);
    end
    for j=1:M-1
        U_0(j+1,0+1) = subs(subs(mu_left, x, 0), y, j*hy);
        U_0(j+1,N+1) = subs(subs(mu_right, x, lx), y, j*hy);
    end
    U_O_fix = U_0;
    U_ex_sub = u_exact_sub(N, M);
    
    flag = 4;
    [f_h, ~, norm2, ~, ro] = Norms14(flag, U_ex_sub, N, M, U_0, ksi); 
    disp('3. Оценка количества итераций');
    norm3 = ceil(log(2/eps)/(2*sqrt(2)*eta^(1/4)));
    
    U = U_0;
    condition = 5;
    j_tab = 0;
    br = false;
    tabs = zeros(1,8);
    while true
        if not((condition>=eps)&&(k<=k_max))
            br = true;
        end
        k = k + 1;
        k_now = mod(k, n) + 1;
        
        LU_0 = L_approx(U_0, N, M);
        Fi = LU_0 + f_h;
        temp1 = zeros(N-1,1);
        Fi = [temp1 Fi temp1];
        temp1 = zeros(1,N+1);
        Fi = [temp1; Fi; temp1];
        
    %решаем первую систему
        w = zeros(M+1,N+1);
        for j=1:M-1
            for i=1:N-1
                temp1 = k1*subs(p,x,i*hx-hx/2)*w(j+1,i);
                temp2 = k2*q*w(j,i+1);
                temp3 = 1 + k1*subs(p,x,i*hx-hx/2) + k2*q;
                w(j+1,i+1) = (temp1 + temp2 + Fi(j+1,i+1))/temp3;
            end
        end
        
        %решаем вторую систему
        wk = zeros(M+1,N+1);
        for j=(M-1):(-1):1
            for i=(N-1):(-1):1
                temp1 = k1*subs(p,x,i*hx+hx/2)*wk(j+1,i+2);
                temp2 = k2*q*wk(j+2,i+1);
                temp3 = 1 + k1*subs(p,x,i*hx+hx/2) + k2*q;
                wk(j+1,i+1) = (temp1 + temp2 + w(j+1,i+1))/temp3;
            end
        end
        
        U = U_0 + tau(k_now)*wk;
        
        %условие выхода
            L_h_uk = L_approx(U, N, M);
            norm_k = norm_C(L_h_uk+f_h);
            condition = norm_k/norm2;
            
        if rem(k, pk)==0
            [norm52, norm54, norm55, norm56, norm57]= Norms5(k, U, U_ex_sub, U_0, U_O_fix, N, M, f_h, condition, ro);          
        end
        
        if ((rem(k, pk)==1)||(pk==1))&&(k~=1)
            j_tab = j_tab+1;
            norm58 = norm_C(U-U_0)/temp(6);
            tabs(j_tab,:) = [temp, norm58]; 
        end
        
        if rem(k, pk)==0
            temp = [k, norm52, condition, norm54, norm55, norm56, norm57];
        end
        if br
            break;
        end
        U_0 = U;
        
    end
    
    cnames = {'k', '||F-AU^(k)||', 'rel.d.','||U^(k)-u*||', 'rel.error', '||U^(k)-U^(k-1)||', 'apost.est.', 'p_k'};
    uitable('Parent', figure('Name', 'Table1', 'Position', [400 150 700 500]), 'Position', [50 50 620 400], 'Data', tabs, 'ColumnName', cnames, 'RowName', ([]));

    Norms67(N, M, U, U_ex_sub);
end