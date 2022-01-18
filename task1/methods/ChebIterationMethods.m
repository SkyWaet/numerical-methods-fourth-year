function U = ChebIterationMethods(N, M, par)
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
global f;
global eps;
global pk;
global k_max
syms x;
syms y;

hx = lx/N;
hy = ly/M;
n = 2^par;

    disp('Метод итераций с оптимальным Чебышевским параметром:');  

%нулевое приближение
    U_0 = zeros(M+1,N+1);
    for i=0:N
        U_0(0 +1, i+1) = subs(subs(mu_bot, y, 0), x, i*hx);
        U_0(M+1, i+1) = subs(subs(mu_top, y, ly), x, i*hx);
    end
    for j=1:M-1
        U_0(j+1,0+1) = subs(subs(mu_left, x, 0), y, j*hy);
        U_0(j+1,N+1) = subs(subs(mu_right, x, lx), y, j*hy);
    end
    U_0_fix = U_0; 
    U_ex_sub = u_exact_sub(N, M);
    
    
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
    
% вычисление оптимального параметра
    delta = c1*4/(hx^2)*(sin(pi*hx/(2*lx)))^2 + d12*4/(hy^2)*(sin(pi*hy/(2*ly)))^2;
    Delta = c2*4/(hx^2)*(cos(pi*hx/(2*lx)))^2 + d12*4/(hy^2)*(cos(pi*hy/(2*ly)))^2;
    for k=1:n
        sigma = cos(pi*teta(k)/(2*n));
        tau(k) = 2/(delta+Delta+(Delta-delta)*sigma);
    end
    ksi = delta/Delta;
    
    flag = 3;
    [f_h, ~, norm2, ~, ro] = Norms14(flag, U_ex_sub, N, M, U_0, ksi);

    U = U_0;
    condition = 5;
    j_tab = 0;
    tabs = zeros(1,8);
    br = false;
    k = 0;
    while true
        if not((condition>=eps)&&(k<=k_max)||(mod(k,n)~=0))
            br = true;
        end
        k = k + 1;
        k_now = mod(k, n) + 1;
        
        for j=1:M-1
            for i=1:N-1
                temp1 = subs(p,x,i*hx+hx/2)*(U_0(j+1,i+2)-U_0(j+1,i+1))/(hx^2);
                temp2 = subs(p,x,i*hx-hx/2)*(U_0(j+1,i+1)-U_0(j+1,i))/(hx^2);
                temp3 = q*(U_0(j+2,i+1)-2*U_0(j+1,i+1)+U_0(j,i+1))/(hy^2);
                temp4 = subs(subs(f,x,i*hx),y,j*hy);
                U(j+1,i+1) = U_0(j+1,i+1) + tau(k_now)*(temp1-temp2+temp3+temp4);
            end
        end

    %условие выхода
        L_h_uk = L_approx(U, N, M);
        norm_k = norm_C(L_h_uk+f_h);
        condition = norm_k/norm2;

        if rem(k, pk)==0
            [norm52, norm54, norm55, norm56, norm57]= Norms5(k, U, U_ex_sub, U_0, U_0_fix, N, M, f_h, condition, ro);          
        end

        if ((rem(k, pk)==1)||(pk==1))&&(k~=1)
            j_tab = j_tab+1;
            norm58 = norm_C(U-U_0)/temp(6);
            tabs(j_tab,:) = [temp, norm58];       
        end
        
        if rem(k, pk)==0
            temp = [k, norm52, condition, norm54, norm55, norm56, norm57]; 
        end
        
        U_0 = U;
        
        if br
            break;
        end

    end
    
    cnames = {'k', '||F-AU^(k)||', 'rel.d.','||U^(k)-u*||', 'rel.error', '||U^(k)-U^(k-1)||', 'apost.est.', 'p_k'};
    uitable('Parent', figure('Name', 'Table1', 'Position', [400 150 700 500]), 'Position', [50 50 620 400], 'Data', tabs, 'ColumnName', cnames, 'RowName', ([]));

    Norms67(N, M, U, U_ex_sub);
end