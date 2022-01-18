function U = ZeidelMethod(N, M)

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


%вычисление оптимального параметра
    delta = c1*4/(hx^2)*(sin(pi*hx/(2*lx)))^2 + d12*4/(hy^2)*(sin(pi*hy/(2*ly)))^2;
    Delta = c2*4/(hx^2)*(cos(pi*hx/(2*lx)))^2 + d12*4/(hy^2)*(cos(pi*hy/(2*ly)))^2;
    ksi = delta/Delta;

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

    
flag = 2;
[f_h, ~, norm2, ~, ro] = Norms14(flag, U_ex_sub, N, M, U_0, ksi);    

    U = U_0;
    condition = 5;
    j_tab = 0;
    tabs = zeros(1,8);
    br = false;
    while true
        if not((condition>=eps)&&(k<=k_max))
            br = true;
        end
        k = k + 1;
        for j=1:M-1
            for i=1:N-1
                temp1 = ( subs(p,x,i*hx-hx/2)*U(j+1,i) + subs(p,x,i*hx+hx/2)*U_0(j+1,i+2) )/(hx^2);
                temp2 = q*( U(j,i+1) + U_0(j+2,i+1) )/(hy^2);
                temp3 = ( subs(p,x,i*hx-hx/2) + subs(p,x,i*hx+hx/2) )/(hx^2);
                temp4 = 2*q/(hy^2);
                U(j+1,i+1) = ( temp1 + temp2 + subs(subs(f,x,i*hx),y,j*hy) )/( temp3 + temp4 );
            end
        end

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
        
        if rem(k, pk)==0 temp = [k, norm52, condition, norm54, norm55, norm56, norm57]; end
        U_0 = U;
        
        if br
            break;
        end

    end
    
    cnames = {'k', '||F-AU^(k)||', 'rel.d.','||U^(k)-u*||', 'rel.error', '||U^(k)-U^(k-1)||', 'apost.est.', 'p_k'};
    uitable('Parent', figure('Name', 'Table1', 'Position', [400 150 700 500]), 'Position', [50 50 620 400], 'Data', tabs, 'ColumnName', cnames, 'RowName', ([]));

Norms67(N, M, U, U_ex_sub)

end