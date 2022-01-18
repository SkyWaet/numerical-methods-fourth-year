function U = DirectionVariablesMethod(N, M)
global lx;
global ly;
global c1;
global c2;
global d12;
global mu_bot;
global mu_top;
global mu_left;
global mu_right;
global eps;
global pk;
global k_max
syms x;
syms y;

hx = lx/N;
hy = ly/M;

disp('Схема переменных направлений:');

%вычисление tau
    delta1 = c1*4/(hx^2)*(sin(pi*hx/(2*lx)))^2;
    Delta1 = c2*4/(hx^2)*(cos(pi*hx/(2*lx)))^2;
    delta2 = d12*4/(hy^2)*(sin(pi*hy/(2*ly)))^2;
    Delta2 = d12*4/(hy^2)*(cos(pi*hy/(2*ly)))^2;
    delta = min(delta1, delta2);
    Delta = max(Delta1, Delta2);
    tau = 2/sqrt(delta*Delta);
    
%первое приближение
    k = 0;
    U_0 = zeros(M+1,N+1);
    U_inter = U_0;
    U = U_0;
    for i=0:N
        U_0(0+1, i+1) = subs(mu_bot, x, i*hx);
        U_0(M+1, i+1) = subs(mu_top, x, i*hx);
    end
    for j=1:M-1
        U_0(j+1,0+1) = subs(mu_left, y, j*hy);
        U_0(j+1,N+1) = subs(mu_right, y, j*hy);
    end
    U_O_fix = U_0;
    U_ex_sub = u_exact_sub(N, M);
    
    flag = 5;
    [f_h, ~, norm2, ~, ro] = Norms14(flag, U_ex_sub, N, M, U_0, 0);
    
    condition = 5;
    j_tab = 0;
    tabs = zeros(1,8);
    br = false;
    while true
        if not((condition>=eps)&&(k<=k_max))
            br = true;
        end
        
        k = k + 1;
        L2 = Lambda2(U_0, N, M);
    %решение на промежуточном слое
        for i=0:N
            Uj(i+1) = subs(mu_bot, x, i*hx);
        end
        U_inter(1,:) = Uj;
    
        for j=1:M-1
            matrix = ABCGst1(N, M, j, U_0, L2, tau);
            Uj = sweepMethod(matrix)';
            U_inter(j+1,:) = Uj;
        end
    
        for i=0:N
            Uj(i+1) = subs(mu_top, x, i*hx);
        end
        U_inter(M+1,:) = Uj;
    
    %решение на целом слое
        for j=0:M
            Ui(j+1) = subs(mu_left, y, j*hy);
        end
        U(:,1) = Ui;

        L1 = Lambda1(U_inter, N, M);
        for i=1:N-1
            matrix = ABCGst2(N, M, i, U_inter, L1, tau);
            Ui = sweepMethod(matrix);
            U(:, i+1) = Ui;
        end
    
        for j=0:M
            Ui(j+1) = subs(mu_right, y, j*hy);
        end
        U(:, N+1) = Ui;
        
        %условие выхода
        L_h_uk = Lambda1(U, N, M) + Lambda2(U, N, M);
        norm_k = norm_C(L_h_uk+f_h);
        condition = norm_k/norm2;
           
        if rem(k, pk)==0
            [norm52, norm54, norm55, norm56, norm57]= Norms5(k, U, U_ex_sub, U_0, U_O_fix, N, M, f_h, condition, ro);          
        end

        if ((rem(k, pk)==1)||(pk==1))&&(k~=1)
            j_tab = j_tab+1;
            norm58 = norm_C(U-U_0)/temp(6);
%            norm58 = 0;
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
    
    %таблица
    cnames = {'k', '||F-AU^(k)||', 'rel.d.','||U^(k)-u*||', 'rel.error', '||U^(k)-U^(k-1)||', 'apost.est.', 'p_k'};
    uitable('Parent', figure('Name', 'Table1', 'Position', [400 150 700 500]), 'Position', [50 50 620 400], 'Data', tabs, 'ColumnName', cnames, 'RowName', ([]));

    Norms67(N, M, U, U_ex_sub)
end

function matrixABCGst = ABCGst1(N, M, j, U, L2, tau)
global f;
global mu_left;
global mu_right;
global lx;
global ly;
global p;
syms x;
syms y;
hx = lx/N;
hy = ly/M;

    A = zeros(1,N+1);
    B = A;
    C = A;
    G = A;
    s = A;
    t = A;
    
    B(1) = -1;
    G(1) = subs(mu_left, y, j*hy);
    for i=1:N-1
        A(i+1) = subs(p,x,i*hx-hx/2)*tau/(2*(hx^2));
        B(i+1) = ( subs(p,x,i*hx+hx/2)+subs(p,x,i*hx-hx/2) )*tau/(2*(hx^2)) + 1;
        C(i+1) = subs(p,x,i*hx+hx/2)*tau/(2*(hx^2));
        G(i+1) = -U(j+1,i+1) - tau/2*(L2(j,i) + subs(subs(f,x,i*hx),y,j*hy));
    end
    B(N+1) = -1;
    G(N+1) = subs(mu_right, y, j*hy);

    s(1) = C(1)/B(1);
    t(1) = -G(1)/B(1);
    for i = 1:N
         s(i+1) = C(i+1)/(B(i+1)-A(i+1)*s(i));
         t(i+1) = (A(i+1)*t(i) - G(i+1))/(B(i+1) - A(i+1)*s(i));
    end

    matrixABCGst = [A' B' C' G' s' t'];

end

function matrixABCGst = ABCGst2(N, M, i, U, L1, tau)
global f;
global mu_bot;
global mu_top;
global lx;
global ly;
global q;
syms x;
syms y;
hx = lx/N;
hy = ly/M;

    A = zeros(1,M+1);
    B = A;
    C = A;
    G = A;
    s = A;
    t = A;
    
    B(1) = -1;
    G(1) = subs(mu_bot, x, i*hx);
    for j=1:M-1
        A(j+1) = q*tau/(2*(hy^2));
        B(j+1) = q*tau/(hy^2) + 1;
        C(j+1) = q*tau/(2*(hy^2));
        G(j+1) = -U(j+1,i+1) - tau/2*(L1(j,i) + subs(subs(f,x,i*hx),y,j*hy));
    end
    B(M+1) = -1;
    G(M+1) = subs(mu_top, x, i*hx);

    s(1) = C(1)/B(1);
    t(1) = -G(1)/B(1);

    for j = 1:M
         s(j+1) = C(j+1)/(B(j+1)-A(j+1)*s(j));
         t(j+1) = (A(j+1)*t(j) - G(j+1))/(B(j+1) - A(j+1)*s(j));
    end

    matrixABCGst = [A' B' C' G' s' t'];

end