function Norms67(N, M, U, U_ex_sub)

disp('6. Приближенное решение на крупной сетке:');
    stepx = N/5;
    stept = M/5;
    U(1:stept:M+1, 1:stepx:N+1)
    
disp('7. Таблица точного решения на крупной сетке:');
    U_ex_sub(1:stept:M+1, 1:stepx:N+1)
    
end