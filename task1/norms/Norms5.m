function [norm52, norm54, norm55, norm56, norm57]= Norms5(~, U, U_ex_sub, U_0, U_O_fix, N, M, f_h, ~, ro)
        L_h_k = L_approx(U, N, M);
        norm52 = norm_C(L_h_k+f_h);
        norm54 = norm_C(U-U_ex_sub);
        norm55 = norm54/norm_C(U_O_fix-U_ex_sub);
        norm56 = norm_C(U-U_0);
        norm57 = ro/(1-ro)*norm56;
end