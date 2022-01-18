function Y = sweepMethod(matrixABCGst)

    n = size(matrixABCGst,1) - 1;
    Y = zeros(n+1,1);

    Y(n+1) = matrixABCGst(n+1,6);

    for i=n:-1:1
         Y(i) = matrixABCGst(i,5)*Y(i+1) + matrixABCGst(i, 6);
    end

end