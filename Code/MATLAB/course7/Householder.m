function M = Householder(M, is_disp)
    if nargin<2, is_disp=false; end
    n = length(M);
    for k = 1:n-2
        ap = -sign(M(k+1,k)) * norm(M(k+1:n,k), 2);
        r  = sqrt(ap^2/2 - ap*M(k+1,k)/2);
        w  = zeros(n, 1);
        w(k+1) = (M(k+1,k)-ap) / (2*r);
        w(k+2:n) = M(k+2:n,k) / (2*r);
        P = eye(n) - 2*w*w';
        M = P*M*P;
        if is_disp
            disp('alpha: '); disp(ap);
            disp('r: '); disp(r);
            disp('w: '); disp(w);
            disp('P: '); disp(P);
            disp('M: '); disp(M);
            disp('------------');
        end
    end
end
