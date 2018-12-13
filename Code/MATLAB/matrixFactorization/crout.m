function [L,U] = crout(a, c, d)
    n  = length(a);

    L = zeros(n);
    U = zeros(n);
    p = 1:n;
    q = 1:n-1;

    p(1) = a(1);
    for i = 1:n-1
        q(i) = c(i)/p(i);
        p(i+1) = a(i+1) - d(i)*q(i);
    end

    for i = 1:n
        L(i,i) = p(i);
        U(i,i) = 1;
    end

    for i = 1:n-1
        L(i+1,i) = d(i);
        U(i,i+1) = q(i);
    end
end
