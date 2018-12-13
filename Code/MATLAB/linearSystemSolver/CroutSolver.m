function x = CroutSolver(A, b)
    % Crout method for solving system of linear equations.
    [L,Ll,U,Uu] = Crout(diag(A), diag(A,1), diag(A,-1));
    n = length(L);
    % Solve y
    y=L; y(1)=b(1)/y(1);
    for i = 2:n
        % Ll(i-1)y(i-1) + L(i)y(i) = b(i)
        y(i) = ( b(i)-Ll(i-1)*y(i-1) ) / L(i);
    end
    x=U; x(end)=y(end)/x(end);
    for i = n-1:-1:1
        % U(i)x(i) + Uu(i)x(i+1) = y(i)
        x(i) = ( y(i)-Uu(i)*x(i+1) ) / U(i);
    end
    disp(['Answer: ', num2str(x)]);
end



function [L,Ll,U,Uu] = Crout(d, u, l)
    % L = diag(L) + diag(Ll, -1)
    % U = diag(U) + diag(Uu,  1)
    n = length(d);

    L  = d;
    Ll = l;
    U  = ones(1, n);
    Uu = 1 : n-1;

    for i = 1:n-1
        Uu(i) = u(i) / L(i);
        L(i+1) = d(i+1) - l(i)*Uu(i);
    end
end
