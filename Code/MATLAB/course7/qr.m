function [Q,R] = qr(A, is_disp)
    % https://github.com/blackxparade/QR-algorithm.
    if nargin<2, is_disp=false; end

    [m,n] = size(A);
    Q     = eye(n);

    for i = 1:n-1
        B = eye(n);
        B(i:end, i:end) = HH( norm(A(i:end, i), 2)*eye(n+1-i,1) - A(i:end, i) );
        Q = Q * B;
        A = B * A;
        if is_disp
            disp('Q: '); disp(B);
            disp('R: '); disp(A);
            disp('------------');
        end
    end
    R = A;
end



function A = HH(v)
    % Householder matrix formula:
    %       2
    % I - ----- h*h'
    %      h'*h
    v = v(:);
    I = eye(length(v));
    A = I - ( (2/(v' * v)) * (v * v') );
end
