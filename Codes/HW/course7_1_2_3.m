%% Init
len = 99;
A = diag(2*ones(len,1)) - ...
    diag(ones(len-1,1),1) - ...
    diag(ones(len-1,1),-1);
A = A / (100^2);
b = sin( (1:len) * pi / (len+1) )';

x0 = ones(len, 1);
judge = @(x,y)norm(x-y,2)/norm(y,2);
eps = 1e-8;
max_loop = 65536; % 2^2^2^2^2 or 2^(2^(2^2)).
show_steps = false;

%% Thomas' / Crout's method
x = CroutSolver(A, b); %#ok

%% Jacobi iteration
x = Jacobi(A, b, x0, eps, judge, [], max_loop, show_steps); %#ok

%% Gauss-Seidel method
x = GaussSeidel(A, b, x0, eps, judge, [], max_loop, show_steps); %#ok

%% SOR with omega=1.8
omega = 1.8 ;
x = SOR(A, b, omega, x0, eps, judge, [], max_loop, show_steps); %#ok

%% Conjugate Gradient
x = ConjugateGradient(A, b, x0, eps, judge, [], max_loop, show_steps);

%% Done

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



function y = Jacobi(A, b, x0, eps, judge, error_gen, max_loop, show_steps)
    % Jacobi method for solving system of linear equations.
    if nargin<8, show_steps=false;          end
    if nargin<7, max_loop=1024;             end
    if nargin<6, error_gen=@(x)x;           end
    if nargin<5, judge=@(x,y)norm(y-x,inf); end
    if nargin<4, eps=1e-6;                  end
    if nargin<3, x0=zeros(size(b));         end
    if isempty(show_steps), show_steps=false;          end
    if isempty(max_loop),   max_loop=1024;             end
    if isempty(error_gen),  error_gen=@(x)x;           end
    if isempty(judge),      judge=@(x,y)norm(y-x,inf); end
    if isempty(eps),        eps=1e-6;                  end
    if isempty(x0),         x0=zeros(size(b));         end
    % Main body
    X = x0;
    L = length(A);
    n = 0;
    while judge(X,x0) >= eps || n < 2
        if show_steps, disp([num2str(n), ': ', num2str(X')]); end
        if n>=max_loop, break; end
        x0 = X;
        for i = 1:L
            X(i) = - A(i,:)*x0 + b(i);
            X(i) = X(i) / A(i,i) + x0(i);
        end
        n  = n + 1;
    end
    y = X;
    disp(['The number of iterations is ', num2str(n)]);
    disp(['Error: ', num2str(judge(X,x0))]);
end



% GaussSeidel function
% Kevin(Yinuo) Huang CID:01051134 16:13-18:20 29/03/2016
% Modified by Iydon @ 09:13-09:30 29/11/2018
% Ax = b

function x = GaussSeidel(A, b, x0, eps, judge, error_gen, max_loop, show_steps)
    % __init__
    if nargin<8, show_steps=false;          end
    if nargin<7, max_loop=1024;             end
    if nargin<6, error_gen=@(x)x;           end
    if nargin<5, judge=@(x,y)norm(y-x,inf); end
    if nargin<4, eps=1e-6;                  end
    if nargin<3, x0=zeros(size(b));         end
    if isempty(show_steps), show_steps=false;          end
    if isempty(max_loop),   max_loop=1024;             end
    if isempty(error_gen),  error_gen=@(x)x;           end
    if isempty(judge),      judge=@(x,y)norm(y-x,inf); end
    if isempty(eps),        eps=1e-6;                  end
    if isempty(x0),         x0=zeros(size(b));         end
    % A and b are input matrix, err is the eps, NumOfIter is the
    % number of iterations
    % x is the output solution
    % take the size of two matrix
    [n,m] = size(A);
    [u,~] = size(b);

   % check if A is a square matrix and weather dimensions of A and b match or not
    if n ~= m
        error('Matrix A must be a square matrix');
    elseif n ~= u
        error('The number of rows of A must be the same as that of b');
    end

    counter = 0;     % actual number of iterations used
    x  = zeros(n,1); % initialise output matrix x

   % check for conditions of GaussSeidel, see if the matrix is strictly
   % diagonally dominant or not
    for i = 1:n
        s = 0;
        for j = 1:n
            if i ~= j 
                s = s + abs(A(i,j));
            end
        end
        
        % invalid when A(1,1)(2,2)...(n,n)>sum of the related row
        % OR any diagnal entry is zero which means A is not a strictly
        % diagonally dominant matrix or a positive definite matrix
        if s < abs(A(i,j))
            fprintf('The conditions of Guass-Seidel have not met\nPlease use other methods to find solution matrix x\n');
            % return;
        end
    end

    while counter <= max_loop % loop ends when exceed max no. of iterations
        % display
        if show_steps, disp([num2str(counter), ': ', num2str(x0')]); end 
        % Gauss-Seidel Iteration
        for i = 1:n
            I = [1:i-1 i+1:n];
            x(i) = error_gen( (b(i)-A(i,I)*x(I))/A(i,i) );
        end
        % calculate error and compare with eps entered
        % break if good enough
        err = judge(x0,x);
        if err < eps
           break;
        end
        x0 = x;              % assign the new x to x0
        counter = counter+1; % no. of iterations
    end

    disp(['The number of iterations is ', num2str(counter-1)]);
    disp(['Error: ', num2str(err)]);
    if counter >= max_loop
        fprintf('Maximum number of iterations has exceeded\n');
    end
end



function y = SOR(A, b, omega, x0, eps, judge, error_gen, max_loop, show_steps)
    % SOR method for solving system of linear equations.
    if nargin<9, show_steps=false;          end
    if nargin<8, max_loop=1024;             end
    if nargin<7, error_gen=@(x)x;           end
    if nargin<6, judge=@(x,y)norm(y-x,inf); end
    if nargin<5, eps=1e-6;                  end
    if nargin<4, x0=zeros(size(b));         end
    if isempty(show_steps), show_steps=false;          end
    if isempty(max_loop),   max_loop=1024;             end
    if isempty(error_gen),  error_gen=@(x)x;           end
    if isempty(judge),      judge=@(x,y)norm(y-x,inf); end
    if isempty(eps),        eps=1e-6;                  end
    if isempty(x0),         x0=zeros(size(b));         end
    % Main body
    X = x0;
    L = length(A);
    n = 0;
    while judge(X,x0) >= eps || n < 2
        if show_steps, disp([num2str(n), ': ', num2str(X')]); end
        if n>=max_loop, break; end
        x0 = X;
        for i = 1:L
            X(i) = b(i) - A(i,1:i-1)*X(1:i-1) - A(i,i+1:L)*x0(i+1:L);
            X(i) = (1-omega)*x0(i) + omega*X(i)/A(i,i);
        end
        n  = n + 1;
    end
    y = X;
    disp(['The number of iterations is ', num2str(n)]);
    disp(['Error: ', num2str(judge(X,x0))]);
end



function y = ConjugateGradient(A, b, x0, eps, judge, error_gen, max_loop, show_steps)
    % Conjugate_gradient method for solving system of linear equations.
    % Modified from [GitHub](https://github.com/Cassie1995/Conjugate-gradient-method-).
    if nargin<8, show_steps=false;          end
    if nargin<7, max_loop=1024;             end
    if nargin<6, error_gen=@(x)x;           end
    if nargin<5, judge=@(x,y)norm(y-x,inf); end
    if nargin<4, eps=1e-6;                  end
    if nargin<3, x0=zeros(size(b));         end
    if isempty(show_steps), show_steps=false;          end
    if isempty(max_loop),   max_loop=1024;             end
    if isempty(error_gen),  error_gen=@(x)x;           end
    if isempty(judge),      judge=@(x,y)norm(y-x,inf); end
    if isempty(eps),        eps=1e-6;                  end
    if isempty(x0),         x0=zeros(size(b));         end
    % Main body
    r0     = error_gen( b - A*x0 );
    p0     = error_gen( r0 );
    alpha0 = error_gen( dot(r0,p0) / dot((A*p0),p0) );
    x1     = error_gen( x0 + alpha0*p0 );
    r1     = error_gen( b - A*x1 );
    n      = 0;
    while judge(x0,x1) >= eps
        if show_steps, disp([num2str(n), ': ', num2str(x0')]); end
        if n>=max_loop, break; end
        beta0  = error_gen( -dot(r1,A*p0) / dot(p0,A*p0) );
        p1     = error_gen( r1 + beta0*p0 );
        alpha1 = error_gen( dot(r1,p1) / dot((A*p1),p1) );
        x2     = error_gen( x1 + alpha1*p1 );
        r2     = error_gen( b - A*x2 );
        r1     = r2;
        p0     = p1;
        X      = x2;
        x0     = x1;
        x1     = x2;
        n      = n + 1;
    end
    y = X;
    disp(['The number of iterations is ', num2str(n)]);
    disp(['Error: ', num2str(judge(x0,x1))]);
end
