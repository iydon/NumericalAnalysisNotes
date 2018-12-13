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
