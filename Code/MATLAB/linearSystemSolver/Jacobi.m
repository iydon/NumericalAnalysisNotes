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
