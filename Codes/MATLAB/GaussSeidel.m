% GaussSeidel function
% Kevin(Yinuo) Huang CID:01051134 16:13-18:20 29/03/2016
% Modified by Iydon @ 09:13-09:30 29/11/2018
% Ax = b

function x = GaussSeidel(A, b, x0, eps, judge, error_gen, max_loop)
    % __init__
    if nargin<7, max_loop=1024;         end
    if nargin<6, error_gen=@(x)x;       end
    if nargin<5, judge=@(x)norm(x,inf); end
    if nargin<4, eps=1e-6;              end
    if nargin<3, x0=zeros(size(b));     end
    if isempty(max_loop),  max_loop=1024;         end
    if isempty(error_gen), error_gen=@(x)x;       end
    if isempty(judge),     judge=@(x)norm(x,inf); end
    if isempty(eps),       eps=1e-6;              end
    if isempty(x0),        x0=zeros(size(b));     end
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

    counter = 1;     % actual number of iterations used
    x0 = zeros(n,1); % we need a zero matrix to build in
    x  = zeros(n,1); % initialise output matrix x

   % check for conditions of GaussSeidel, see if the matrix is strictly
   % diagonally dominant or not
    for i = 1:n
        s = 0;
        for j = 1:n
            if i ~= j 
                s = s+abs(A(i,j));
            end
        end
        
        % invalid when A(1,1)(2,2)...(n,n)>sum of the related row
        % OR any diagnal entry is zero which means A is not a strictly
        % diagonally dominant matrix or a positive definite matrix
        if s < abs(A(i,j))
            fprintf('The conditions of Guass-Seidel have not met\nPlease use other methods to find solution matrix x\n');
            return;
        end
    end

    while counter < max_loop % loop ends when exceed max no. of iterations 
        % Gauss-Seidel Iteration
        for i = 1:n
            I = [1:i-1 i+1:n];
            x(i) = error_gen( (b(i)-A(i,I)*x(I))/A(i,i) );
        end
        % calculate error and compare with eps entered
        % break if good enough
        if judge(x-x0) < eps
           break;
        end
        x0 = x;              % assign the new x to x0
        counter = counter+1; % no. of iterations
    end

    disp(['The number of iterations is ', num2str(counter)]);
    if counter >= max_loop
        fprintf('Maximum number of iterations has exceeded\n');
    end
end
