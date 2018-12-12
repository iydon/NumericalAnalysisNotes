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