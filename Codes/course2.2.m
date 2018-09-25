f = @(x) cos(x)-x;
Bisection(f, [0,1])
g = @(x) cos(x);
fixed_point(g, 0)
h = @(x) exp(x)-7;
Newton_method(h, 0)

function mid = Bisection(f, interval, max_step, epsilon)
    if nargin<4, epsilon=1e-6; end
    if nargin<3, max_step=128; end
    mid_last = interval(1);
    if prod(f(interval))<0
        for i=1:max_step
            mid = mean(interval);
            if abs(mid-mid_last)<epsilon || abs(f(mid))<epsilon
                disp(['Step: ', num2str(i)]);
                disp(['Zero: ', num2str(mid)]);
                return
            else
                if prod(f([mid,interval(1)]))<0
                    interval(2) = mid;
                else
                    interval(1) = mid;
                end
            end
            mid_last = mid;
        end
    end
end

function new_val = fixed_point(f, start, max_step, epsilon)
    if nargin<4, epsilon=1e-6; end
    if nargin<3, max_step=128; end
    new_val = f(start);
    for i=1:max_step
        old_val = new_val;
        new_val = f(old_val);
        if abs(old_val-new_val)<epsilon
            disp(['Step: ', num2str(i)]);
            disp(['Fixed-point: ', num2str(new_val)]);
            return
        end
    end
end

function val = Newton_method(f, start, max_step, epsilon)
    if nargin<4, epsilon=1e-6; end
    if nargin<3, max_step=128; end
    df = matlabFunction(diff(sym(f)));
    newton_fun = @(x) x - f(x)./df(x);
    val = fixed_point(newton_fun, start, max_step, epsilon);
end