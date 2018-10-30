%{
Mid-point:      1+h¦Ë+(h¦Ë)^2/2
Modified Euler: 1+h¦Ë+(h¦Ë)^2/2
Runge-Kutta 4:  1+h¦Ë+(h¦Ë)^2/2+(h¦Ë)^3/6+(h¦Ë)^4/24
Let HL = h¦Ë.
%}
t = linspace(-5,5,256);
[X,Y] = meshgrid(t);

syms x
for i = 2:10
    f = matlabFunction(taylor(exp(x),'Order',i));
    stab_plot(X,Y,f,i);
end
%{
f = @(HL) 1+HL;
stab_plot(X,Y,f);
f = @(HL) 1+HL+HL.^2/2;
stab_plot(X,Y,f);
f = @(HL) 1+HL+HL.^2/2+HL.^3/6+HL.^4/24;
stab_plot(X,Y,f);
%}

function stab_plot(X,Y,f,i)
    figure;
    Z = abs(f(X+1i*Y));
    Z(Z>1) = nan;
    contour(X,Y,Z);
    axis equal; colorbar;
    title(['Order: ',num2str(i)])
    shading interp;
    saveas(gcf,[num2str(i-1),'.jpg'])
end