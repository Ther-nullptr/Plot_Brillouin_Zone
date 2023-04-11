a = 5.43 * 10^(-10);
m_0 = 9.1 * 10^(-31);
h = 6.63 * 10^(-34);

x = -a/2:a/100:a/2;
V = 10^(-19) * cos(2*pi/a*x) .* (x < a/4 & x > -a/4);

% plot x and V
plot(x, V);
xlabel('x/m');
ylabel('V(x)/J');
title('potential Energy of half cosine field');

