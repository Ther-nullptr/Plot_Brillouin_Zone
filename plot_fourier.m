a = 5.43 * 10^(-10);
m_0 = 9.1 * 10^(-31);
h = 6.63 * 10^(-34);

x = -a/2:a/100:a/2;
V = @(x) 10^(-19) * cos(2*pi/a*x) .* (x < a/4 & x > -a/4);

N = 8;

% find the fourier series coefficients of V
a_n = zeros(1,2*N + 1);

for n = -N:1:N
    a_n(n + N + 1) = real(1/a * integral(@(x) exp(-1i*2*pi*n*x/a) .* V(x), -a/2, a/2));
end

% plot the fourier series coefficients
figure(1);
stem(-N:N, a_n, 'o');
title('a_n');
xlabel('n');
ylabel('V_n');