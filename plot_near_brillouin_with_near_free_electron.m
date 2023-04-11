a = 5.43 * 10^(-10);
m_0 = 9.1 * 10^(-31);
h = 6.63 * 10^(-34);
hbar = h/(2*pi);

x = -a/2:a/100:a/2;
V = @(x) 10^(-19) * cos(2*pi/a*x) .* (x < a/4 & x > -a/4);

N = 8;

% find the fourier series coefficients of V
V_n = zeros(1,2*N + 1);

for n = -N:1:N
    V_n(n + N + 1) = real(1/a * integral(@(x) exp(-1i*2*pi*n*x/a) .* V(x), -a/2, a/2));
end

% plot the brillouin zone (theory)
k_vec = [];
for i = -N/2-1:1:N/2+1
    k_vec = [k_vec, (-pi/(20*a)+pi*i/a):pi/(a*100):(pi/(20*a)+pi*i/a)];
end

diag_vec = zeros(1, N + 1);
base_mat = zeros(N + 1,N + 1);
eigen_vec = zeros(N + 1, length(k_vec));

for n = -N:1:N
    series_vec = ones(1, N + 1 - abs(n)) * V_n(n + N + 1);
    base_mat = base_mat + diag(series_vec, -n);
end

for k = 1:length(k_vec)
    diag_vec = hbar ^ 2 * (k_vec(k) + (-N/2:1:N/2)*2*pi/a).^2 / (2*m_0);
    mat = base_mat + diag(diag_vec);
    eigen_vec(:,k) = eig(mat);
end

% plot the brillouin zone (use near free electron approximation)
V_n_new = V_n;
V_n_new(N + 1) = [];
n_vec = -N:1:N;
n_vec(N + 1) = [];
V_bar = V_n(N + 1);
E_0 = hbar^2 * k_vec.^2 / (2*m_0);
E_1 = 0;
E_2 = zeros(1, length(k_vec));
for n = 1:2*N
    E_2 = E_2 + V_n_new(n)^2 ./ (hbar^2 * (k_vec.^2 - (k_vec + 2*pi*n_vec(n)/a).^2) / (2*m_0));
end

% plot the eigenvalues
figure(1);
for n = 1:N/2 + 1
    plot(k_vec, eigen_vec(n,:)); 
    hold on;
end

stem(k_vec, V_bar + E_0 + E_1 + E_2, 'r');
xlabel('k');
ylabel('E/J');
title('Brillouin Zone, |k|<pi/5a');
