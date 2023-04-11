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

% plot the brillouin zone
k_vec = -pi/(5*a):pi/(a*100):pi/(5*a);
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

% plot the eigenvalues
figure(1);
for n = 1:N + 1
    if mod(n,2) == 0
        figure(floor(n/2) +1);
    else
        hold on;
        xlabel('k');
        ylabel('E/J');
        title('Brillouin Zone, |k|<pi/5a');
    end
    plot(k_vec, eigen_vec(n,:)); 
end

