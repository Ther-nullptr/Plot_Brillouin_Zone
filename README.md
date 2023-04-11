# 固体物理大作业——绘制布里渊区图景

[TOC]

## 0 代码组织

```bash
.
├── plot_brillouin_zone.m # 绘制布里渊区
├── plot_field.m # 绘制势场
├── plot_fourier.m # 绘制傅里叶系数
├── plot_gap.m # 计算带隙
├── plot_near_brillouin.m # 绘制在原点附近的布里渊区
├── plot_near_brillouin_with_near_free_electron.m # 使用近自由电子法绘制布里渊区
└── README.md # 文档
```



## 1 画出势能分布曲线

势能分布为：

$$
V(x) = \begin{cases}
1\times10^{-19}\cdot\cos\frac{2\pi}{a}x, -\frac{\pi}{2}\le\frac{2\pi}{a}x\le\frac{\pi}{2} \\
0, elsewhere
\end{cases}
$$

编写脚本`plot_field.m`:

```matlab
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
```

结果如下：

![image-20230410205319356](https://s2.loli.net/2023/04/10/qicQlTCp1WnwYXo.png)

## 2 至少画出前4条能带并对比带隙和能带曲线

要利用特征根法求解布里渊区曲线，需要求以下方程的解：

$$
\det\begin{bmatrix}
\cdots & \cdots & \cdots & \cdots & \cdots \\
\cdots & \frac{\hbar^2}{2m_0}(k-\frac{2\pi}{a})^2+V_0-E & V_{-1} & V_{-2} & \cdots \\
\cdots & V_1 &  \frac{\hbar^2}{2m_0}k^2+V_0-E & V_{-1} & \cdots \\
\cdots & V_2 &  V_1 &  \frac{\hbar^2}{2m_0}(k+\frac{2\pi}{a})^2+V_0-E & \cdots \\
\cdots & \cdots & \cdots & \cdots & \cdots \\
\end{bmatrix} = 0
$$

给定k，可以得到N个E的解，N的数目由周期性势场的傅里叶展开决定。

对此，我们先对势场进行傅里叶展开，得到各傅里叶系数（令$N=10$）：

$$
V_n = \frac{1}{P}\int_P \Re{V(x)}\cdot e^{-i\frac{2\pi nx}{P}} dx \\
=\frac{1}{P}\int_P V(x) \cdot e^{-i\frac{2\pi nx}{P}} dx
$$

```matlab
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
```

结果如下：

![image-20230410212516771](https://s2.loli.net/2023/04/10/eV3PSFE6hbLKYrR.png)

将得到的Fourier级数储存，构造求解矩阵：

```matlab
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
k_vec = -pi/a:pi/(a*100):pi/a;
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
figure(2)
hold on
for n = 1:N + 1
    plot(k_vec, eigen_vec(n,:));
end

xlabel('k');
ylabel('E/J');
```

> 如何快速从对角方向构建矩阵：
>
> ```matlab
> mat = zeros(N + 1, N + 1);
> 
> for n = -N:1:N
>     	series = ones(1, N + 1 - abs(n)) * n;
>     	mat = mat + diag(series, -n);
> end
> 
> disp(mat);
> ```
>
> ```bash
> >> test
>  0    -1    -2    -3    -4    -5    -6    -7    -8
>  1     0    -1    -2    -3    -4    -5    -6    -7
>  2     1     0    -1    -2    -3    -4    -5    -6
>  3     2     1     0    -1    -2    -3    -4    -5
>  4     3     2     1     0    -1    -2    -3    -4
>  5     4     3     2     1     0    -1    -2    -3
>  6     5     4     3     2     1     0    -1    -2
>  7     6     5     4     3     2     1     0    -1
>  8     7     6     5     4     3     2     1     0
> ```

布里渊区绘景如下：

![image-20230410225353842](https://s2.loli.net/2023/04/10/PosH3nIRKAYdWVL.png)

要计算能带宽度，即计算以下间隙的宽度：

![image-20230410232225338](https://s2.loli.net/2023/04/11/cLX5Dbgr8eUIvB6.png)

对此计算 $k=0$ 和 $k=\frac{\pi}{a}$ 处的能量值：

```matlab
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
k_vec = 0:pi/a:pi/a;
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

middle_gap = eigen_vec(:,1);
side_gap = eigen_vec(:,2);
middle_gap_diff = diff(middle_gap);
side_gap_diff = diff(side_gap);
gap_vec = zeros(1, N);

for i = 1:N
    if (mod(i, 2) == 0)
        gap_vec(i) = side_gap_diff(i-1);
    else
        gap_vec(i) = middle_gap_diff(i+1);
    end
end

disp(gap_vec)
```

计算结果如下：

```bash
>> plot_gap
   1.0e-19 *

    0.4934    0.2271    0.0072    0.0416    0.0008    0.0179    0.0010    0.0099
```

根据理论推导可知，带隙的大小为 $2|V_n|$ 。将理论计算值和实际计算值进行对比：

![image-20230411092503473](https://s2.loli.net/2023/04/11/XuafeKM6ygv9tj2.png)

```bash
>> plot_gap
   1.0e-19 *

    0.4934    0.2271    0.0072    0.0416    0.0008    0.0179    0.0010    0.0099

   1.0e-19 *

    0.5000    0.2121    0.0000    0.0423    0.0000    0.0181    0.0000    0.0100
```

可见，近自由电子得到的带隙与实际带隙相差较小。

首先绘制出布里渊附近$\Delta k = \pm \frac{1}{10}\cdot \frac{2\pi}{a}$的能带曲线：

![image-20230411093621903](https://s2.loli.net/2023/04/11/4ZDFkaiV17AlxH8.png)

为了更好地显示近似效果，我们取用第一和第二第三条能带的曲线：

![image-20230411101628283](https://s2.loli.net/2023/04/11/T5p18vEbqIRZ9Xn.png)

![image-20230411101643210](https://s2.loli.net/2023/04/11/JHXcylrBzb8hTYN.png)

可见能量分布近似于抛物线。

近自由电子下，零级解：

$$
E_k^0 = \frac{\hbar^2 k^2}{2m_0} +\bar V
$$

本征值的零级修正：

$$
E_k^{(0)} = \frac{\hbar^2k^2}{2m_0}
$$

本征值的一级修正：

$$
E_k^{(1)} = 0
$$

本征值的二级修正：

$$
E_k^{(2)} = \sum_{n\ne0} \frac{|V_n|^2}{\frac{\hbar^2}{2m_0}[k^2 - (k+\frac{2\pi}{a}n)^2]}
$$

则：

$$
E = \bar V + E_k^{(0)} + E_k^{(1)} + E_k^{(2)}
$$

使用近似比较在 $\Delta k = \pm \frac{1}{10} \frac{2\pi}{a}$ 处的效果：

![image-20230411141640283.png](https://s2.loli.net/2023/04/11/E8RywQ2fALK3oW6.png)

可见拟合效果良好。
