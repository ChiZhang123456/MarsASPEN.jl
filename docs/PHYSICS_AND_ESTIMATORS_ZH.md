# MarsASPEN 物理模型与三维估计量

本文档说明 MarsASPEN 中粒子源、随机碰撞、能量和电荷态更新、宏粒子权重、
三维密度、通量、反应率、Ly-alpha 发射率和能量沉积率的定义。所有公式均
对应当前 Julia 实现。

## 1. 坐标系和计算区域

模型使用 Mars Solar Orbital，MSO，坐标系：

* \(+X\) 从火星指向太阳。
* \(+Z\) 指向火星北极方向。
* \(+Y\) 完成右手坐标系。

球坐标由火星中心 Cartesian 坐标计算：

$$
r=\sqrt{x^2+y^2+z^2},
$$

$$
\lambda=\operatorname{atan2}(y,x),
\qquad
\phi=\sin^{-1}\left(\frac{z}{r}\right),
$$

$$
h=\frac{r}{1000}-R_{\mathrm{Mars}},
\qquad
R_{\mathrm{Mars}}=3388.25\ \mathrm{km}.
$$

经度输出为 \(0^\circ\) 至 \(360^\circ\)，纬度为
\(-90^\circ\) 至 \(90^\circ\)。

标准物理边界为：

$$
h_{\min}=80\ \mathrm{km},
\qquad
h_{\max}=600\ \mathrm{km}.
$$

粒子从 600 km 日侧半球注入。粒子离开 600 km 上边界、越过 80 km
下边界，或者能量低于 10 eV 时结束。默认不限制实际碰撞次数。

## 2. 日侧均匀源

### 2.1 位置采样

在半径

$$
r_{\mathrm{inj}}
=
\left(R_{\mathrm{Mars}}+600\ \mathrm{km}\right)1000
$$

的日侧球面上，对

$$
\mu=\cos(\mathrm{SZA})\sim U(0,1),
\qquad
\varphi\sim U(0,2\pi)
$$

进行采样。初始位置为

$$
x=r_{\mathrm{inj}}\mu,
$$

$$
y=r_{\mathrm{inj}}\sqrt{1-\mu^2}\cos\varphi,
$$

$$
z=r_{\mathrm{inj}}\sqrt{1-\mu^2}\sin\varphi.
$$

这种采样保证每单位球面面积的采样概率相同。所有初始位置满足
\(x\geq0\) 和 \(0\leq\mathrm{SZA}\leq90^\circ\)。

### 2.2 速度采样

标准太阳风速度分布为三维漂移 Maxwellian：

$$
f(\boldsymbol v)
\propto
\exp\left[
-\frac{m\left|\boldsymbol v-\boldsymbol U\right|^2}
{2kT}
\right],
$$

其中

$$
\boldsymbol U=(-400,0,0)\ \mathrm{km\ s^{-1}},
\qquad
kT=10\ \mathrm{eV}.
$$

速度在全局 MSO 坐标中采样。模型不会把
\(\boldsymbol U\) 旋转到每个位置的局地径向方向。局地径向速度逐粒子计算：

$$
V_r
=
\frac{\boldsymbol v\cdot\boldsymbol r}{|\boldsymbol r|}.
$$

因此同一全局速度在不同 SZA 位置具有不同的 \(V_r\)。

## 3. 宏粒子权重

### 3.1 Importance weight

物理速度分布为 \(f\)，实际抽样分布为 \(f_s\)。每个样本的无量纲
importance weight 为：

$$
w_i
=
\frac{f(\boldsymbol v_i;\boldsymbol U,T)}
{f_s(\boldsymbol v_i;\boldsymbol U,T_s)}.
$$

当 `sampling_temperature_factor=1` 时，\(T_s=T\)，因此所有
\(w_i=1\)。

### 3.2 Density weight

宏粒子的物理密度权重为：

$$
W_{n,i}
=
n_{\mathrm{source}}
\frac{w_i}{\sum_{p=1}^{N_{\mathrm{MC}}}w_p}.
$$

单位为：

$$
[W_{n,i}]=\mathrm{m}^{-3}.
$$

并且：

$$
\sum_i W_{n,i}=n_{\mathrm{source}}.
$$

\(W_n\) 本身不乘 \(V_r\)，也不乘总速度
\(|\boldsymbol v|\)。速度只在计算通量或稳态粒子注入率时使用。

### 3.3 三维稳态粒子率

三维 residence time 和 event rate 估计量需要一个宏粒子注入率。代码使用
日侧球面面积：

$$
A_{\mathrm{day}}=2\pi r_{\mathrm{inj}}^2,
$$

并计算：

$$
\dot N_i
=
A_{\mathrm{day}}W_{n,i}|V_{r,i}|.
$$

单位为 s\(^{-1}\)。这是三维稳态累计使用的派生量，不会改变或覆盖
\(W_{n,i}\)。

需要注意，这一归一化把抽样的源面密度转换为穿越日侧球面的粒子率。若将来
改用体源、有限面积平面源或指定微分通量源，需要同时重新定义
\(\dot N_i\)，不能只修改初始位置。

## 4. 自由飞行和碰撞概率

### 4.1 碰撞系数

对于位置 \(\boldsymbol r\)、能量 \(E\) 和当前电荷态，局地总碰撞系数为：

$$
\alpha(\boldsymbol r,E)
=
\sum_j n_j(\boldsymbol r)
\sum_k\sigma_{j,k}(E).
$$

其中 \(j\) 表示 CO2、O、N2，\(k\) 表示 state change、target
ionization、Ly-alpha 和 elastic。单位为：

$$
[\alpha]=\mathrm{m}^{-1}.
$$

### 4.2 Optical depth sampling

每个自由飞行段开始时抽取：

$$
U_\tau\sim U(0,1),
\qquad
\tau_{\mathrm{collision}}=-\ln U_\tau.
$$

沿路径累计：

$$
\tau(s)=\int_0^s\alpha(\boldsymbol r(s'),E(s'))\,ds'.
$$

当

$$
\tau(s)\geq\tau_{\mathrm{collision}}
$$

时发生一次碰撞。每次碰撞以后重新抽取新的
\(\tau_{\mathrm{collision}}\)。

### 4.3 自适应步长

一步的候选路径长度同时满足：

$$
\Delta s
\leq
\frac{\mathrm{safety\ factor}}{\alpha}
$$

和

$$
\Delta s\leq\Delta s_{\max}.
$$

标准三维 example 使用 \(\Delta s_{\max}=1000\) m。为了避免把一段轨迹
跨过多个高度网格，代码要求 `max_step_m` 不大于最小高度 bin 宽度。

## 5. 反应选择

碰撞发生后，目标和反应通道按照相对贡献抽样：

$$
P(j,k\mid\mathrm{collision})
=
\frac{n_j\sigma_{j,k}}
{\sum_{j'}\sum_{k'}n_{j'}\sigma_{j',k'}}.
$$

内部反应编号为：

| 编号 | 反应 |
|---:|---|
| 1 | State change |
| 2 | Atmospheric target ionization |
| 3 | H Ly-alpha production |
| 4 | Elastic collision |

对于 H+，state change 代表 charge exchange，H+ 变为 H-ENA。对于
H-ENA，state change 代表 electron stripping，H-ENA 变为 H+。

## 6. 散射角和速度更新

极角由 inverse CDF lookup table 抽样：

$$
U_\theta\sim U(0,1),
\qquad
\theta_{\mathrm{LAB}}=\theta(U_\theta).
$$

方位角独立抽样：

$$
\varphi_{\mathrm{scat}}\sim U(0,2\pi).
$$

因此即使碰撞能量和反应通道相同，散射方向仍然是随机的。

对于非弹性反应：

$$
E_{\mathrm{after}}
=
\max(E_{\mathrm{before}}-\Delta E_{j,k},0).
$$

对于弹性碰撞，代码使用实验室坐标系中的二体碰撞运动学，根据目标质量和
\(\theta_{\mathrm{LAB}}\) 计算碰撞后的速度模长。每次碰撞以后立即更新：

1. 能量。
2. 速度模长。
3. 速度方向。
4. 电荷态。
5. 下一段自由飞行的碰撞截面和碰撞系数。

## 7. 三维网格体积

经度、纬度和高度网格的精确球壳体积为：

$$
V_{\mathrm{cell}}
=
\frac{r_2^3-r_1^3}{3}
\Delta\lambda
\left[
\sin(\phi_2)-\sin(\phi_1)
\right].
$$

这里 \(\Delta\lambda\) 使用弧度，\(r_1\) 和 \(r_2\) 使用 m。因此
\(V_{\mathrm{cell}}\) 的单位为 m\(^3\)。

## 8. Residence time number density

粒子在一个网格内的驻留时间为：

$$
\Delta t_i=\frac{\Delta s_i}{|\boldsymbol v_i|}.
$$

数密度估计量为：

$$
n
=
\frac{1}{V_{\mathrm{cell}}}
\sum_i\dot N_i\Delta t_i.
$$

单位为 m\(^{-3}\)。H-ENA 和 H+ 根据自由飞行段当前的电荷态分别累计。

## 9. Flux estimators

### 9.1 Total scalar flux

$$
F_{\mathrm{total}}
=
\frac{1}{V_{\mathrm{cell}}}
\sum_i\dot N_i\Delta s_i.
$$

单位为 m\(^{-2}\) s\(^{-1}\)。该量使用路径长度，不区分方向。

### 9.2 Signed radial flux

$$
F_r
=
\frac{1}{V_{\mathrm{cell}}}
\sum_i\dot N_iV_{r,i}\Delta t_i.
$$

正值表示净向外，负值表示净向内。

### 9.3 Upward and downward radial flux

$$
F_{\mathrm{up}}
=
\frac{1}{V_{\mathrm{cell}}}
\sum_i\dot N_i\max(V_{r,i},0)\Delta t_i,
$$

$$
F_{\mathrm{down}}
=
\frac{1}{V_{\mathrm{cell}}}
\sum_i\dot N_i\max(-V_{r,i},0)\Delta t_i.
$$

它们满足：

$$
F_r=F_{\mathrm{up}}-F_{\mathrm{down}}.
$$

## 10. 反应率和电离率

一次实际 Monte Carlo 碰撞事件对所在网格的贡献为：

$$
q_{\mathrm{event},i}
=
\frac{\dot N_i}{V_{\mathrm{cell}}}.
$$

代码根据碰撞前 projectile charge state、target 和 reaction channel 分开
累计。目标 \(j\) 的总电离率为：

$$
q_{\mathrm{ion},j}
=
\sum_{\mathrm{H},\mathrm{H^+}}
q_{\mathrm{event}}
\quad
\text{for ionization events on target }j.
$$

单位为 m\(^{-3}\) s\(^{-1}\)。

`raw_monte_carlo_event_count` 只记录模拟事件个数，不乘宏粒子率，不除以
网格体积，因此不能直接解释为物理电离率。

## 11. H Ly-alpha volume emission rate

每个 Ly-alpha channel 事件贡献一个光子发射事件率。局地体积发射率为：

$$
\epsilon_{\mathrm{Ly\alpha}}
=
\sum_{\mathrm{Ly\alpha\ events}}
\frac{\dot N_i}{V_{\mathrm{cell}}}.
$$

单位为 photons m\(^{-3}\) s\(^{-1}\)。

该量不是 Rayleigh。若沿 line of sight 坐标 \(s\) 积分，并假设各向同性、
optically thin emission，则：

$$
I_{\mathrm{R}}
=
10^{-10}
\int\epsilon_{\mathrm{Ly\alpha}}(s)\,ds,
$$

其中 \(s\) 使用 m，结果单位为 Rayleigh。共振散射、吸收和多次散射需要
独立的辐射转移处理，当前 transport kernel 不包含这些过程。

## 12. 能量沉积率

碰撞能量转移为：

$$
Q_{\mathrm{collision}}
=
\frac{1}{V_{\mathrm{cell}}}
\sum_{\mathrm{events}}
\dot N_i
\Delta E_i
q_e,
$$

其中 \(q_e=1.602176634\times10^{-19}\) J eV\(^{-1}\)。单位为
W m\(^{-3}\)。

当粒子能量低于 10 eV 时，剩余能量作为 thermalization term 沉积在最终
网格：

$$
Q_{\mathrm{thermal}}
=
\frac{\dot N_iE_iq_e}{V_{\mathrm{cell}}}.
$$

总能量沉积率为：

$$
Q_{\mathrm{total}}
=
Q_{\mathrm{collision}}+Q_{\mathrm{thermal}}.
$$

## 13. 终止条件

`ParticleSummary.stop_code` 定义如下：

| Code | 条件 |
|---:|---|
| 1 | 能量低于 `min_energy_ev` |
| 2 | 越过下边界 |
| 3 | 越过上边界 |
| 4 | 触发数值 step safety limit |
| 5 | 达到用户显式设置的 `max_collisions` |

默认 `max_collisions=nothing`，因此标准模拟不应出现 code 5。生产运行结束后
必须检查 `stop_counts`。code 4 或 code 5 非零通常表示结果需要进一步诊断。

## 14. 解释密度和转换概率时的注意事项

H-ENA 密度比例不能直接等同于单次 charge exchange probability。密度是
residence time estimator，受以下因素共同控制：

1. H+ 转换为 H-ENA 的概率。
2. H-ENA 再次 stripping 为 H+ 的概率。
3. 散射以后向上或向下运动的方向。
4. 每个状态在网格中的驻留时间。
5. 速度变化。
6. 上边界逃逸和低能热化。

例如，600 km 附近的 H-ENA 可以主要来自低空 charge exchange 后向上逃逸
的粒子，而不是来自初始注入源。判断来源时应同时比较 H-ENA 上行和下行
radial flux。
