# MarsASPEN.jl 代码结构与使用说明

本文档说明 MarsASPEN.jl 中每个主要代码文件的职责、输入输出、调用关系、
物理单位和常用运行方式。它同时说明哪些输出是原始 Monte Carlo 统计量，
哪些输出可以进一步转换为物理通量、数密度、体积发射率或亮度。

## 1. 总体计算流程

一次模拟依次执行以下步骤：

1. `load_model` 读取一组 MGITM 冷大气、对应的 MAMPS hot O，以及 H 和 H⁺
   碰撞截面。
2. `MonteCarloConfig` 定义粒子数、初始高度、初始速度、随机种子和终止条件。
3. `MonteCarloWeight` 独立定义采样温度、源区数密度和宏粒子权重。
4. `run_particle_core` 为每个粒子产生独立随机数流，并从 600 km 开始输运。
5. 每个自由飞行段计算局地中性密度、总碰撞系数和到下一次碰撞的距离。
6. 碰撞发生时，根据各通道的 $n_j\sigma_{j,k}$ 权重选择靶粒子和反应。
7. 立即更新能量、速度方向和电荷态，然后重新计算下一段的碰撞概率。
8. ensemble driver 将单粒子结果汇总为 summary、反应高度 histogram、
   高度能量 histogram，或详细逐步轨迹。

核心调用关系为：

```text
load_model
  ├─ load_atmospheres
  └─ load_cross_sections

run_ensemble / run_binned_ensemble / run_phase_space_ensemble
  └─ run_particle_core
       ├─ local_state
       │    ├─ transport_density3
       │    ├─ hot_o_density
       │    └─ sigma_at
       ├─ choose_event
       └─ rotate_velocity
```

## 2. 单位和编号约定

内部计算尽量使用 SI 单位：

| 变量 | 单位 |
|---|---|
| Cartesian position | m |
| Altitude interface | km |
| Velocity | m s⁻¹ |
| Energy | eV |
| Number density | m⁻³ |
| Cross section | m² |
| Temperature | K |
| Path length | m |

电荷态编号：

| Code | Projectile |
|---:|---|
| 0 | neutral H ENA |
| 1 | H⁺ |

靶粒子编号：

| Code | Target |
|---:|---|
| 1 | CO₂ |
| 2 | O |
| 3 | N₂ |

反应编号：

| Code | Reaction |
|---:|---|
| 1 | state change，H 与 H⁺ 之间转换 |
| 2 | target ionization |
| 3 | H Ly-alpha production |
| 4 | elastic collision |

终止编号：

| Code | Meaning |
|---:|---|
| 1 | projectile energy below `min_energy_ev` |
| 2 | altitude below lower boundary |
| 3 | altitude above upper boundary |
| 4 | numerical step safety limit |
| 5 | collision safety limit |

## 3. `src` 核心 Julia 文件

### 3.1 `src/MarsASPEN.jl`

这是 package 入口文件。它只负责：

1. 定义 `MarsASPEN` module。
2. 加载 MAT、Random、Statistics 和 Threads。
3. export 用户需要调用的公开接口。
4. 按依赖顺序 include 其他源文件。

此文件不包含具体物理算法。

### 3.2 `src/types.jl`

保存物理常数、数据结构和运行配置。

主要结构：

* `Atmosphere`：MGITM 经度、纬度、高度、五种组分的 log density 和温度。
* `HotAtmosphere`：MAMPS hot O 网格和 log density。
* `CrossSections`：统一能量网格、截面、固定能量损失和散射角表。
* `AspenModel`：将冷大气、hot O 和碰撞数据库组合为完整模型。
* `MonteCarloConfig`：仅定义粒子初态、数值参数和终止条件。
* `MonteCarloWeight`：定义重要性采样、源区密度和宏粒子权重。
* `ParticleSummary`：每个粒子的最终能量、高度、碰撞数和终止原因。
* `HistoryEvent`：详细轨迹模式下的一行事件记录。

当前默认粒子为 600 km 处、速度 400 km s⁻¹、朝向 Mars 的 neutral H ENA。
对应初始能量约为 835 eV。

权重不再保存在 `MonteCarloConfig` 中。无物理源密度时，
`MonteCarloWeight().unit_particle_weight` 默认为 1。

### 3.3 `src/initialization.jl`

负责模型初始化。

`normalize_f107` 将 `solar_min`、`solar_moderate`、`solar_max` 或数字形式统一为
70、130、200。

`available_atmosphere_cases` 返回 12 组可用的 $L_s\times\mathrm{F10.7}$ 条件。

`load_atmospheres` 读取：

* `gitm_lsXXX_fXXX.mat`
* `mamps_lsXXX_fXXX.mat`

密度在读入时转换为 log density，保证空间插值和向下外推始终为正值。

`load_model` 是用户主要入口。它调用 `load_atmospheres` 和
`load_cross_sections`，最后返回一个 `AspenModel`。

### 3.4 `src/atmosphere.jl`

负责所有中性大气插值和外推。

经度使用周期边界，纬度使用线性插值，密度在 log density 空间插值。

MGITM 原始高度范围为 98.75 至 251.25 km。

* 98.75 至 80 km：每个经纬度柱使用最低两层进行 log-linear 外推。
* 80 km 以下：固定为 80 km 的值。
* 251.25 km 以上：使用局地顶层温度和各组分质量进行静力指数外推。

MAMPS hot O 仅在原始高度范围内插值。其最低高度是 100 km，因此 100 km
以下 hot O 为 0。

主要函数：

* `density3`：返回 CO₂、O、O₂、N₂、CO。
* `transport_density3`：只返回有碰撞截面的 CO₂、O、N₂。
* `hot_o_density`：返回 MAMPS hot O。
* `temperature3`：返回 MGITM 中性温度。
* `neutral_density`：公开的 named tuple 接口。
* `neutral_density_xyz`：从 Cartesian position 查询中性大气。

### 3.5 `src/cross_sections.jl`

负责碰撞截面读取、能量插值和反应选择。

原始截面为 cm²，读入后统一乘以 `1e-4` 转换成 m²。

内部截面数组为：

```text
sigma[charge + 1, target, reaction, energy]
```

`load_cross_sections` 将 H 和 H⁺ 的不同原始列顺序转换为统一反应顺序：

```text
state change, ionization, Ly-alpha, elastic
```

`local_state` 计算局地总碰撞系数：

```math
\alpha
=
\sum_j n_j
\sum_k \sigma_{j,k}.
```

$\alpha$ 的单位是 m$^{-1}$，表示单位路径长度的碰撞概率。

`choose_event` 使用所有通道的 $n_j\sigma_{j,k}$ 作为权重，随机选择靶粒子和反应。

### 3.6 `src/monte_carlo_weight.jl`

这个文件负责漂移 Maxwellian 的重要性采样和物理宏粒子权重。物理速度
分布为 $f(\boldsymbol v;\boldsymbol U,T)$，采样时可以使用更宽的分布
$f_s(\boldsymbol v;\boldsymbol U,T_s)$。每个
随机粒子的无量纲重要性权重为：

```math
W_i
=
\frac{
f(\boldsymbol v_i;\boldsymbol U,T)
}{
f_s(\boldsymbol v_i;\boldsymbol U,T_s)
}.
```

给定源区数密度 `n_source` 后，粒子密度权重为：

```math
W_{n,i}
=
n_{\mathrm{source}}
\frac{W_i}{\sum_{p=1}^{N_{\mathrm{MC}}}W_p}.
```

因此所有粒子的密度权重之和严格等于输入数密度。初始宏粒子权重就是
$W_{n,i}$，单位是 m$^{-3}$，不乘任何速度。需要计算通量时，再根据所需
分量使用：

```math
F_{r,i}=W_{n,i}V_{r,i},
\qquad
F_{\mathrm{total},i}
=
W_{n,i}
\sqrt{V_{x,i}^{2}+V_{y,i}^{2}+V_{z,i}^{2}}.
```

$F_r$ 和 $F_{\mathrm{total}}$ 的单位均为 m$^{-2}$ s$^{-1}$。径向速度
的正负用于区分上行和下行。两个方向的粒子都保留在 Maxwellian 样本中。

标准调用方式为：

```julia
config = MonteCarloConfig(
    n_particles=1_000_000,
    initial_speed_m_s=400_000.0,
    initial_temperature_ev=10.0,
)
weighting = MonteCarloWeight(
    sampling_temperature_factor=5.0,
    source_number_density_m3=5.0e6,
)
result = run_directional_flux_ensemble(
    model, config; weighting=weighting,
)
```

### 3.7 `src/spatial_grid.jl`

该文件负责完整的三维经纬度和高度网格统计。水平网格与MGITM一致，高度网格可由调用者设置，标准example使用1 km。它根据宏粒子的真实粒子率权重和每段轨迹的驻留时间计算局地数密度，根据路径长度计算总通量，根据径向位移计算上行、下行和带符号径向通量。

碰撞发生后，该文件按照碰撞位置、入射粒子电荷态、背景目标种类和反应种类累计体积反应率。它还负责把三维结果拆分写入grid、moments、reactions和energy四个MAT文件，避免把所有大型数组集中到单一文件中。

### 3.8 `src/transport.jl`

这是单粒子输运核心。

对每个粒子，先抽取下一次碰撞的光学深度：

```math
\tau_{\mathrm{collision}}=-\ln U.
```

然后沿轨迹累计：

```math
\tau=\int\alpha\,\mathrm{d}s.
```

当累计光学深度达到抽样阈值时发生碰撞。

空间步长同时受两个条件限制：

1. `safety_factor / alpha`
2. `max_step_m`

碰撞后立即更新能量和电荷态。下一段自由飞行重新计算截面和碰撞概率。

`rotate_velocity` 在以原速度为极轴的局地坐标中应用散射角，并保持更新后的
速度模长。

`run_particle_core` 支持三种可选统计方式：

* `record=true`：保存所有 step 和 collision。
* `reaction_counts`：只累计不同高度和反应类型的事件数。
* `path_length_m`：累计高度、能量和电荷态中的
  $\mathtt{particle\_weight}\times\mathrm{d}s$。

### 3.9 `src/ensembles.jl`

负责多粒子并行运行。

每个粒子的随机数由 `(seed, particle_id)` 独立确定，因此线程调度不会改变
模拟结果。

`run_ensemble` 只返回每个粒子的 compact summary。

`run_binned_ensemble` 返回：

```text
reaction_counts[altitude_bin, reaction]
```

这是碰撞事件数，不是体积反应率。

`run_phase_space_ensemble` 返回：

```text
path_length_m[altitude_bin, energy_bin, charge_state]
```

每段贡献为 $\mathtt{particle\_weight}\times\mathrm{d}s$。使用路径长度而不是
step count，可以避免自适应步长造成的统计偏差。

`run_detailed_ensemble` 保存完整轨迹，只适合少量粒子。

为了避免线程锁，histogram 模式为每个 Julia thread 分配独立矩阵，模拟结束
后再统一相加。

### 3.10 `src/io.jl`

`write_detailed_mat` 将变长粒子历史展平为列数组，并写入压缩 MAT v7.3。

`particle_start_index_1based` 和 `particle_event_count` 用于恢复每个粒子的
事件范围。

事件类型：

| Code | Event |
|---:|---|
| 0 | initial state |
| 1 | transport step |
| 2 | collision |
| 3 | final state |

## 4. Julia 运行脚本

### 4.1 `scripts/run_detailed.jl`

为少量粒子保存完整轨迹 MAT 文件。

```powershell
julia --project=. -t auto scripts/run_detailed.jl 10 output/detailed.mat
```

不要用此模式保存 100 万粒子，因为完整 step history 可能达到数百 GB。

### 4.2 `examples/sample_dayside_injection_10000000.jl`

在 600 km 高度的完整 MSO 日侧半球上按球面面积均匀采样 10 万个初始位置。
太阳风速度在 MSO 坐标系中取 $[-400,0,0]$ km s$^{-1}$，而不是在每个位置
旋转到局地径向方向。该脚本只保存注入位置、速度、SZA 和 $W_n$，不运行
输运。

```powershell
julia --project=. -t auto examples/sample_dayside_injection_10000000.jl
```

### 4.3 `examples/run_dayside_3d_10000000.jl`

这是当前的大粒子数标准示例。它使用均匀日侧半球源运行 10 万个 H⁺
macro particles，并在经度、纬度和 1 km 高度网格中累计粒子矩、反应率和
能量转移率。

```powershell
julia --project=. -t auto examples/run_dayside_3d_10000000.jl
```

输出分成 grid、moments、reactions 和 energy 四个 MAT v7.3 文件。详细的
macro-particle 归一化、单位和变量说明见 `examples/README.md`。

## 5. Python 分析代码

### 5.1 `analysis/marsaspen_analysis/io.py`

读取普通 MAT 或 MAT v7.3，并提供：

* `load_history_mat`
* `particle_history`
* `reaction_events`

### 5.2 `analysis/marsaspen_analysis/plotting.py`

分析详细轨迹文件，输出：

* 多粒子高度、能量、速度和电荷态 overview。
* Cartesian trajectory。
* 每个粒子的高度时间图及反应标记。
* 所有碰撞事件的 CSV 表。

### 5.3 `analysis/scripts/plot_neutral_atmosphere.py`

绘制 MGITM 冷组分、MAMPS hot O、total O 和中性温度。

### 5.4 `analysis/scripts/plot_detailed_mat.py`

这是详细轨迹绘图的命令行入口。

### 5.5 Uniform-dayside 三维绘图

`examples/plot_dayside_3d_maps_120km.py` 绘制120 km经纬度分布。
`examples/plot_dayside_3d_altitude_profiles.py` 绘制100至300 km高度剖面。
`examples/plot_dayside_3d_sza_altitude.py` 绘制SZA和高度二维分布。

### 5.6 `analysis/tests/test_io.py`

测试 MAT v7.3 reader、单粒子筛选和反应事件筛选。

## 6. 原始统计量与物理量

三维均匀日侧模拟必须区分以下物理量：

1. 原始Monte Carlo event count，它只用于统计收敛性检查。
2. Residence-time estimator，它给出网格内数密度。
3. Track-length estimator，它给出total scalar flux。
4. Event-rate estimator，它给出电离率和Ly-alpha体发射率。

宏粒子首先携带与速度无关的密度权重：

```math
W_{n,i}
=
n_{\mathrm{source}}
\frac{W_i}{\sum_p W_p},
\qquad
[W_{n,i}]=\mathrm{m}^{-3}.
```

因此，$W_{n,i}$ 本身不乘 $V_r$ 或 $|\boldsymbol V|$。对于高度面通量，
代码直接根据所需方向计算 $W_{n,i}V_{r,i}$，或者根据速度模长计算
$W_{n,i}|\boldsymbol V_i|$。

三维稳态 residence-time estimator 需要把源面密度权重转换为穿过日侧球面
的粒子率。这个派生量只用于三维稳态累计：

```math
\dot N_i
=
A_{\mathrm{day}}W_{n,i}|V_{r,i}|,
\qquad
A_{\mathrm{day}}=2\pi r_{\mathrm{inj}}^2,
\qquad
[\dot N_i]=\mathrm{s}^{-1}.
```

这里使用绝对值，因此 inward 和 outward 样本都保留。$\dot N_i$ 不会写回
或改变 $W_{n,i}$。

对于网格体积 $V_{\mathrm{cell}}$、驻留时间 $\Delta t_i$ 和路径长度
$\Delta s_i=|\boldsymbol V_i|\Delta t_i$，三个主要估计量为：

```math
n
=
\frac{\sum_i\dot N_i\Delta t_i}{V_{\mathrm{cell}}},
```

```math
F_{\mathrm{total}}
=
\frac{\sum_i\dot N_i\Delta s_i}{V_{\mathrm{cell}}},
```

```math
q_{\mathrm{event}}
=
\frac{\sum_{\mathrm{events}}\dot N_i}{V_{\mathrm{cell}}}.
```

它们的单位依次为 m$^{-3}$、m$^{-2}$ s$^{-1}$ 和 m$^{-3}$ s$^{-1}$。

## 7. 修改代码后的验证

建议每次修改物理或数值算法后执行：

```powershell
julia --project=. -t auto test/runtests.jl
C:\Users\Win\.conda\envs\mars\python.exe -m pytest analysis/tests
```

需要重点检查：

* 固定 seed 的结果是否保持可重复。
* 密度和截面单位是否一致。
* 能量变化后是否重新计算截面。
* charge state 是否在 state-change 反应后正确更新。
* 三维输出使用的是event、residence-time还是track-length estimator。
* 大规模运行是否触发 collision 或 step safety limit。
