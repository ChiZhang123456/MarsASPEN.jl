# MarsASPEN.jl 代码结构与使用说明

本文档说明 MarsASPEN.jl 中每个主要代码文件的职责、输入输出、调用关系、
物理单位和常用运行方式。它同时说明哪些输出是原始 Monte Carlo 统计量，
哪些输出可以进一步转换为物理通量、数密度、体积发射率或亮度。

## 1. 总体计算流程

一次模拟依次执行以下步骤：

1. `load_model` 读取一组 MGITM 冷大气、对应的 MAMPS hot O，以及 H 和 H⁺
   碰撞截面。
2. `MonteCarloConfig` 定义粒子数、初始高度、初始速度、随机种子、终止条件
   和宏粒子权重。
3. `run_particle_core` 为每个粒子产生独立随机数流，并从 600 km 开始输运。
4. 每个自由飞行段计算局地中性密度、总碰撞系数和到下一次碰撞的距离。
5. 碰撞发生时，根据各通道的 `n × sigma` 权重选择靶粒子和反应。
6. 立即更新能量、速度方向和电荷态，然后重新计算下一段的碰撞概率。
7. ensemble driver 将单粒子结果汇总为 summary、反应高度 histogram、
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
* `MonteCarloConfig`：定义一次模拟的全部粒子和数值参数。
* `ParticleSummary`：每个粒子的最终能量、高度、碰撞数和终止原因。
* `HistoryEvent`：详细轨迹模式下的一行事件记录。

当前默认粒子为 600 km 处、速度 400 km s⁻¹、朝向 Mars 的 neutral H ENA。
对应初始能量约为 835 eV。

`particle_weight` 默认为 1。它是宏粒子权重，不是自动计算得到的物理通量。

### 3.3 `src/initialization.jl`

负责模型初始化。

`normalize_f107` 将 `solar_min`、`solar_moderate`、`solar_max` 或数字形式统一为
70、130、200。

`available_atmosphere_cases` 返回 12 组可用的 `Ls × F10.7` 条件。

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

```text
alpha = sum(n_target × sigma_target,reaction)
```

`alpha` 的单位是 m⁻¹，表示单位路径长度的碰撞概率。

`choose_event` 使用所有通道的 `n × sigma` 作为权重，随机选择靶粒子和反应。

### 3.6 `src/transport.jl`

这是单粒子输运核心。

对每个粒子，先抽取下一次碰撞的光学深度：

```text
tau_collision = -log(U)
```

然后沿轨迹累计：

```text
tau = integral(alpha ds)
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
* `path_length_m`：累计高度、能量和电荷态中的 `particle_weight × ds`。

### 3.7 `src/ensembles.jl`

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

每段贡献为 `particle_weight × ds`。使用路径长度而不是 step count，可以避免
自适应步长造成的统计偏差。

`run_detailed_ensemble` 保存完整轨迹，只适合少量粒子。

为了避免线程锁，histogram 模式为每个 Julia thread 分配独立矩阵，模拟结束
后再统一相加。

### 3.8 `src/io.jl`

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

## 4. `scripts` Julia 运行脚本

### 4.1 `scripts/benchmark.jl`

运行 compact ensemble 并报告速度、平均最终能量、平均碰撞数、平均 step 数
和各终止类型数量。

```powershell
julia --project=. -t auto scripts/benchmark.jl 1000000
```

### 4.2 `scripts/run_detailed.jl`

为少量粒子保存完整轨迹 MAT 文件。

```powershell
julia --project=. -t auto scripts/run_detailed.jl 10 output/detailed.mat
```

不要用此模式保存 100 万粒子，因为完整 step history 可能达到数百 GB。

### 4.3 `scripts/run_reaction_altitude_counts.jl`

使用低内存 histogram 统计不同高度的四类反应事件数。

```powershell
julia --project=. -t auto scripts/run_reaction_altitude_counts.jl 1000000 output/counts.csv 1
```

最后一个参数是高度 bin 宽度，单位 km。

### 4.4 `scripts/run_phase_space_histogram.jl`

统计 H ENA 和 H⁺ 的高度能量分布。

```powershell
julia --project=. -t auto scripts/run_phase_space_histogram.jl 1000000 output/phase.mat 1 100 1
```

最后三个参数分别为：

1. 高度 bin 宽度，km。
2. 对数能量 bin 数。
3. 初始宏粒子权重。

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

### 5.4 `analysis/scripts/plot_reaction_altitude_counts.py`

绘制 state change、ionization、Ly-alpha 或 elastic count 随高度的分布。

### 5.5 `analysis/scripts/plot_phase_space_histogram.py`

读取 phase-space MAT 文件，绘制 H ENA 和 H⁺ 的高度能量二维分布。颜色表示：

```text
sum(particle_weight × ds)
```

### 5.6 `analysis/scripts/plot_detailed_mat.py`

这是详细轨迹绘图的命令行入口。

### 5.7 `analysis/tests/test_io.py`

测试 MAT v7.3 reader、单粒子筛选和反应事件筛选。

## 6. 原始统计量与物理量

必须区分以下三类量：

1. Reaction count：反应事件数。
2. Path-length histogram：`sum(w_i ds)`。
3. Physical rate or brightness：需要通量权重、时间、面积、体积或视线积分。

若入射 H 通量为 `F_H`，模拟面积为 `A`，统计时间为 `Delta t`，模拟粒子数为
`N`，可以使用：

```text
w_i = F_H A Delta t / N
```

一维高度模型中的体积反应率还需要除以高度层厚度和对应体积。Ly-alpha
brightness 需要从体积发射率沿观测视线积分，并转换为 Rayleigh。

## 7. 修改代码后的验证

建议每次修改物理或数值算法后执行：

```powershell
julia --project=. -t auto test/runtests.jl
julia --project=. -t auto scripts/benchmark.jl 1000
```

需要重点检查：

* 固定 seed 的结果是否保持可重复。
* 密度和截面单位是否一致。
* 能量变化后是否重新计算截面。
* charge state 是否在 state-change 反应后正确更新。
* histogram 是 event count 还是 path-length estimator。
* 大规模运行是否触发 collision 或 step safety limit。
