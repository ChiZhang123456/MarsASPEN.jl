# MarsASPEN MAT 输出格式

本文档说明详细轨迹 MAT 文件和三维空间网格 MAT 文件的变量、维度、单位和
读取方法。

## 1. 两类输出

MarsASPEN 支持两种输出模式：

1. Detailed trajectory output，保存少量粒子的每一步和每次碰撞。
2. Spatial grid output，直接累计大量粒子的三维物理量，不保存单粒子历史。

详细轨迹适合检查算法和制作单粒子图。三维网格适合 \(10^5\) 至
\(10^7\) 粒子的生产模拟。

## 2. Detailed trajectory MAT

`write_detailed_mat` 写出压缩 MAT v7.3 文件，格式版本为：

```text
aspen_julia_history_v1
```

### 2.1 变长轨迹索引

所有粒子的 event rows 被连接为一维表。以下变量恢复每个粒子的范围：

| Variable | Meaning |
|---|---|
| `particle_start_index_1based` | 每个粒子第一行的 Julia one based index |
| `particle_event_count` | 每个粒子的 event row 数量 |
| `particle_id` | 每行所属的 particle ID |

Python reader 会把 one based index 转换为 NumPy index。

### 2.2 Event type

| Code | Name | Meaning |
|---:|---|---|
| 0 | `initial` | 初始状态 |
| 1 | `transport` | 自由飞行 step 结束 |
| 2 | `collision` | 碰撞和状态更新 |
| 3 | `final` | 最终状态 |

### 2.3 每一行的主要变量

| Variable | Unit | Meaning |
|---|---|---|
| `time_s` | s | 从注入开始的累计时间 |
| `x_m`, `y_m`, `z_m` | m | MSO Cartesian position |
| `vx_m_s`, `vy_m_s`, `vz_m_s` | m s\(^{-1}\) | 当前速度 |
| `altitude_km` | km | 相对 \(R_{\mathrm{Mars}}\) 的高度 |
| `energy_before_ev` | eV | 碰撞前能量 |
| `energy_ev` | eV | 当前或碰撞后能量 |
| `vx_before_m_s`, etc. | m s\(^{-1}\) | 碰撞前速度 |
| `charge_state` | integer | 0 为 H-ENA，1 为 H+ |
| `target_code` | integer | 0 none，1 CO2，2 O，3 N2 |
| `reaction_code` | integer | 0 none，1 state change，2 ionization，3 Ly-alpha，4 elastic |
| `energy_loss_ev` | eV | 本次碰撞 projectile energy change |

### 2.4 每个粒子的 summary

| Variable | Meaning |
|---|---|
| `final_energy_ev` | 最终能量 |
| `final_altitude_km` | 最终高度 |
| `n_steps` | transport step 数量 |
| `n_collisions` | 碰撞数量 |
| `stop_code` | 终止原因 |

## 3. 三维空间输出

`write_spatial_grid_mats` 将结果分成四个文件，避免一个 MAT 文件过大：

```text
<prefix>_grid.mat
<prefix>_moments.mat
<prefix>_reactions.mat
<prefix>_energy.mat
```

格式版本为：

```text
marsaspen_spatial_grid_v1
```

### 3.1 Dimension order

Julia 内部和文件的逻辑维度顺序为：

```text
longitude, latitude, altitude, optional components
```

MAT v7.3 底层使用 HDF5。某些 Python HDF5 reader 会反转维度。项目提供的
`load_history_mat` 会根据坐标长度恢复正确的 NumPy 维度。不要仅根据
HDF5 dataset 的原始 shape 猜测维度含义。

### 3.2 公共坐标和 metadata

四个文件都包含：

| Variable | Unit | Meaning |
|---|---|---|
| `longitude_edges_deg` | deg | \(0^\circ\) 至 \(360^\circ\) 经度边界 |
| `latitude_edges_deg` | deg | \(-90^\circ\) 至 \(90^\circ\) 纬度边界 |
| `altitude_edges_km` | km | 高度边界 |
| `longitude_centers_deg` | deg | 经度中心 |
| `latitude_centers_deg` | deg | 纬度中心 |
| `altitude_centers_km` | km | 高度中心 |
| `coordinate_system` | text | `MSO` |
| `n_particles` | count | Monte Carlo 粒子数 |
| `seed` | integer | 随机种子 |
| `stop_counts` | count | 五种停止原因的数量 |
| `total_collisions` | count | 全部轨迹碰撞数 |
| `total_steps` | count | 全部 transport step 数 |

生产 example 还写入 source density、初始速度、初始高度、上下边界和注入面积。

## 4. Grid file

`<prefix>_grid.mat` 包含：

| Variable | Shape | Unit |
|---|---|---|
| `cell_volume_m3` | latitude, altitude | m\(^3\) |

同一纬度和高度的所有经度 cell 体积相同，因此经度维没有重复存储。

## 5. Moments file

Charge dimension order 为：

```text
0: H_ENA
1: Hplus
```

| Variable | Shape | Unit |
|---|---|---|
| `number_density_by_charge_m3` | charge, altitude, latitude, longitude in normalized Python output | m\(^{-3}\) |
| `total_number_density_m3` | altitude, latitude, longitude | m\(^{-3}\) |
| `total_flux_by_charge_m2_s` | charge plus spatial dimensions | m\(^{-2}\) s\(^{-1}\) |
| `total_flux_m2_s` | spatial dimensions | m\(^{-2}\) s\(^{-1}\) |
| `signed_radial_flux_by_charge_m2_s` | charge plus spatial dimensions | m\(^{-2}\) s\(^{-1}\) |
| `signed_radial_flux_m2_s` | spatial dimensions | m\(^{-2}\) s\(^{-1}\) |
| `upward_radial_flux_by_charge_m2_s` | charge plus spatial dimensions | m\(^{-2}\) s\(^{-1}\) |
| `downward_radial_flux_by_charge_m2_s` | charge plus spatial dimensions | m\(^{-2}\) s\(^{-1}\) |

应满足：

$$
F_r=F_{\mathrm{up}}-F_{\mathrm{down}}.
$$

`total_flux` 是 scalar track length flux，不等于
\(F_{\mathrm{up}}+F_{\mathrm{down}}\)，因为粒子速度可能包含水平分量。

## 6. Reactions file

Target order：

```text
0: CO2
1: O
2: N2
```

Reaction order：

```text
0: state_change
1: ionization
2: Ly_alpha
3: elastic
```

| Variable | Meaning | Unit |
|---|---|---|
| `reaction_rate_m3_s1` | charge, target, channel resolved rate | m\(^{-3}\) s\(^{-1}\) |
| `reaction_rate_by_channel_m3_s1` | charge and target summed | m\(^{-3}\) s\(^{-1}\) |
| `ionization_rate_by_target_m3_s1` | charge summed, target resolved | m\(^{-3}\) s\(^{-1}\) |
| `total_lya_volume_emission_rate_photons_m3_s1` | all charges and targets summed | photons m\(^{-3}\) s\(^{-1}\) |
| `raw_monte_carlo_event_count` | raw channel event count | count |

`raw_monte_carlo_event_count` 只用于收敛性检查，不能当作物理反应率。

## 7. Energy file

| Variable | Meaning | Unit |
|---|---|---|
| `collision_energy_transfer_by_charge_w_m3` | collision loss by charge | W m\(^{-3}\) |
| `collision_energy_transfer_w_m3` | charge summed collision loss | W m\(^{-3}\) |
| `cutoff_thermalization_by_charge_w_m3` | sub 10 eV remaining energy | W m\(^{-3}\) |
| `cutoff_thermalization_w_m3` | charge summed cutoff energy | W m\(^{-3}\) |
| `total_energy_transfer_w_m3` | collision plus cutoff thermalization | W m\(^{-3}\) |

## 8. Python 读取示例

```python
from pathlib import Path
import numpy as np

from marsaspen_analysis import load_history_mat

prefix = Path("examples/output/dayside_hplus_10000000_3d")
moments = load_history_mat(Path(f"{prefix}_moments.mat"))

altitude = np.ravel(moments["altitude_centers_km"])
density = np.asarray(moments["number_density_by_charge_m3"])

hena_density = density[0]
hplus_density = density[1]
```

## 9. 全球球面平均

纬度 cell 的面积权重为：

$$
w_\phi
=
\sin(\phi_{\mathrm{upper}})
-\sin(\phi_{\mathrm{lower}}).
$$

若经度 bin 等宽，则全球球面平均为：

$$
\langle X\rangle(h)
=
\frac{\sum_{\lambda,\phi}X(\lambda,\phi,h)w_\phi}
{\sum_{\lambda,\phi}w_\phi}.
$$

不能直接对纬度 index 做普通平均，因为高纬度 cell 面积较小。

## 10. 文件大小和 Git

一千万粒子的三维 MAT 文件可能达到数十 MB。`examples/output/` 默认不纳入
Git。GitHub 保存绘图代码和 PNG 结果，不保存每次生产运行的全部大型 MAT。

若需要长期归档，应保存：

1. 四个 MAT 文件。
2. Git commit hash。
3. Julia version 和 thread count。
4. 运行命令。
5. 完整 metadata。
6. stdout 和 stderr log。
