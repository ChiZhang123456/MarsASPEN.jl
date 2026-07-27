# MarsASPEN 运行、收敛和验证指南

本文档给出从安装、小规模测试到一千万粒子生产运行的推荐流程，并说明每次
运行以后必须检查的数值和物理量。

## 1. 环境

推荐版本：

* Julia 1.10 或更高。
* Python 3.10 或更高。
* Windows 上的 Python：
  `C:\Users\Win\.conda\envs\mars\python.exe`

初始化 Julia：

```powershell
julia --project=. -e "using Pkg; Pkg.instantiate()"
```

安装 Python analysis：

```powershell
C:\Users\Win\.conda\envs\mars\python.exe -m pip install -e analysis
```

## 2. 修改代码以后先运行测试

Julia：

```powershell
julia --project=. -t 4 -e "using Pkg; Pkg.test()"
```

Python：

```powershell
C:\Users\Win\.conda\envs\mars\python.exe -m pytest analysis/tests
```

Julia tests 检查：

1. 80 km 向下外插。
2. 日侧位置均匀采样。
3. MSO 速度和局地径向速度。
4. 三维 density、flux、reaction 和 energy accumulator。
5. 固定 seed 的确定性。
6. H 和 H+ state change。
7. ionization 和 Ly-alpha estimator。
8. \(W_n\) 和 importance sampling normalization。

## 3. 单粒子验证

先运行 H-ENA：

```powershell
julia --project=. examples/run_h_ena_trajectory.jl
```

再运行 H+：

```powershell
julia --project=. examples/run_hplus_trajectory.jl
```

需要检查：

* Altitude 是否连续。
* Energy 是否只在碰撞时改变。
* State change 是否正确切换 0 和 1。
* Scattering angle 是否在 collision rows 中有值。
* Collision 以后是否重新抽取 optical depth。
* 终止高度和 stop code 是否合理。

## 4. 注入器验证

只采样，不输运：

```powershell
julia --project=. -t 16 examples/sample_dayside_injection_10000000.jl
```

验证标准：

* 所有 \(x\geq0\)。
* 所有位置位于 600 km 球面。
* \(\cos(\mathrm{SZA})\) 近似均匀。
* \(V_x\) 均值接近 \(-400\) km s\(^{-1}\)。
* \(V_y\) 和 \(V_z\) 均值接近 0。
* \(\sum_iW_{n,i}=5\times10^6\) m\(^{-3}\)。

## 5. 小规模三维测试

生产运行前建议依次使用：

```powershell
julia --project=. -t 8 examples/run_dayside_3d_10000000.jl examples/output/test_1000 1000 80
```

```powershell
julia --project=. -t 8 examples/run_dayside_3d_10000000.jl examples/output/test_100000 100000 80
```

检查 stdout：

```text
transport_elapsed_s
write_elapsed_s
total_steps
total_collisions
```

并读取 metadata 中的 `stop_counts`。

## 6. Stop count 审查

五个 stop code 为：

```text
1 low energy
2 lower boundary
3 upper boundary
4 numerical step limit
5 optional collision limit
```

标准配置使用 `max_collisions=nothing`，因此：

```text
stop_counts[5] == 0
```

必须成立。正常生产运行也应满足：

```text
stop_counts[4] == 0
```

如果 code 4 非零，应检查：

1. `max_step_m` 是否过小。
2. 局地密度或截面是否产生异常 \(\alpha\)。
3. 是否有速度接近 0 但能量仍高于 cutoff 的粒子。
4. 边界判断是否正确。

## 7. 一千万粒子生产运行

Windows 上建议固定线程数：

```powershell
julia --project=. -t 16 examples/run_dayside_3d_10000000.jl
```

不建议在已知高核心数环境中盲目使用 `-t auto`。过多线程会增加锁竞争和
垃圾回收压力。固定 8 至 16 线程通常更稳定，实际最优值应通过小规模 benchmark
确定。

运行条件：

```text
N = 10,000,000
altitude = 600 km
geometry = uniform dayside hemisphere
U_MSO = (-400, 0, 0) km/s
kT = 10 eV
n_source = 5 cm^-3
lower boundary = 80 km
upper boundary = 600 km
minimum energy = 10 eV
max collisions = none
```

## 8. Monte Carlo 收敛

建议至少比较：

```text
10^3, 10^4, 10^5, 10^6 particles
```

对每个粒子数比较以下剖面：

1. H+ density。
2. H-ENA density。
3. H+ upward and downward flux。
4. H-ENA upward and downward flux。
5. O ionization rate。
6. CO2 ionization rate。
7. H Ly-alpha VER。
8. Energy deposition rate。

对于 profile \(X_N(h)\)，可以用较大样本作为参考：

$$
\delta_N(h)
=
\frac{|X_N(h)-X_{\mathrm{ref}}(h)|}
{\max(|X_{\mathrm{ref}}(h)|,X_{\mathrm{floor}})}.
$$

二维 map 的空白或低计数区域应同时检查 raw event count，避免把 Monte Carlo
零计数误认为严格的物理零值。

## 9. 守恒和一致性检查

### 9.1 Weight normalization

$$
\sum_iW_{n,i}=n_{\mathrm{source}}.
$$

### 9.2 Radial flux identity

$$
F_r=F_{\mathrm{up}}-F_{\mathrm{down}}.
$$

### 9.3 Nonnegative quantities

以下量应非负：

* Number density。
* Total scalar flux。
* Upward flux。
* Downward flux。
* Reaction rate。
* Ly-alpha VER。
* Energy deposition rate。

### 9.4 Charge state interpretation

H-ENA density 高于 H+ density 不等于所有注入 H+ 都在该高度第一次发生
charge exchange。密度包含驻留时间和返回粒子。应结合上下行 flux、
state change event rate 和单粒子轨迹解释。

## 10. Atmosphere validation

每次替换 GITM 或 AMPS 数据后检查：

1. 经度为 \(0^\circ\) 至 \(360^\circ\)。
2. Subsolar point 对应 longitude \(=0^\circ\)，latitude \(=0^\circ\)。
3. SZA 满足

$$
\mathrm{SZA}
=
\cos^{-1}(\cos\phi\cos\lambda).
$$

4. 所有密度读入后单位为 m\(^{-3}\)。
5. GITM density 在 log space 插值。
6. GITM temperature 线性插值。
7. 80 km 外插连续。
8. AMPS hot O 只在 native altitude range 内使用。

## 11. Cross section validation

检查：

1. 输入单位 cm\(^2\) 是否乘 \(10^{-4}\) 转为 m\(^2\)。
2. H 和 H+ 的列顺序是否分别正确映射。
3. Energy 超出 tabulated range 时是否返回 0。
4. State change 后 charge state 是否更新。
5. Inelastic energy loss 符号是否正确。
6. Elastic energy transfer 是否使用 target mass。

## 12. 图形验证

所有图：

* 非数学文字使用 Arial。
* 单位写在 axis 或 colorbar label。
* Log colorbar 只显示正值。
* H+ 和 H-ENA 使用一致的颜色定义。
* 共用 colorbar 时两个 panel 必须使用同一个 `Normalize` object。
* 经度 map 使用 \(0^\circ\) 至 \(360^\circ\)。
* SZA 图标出 \(90^\circ\) terminator。

## 13. 生产运行记录

建议每次保存一个纯文本记录：

```text
date
git commit
Julia version
thread count
command
particle count
seed
atmosphere Ls and F10.7
source density
source temperature
source bulk velocity
boundaries
energy cutoff
stop_counts
total_steps
total_collisions
runtime
output filenames
```

只有记录完整配置，两个 Monte Carlo 结果才可以严格比较。
