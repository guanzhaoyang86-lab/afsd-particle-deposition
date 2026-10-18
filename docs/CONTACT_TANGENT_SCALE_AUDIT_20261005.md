# 固定 1 nm 接触差分尺度审计（2026-10-05）

结论：在本次保存的单一候选上，未发现 1 nm 扰动使任何有自由轴的粒子跨越
active、friction_active、sliding 或 normal_derivative_limit_branch。
差分尺度确实造成有限的接触块近似差，但证据不足以把它认作频繁 Broyden/Newton 失败的主因。
不开发新切线变体、不改变原 eps，也不将结论推广至其他候选或接触事件。

## 输入和计算边界

- 真父：`outputs/expanded_hot_budget_prefix_20261005/window_001/fine/step_012/accepted.npz`，
  SHA `12ea12433454492d3157d3efd03c6d6178934611c3a484f6a3d2647d76659f3b`，62.5 μs / 14 步。
- 实际候选：`outputs/contact_onset_time_study_20261005/v1/eight/step_001/FULL_UNCOMMITTED_candidate.npz`；
  SHA `dba364866cdb396d4ef093c219a51c997e06c40201eadb3229a4df6c0d3b9f94`，dt=0.9765625 μs。
  同目录 `NOT_RESTARTABLE_candidate.npz` 的 SHA 为
  `469a4ad6dcc499a3fa9c7c5a9c2bcce9d4f9fe43547dd205296edfd5923401d6`。
  918 点，最大增量位移 10.2973 nm；该候选不作为已接受父态。
- 原 `quasistatic_process.py::tangent` 对 x/y/z 分别施加 ±1e−9 m，六次接触映射，
  把中心差分的 3×3 局部块加入材料/空间近似切线。该切线本来就不是完整全场 Jacobian。
- 本审计仅重建完全相同的 `BatchedChainRuleContact`，原 x0、实际 x1、工具位姿、dt 和接触历史不变。
  `check_path=False` 与原 tangent 扰动相同，不能把它当连续接触或实体穿透资格。
  未调用 NNPS、材料、热、trial.evaluate、tangent 或求根；没有过程接受/时间推进。

正式 v3 基线接触的全部数组与实际保存字段逐项相同，最大差为 0；输入数组未被修改。
随后只比较预登记的 h=1e−9、1e−10、1e−11、1e−12 m，每个 h 的三个轴 ± 扰动。

## 分支和导数结果

四档 h 下，z 扰动均使 256 个恰在 BASE 壳边界的固定底面粒子切换 active/friction_active。
三个轴扰动也改变这些固定粒子的 derivative-limit 分支。它们的三个轴均固定，
局部接触块在抽取自由行/自由列后全部移除，不能将这些切换当作自由方棒求根的病因。
任一有自由轴粒子的上述四种分支，在 ± 两侧及相对基线的变化数都为 0。

对真实自由行/自由列投影后的接触块，邻接差分尺度结果为：

| h 粗→细 | 整体块 Frobenius 相对差 | 最大粒子块差 N/m | 最大逐点相对差 |
|---|---:|---:|---:|
| 1 nm →0.1 nm | 4.68235e−4 | 291.4678 | 2.36095e−3 |
| 0.1 nm →0.01 nm | 4.67682e−6 | 2.910778 | 2.35773e−5 |
| 0.01 nm →0.001 nm | 5.66323e−8 | 0.0331227 | 2.68294e−7 |

逐点分母为 max(较细块范数,1 N/m)，三行最大均为 PID 796。
其 1 nm 与 0.001 nm 块差为 294.4117 N/m，即约 0.23847%；
PID 892 对应约 0.12636%。缩小 h 的变化趋势与平滑分支上的中心差分收敛相容，
但四档数据不构成精确导数或全场 Newton 收敛证明。
原时间筛选相关 PID 797 在本候选接触力为 0，四档局部接触块均为 0；
这同样不能代表其后续接触事件或之前失败搜索点。

未保存的完整失败迭代、其他起始状态、材料/密度/热反馈和空间算子未在此审计中覆盖。
因此这里不能判定唯一根因，不能证明更小 h 能让过程通过时间筛选，也不能改变任何旧 FAIL。

## 可复核证据与实际预算

- 正式结果：`outputs/contact_tangent_scale_20261005/v3/summary.json`。
- 全部扰动力、法向/切向合力、分支标志和 3×3 块：同目录 `contact_scale_arrays_NOT_RESTARTABLE.npz`。
- 纯保存数组自由轴投影读回：`outputs/contact_tangent_scale_20261005/free_projection_readback.json`。
- 正式接触-only脚本：`scripts/audit_contact_tangent_scale_v3.py`，SHA
  `d9731fb337491c248735e39a4c169febb20f83fc7faa9ab5895e17470b1d100c`。
- 外部 receipt：`logs/contact_tangent_scale_20261005/v3/receipt.json`，exit 0，内部 0.6638 s，外部 0.8979 s。

两个诊断记录错误原样保留：v1 在布尔数组差值摘要处 TypeError（1 次接触调用、外部 0.5989 s）；
v2 完成 25 次接触调用后，在 `set.update(int)` 汇总处 TypeError（外部 0.7969 s）。
v3 只修正这两项记录操作，未改接触或其他冻结物理源。
三次合计 **51 次 contact-only 调用、外部约 2.2937 s**，低于本次 30 s 总预算；
NNPS、材料、热、完整物理 evaluate 和求根调用均为 0。v1/v2 的 FAIL 不改判。

该支线至此停止，原过程继续使用冻结接触切线。
