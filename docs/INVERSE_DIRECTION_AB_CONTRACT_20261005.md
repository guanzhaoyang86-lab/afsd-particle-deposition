# 同父态逆矩阵方向 A/B 数值诊断（2026-10-05）

状态：脚本/合同准备；真实 A/B 为 NOT_RUN。只有 root 排期后显式执行本合同的诊断。
该诊断不接受过程状态、不推进正式时间，不修复旧时间细化失败，不代表沉积或实验验证。
旧冻结源码及失败输出不改写；PySPH、物理方程、材料、邻域、热和接触模型保持原范围。

## 固定范围和来源

唯一父态是 `outputs/expanded_hot_budget_prefix_20261005/window_001/fine/step_012/accepted.npz`，
SHA256 `12ea12433454492d3157d3efd03c6d6178934611c3a484f6a3d2647d76659f3b`；
918 粒子、1961 自由位移 DOF，独立扩域热态时间 62.5 μs / 全局 14 步。
父态自己的旧 implementation/controller map 先按原值核验，随后新诊断 source map 才增加新文件。
完整 14 条能量历史、质量、持久 ID、初始资格、被动标签及独立起源均继承，不能接旧模型时间。
父态 readback SHA 是 `d30f8f1f30e58a5c66a6dd6d8069c482210d2980b2884ba4b58762f23b403373`。
新同父时间研究旧 FAIL 的 readback SHA 是 `89f999f398a9f3de70535500121bbba3a359ab3b39cbd351a0e8b00960e88658`；
其所有证据继续绑定，任何 A/B PASS 都不能改判它。

每根固定 dt=7.8125 μs，从同一真实父 NPZ 分别重载原序粒子及所有历史，生成完全相同初猜。
候选终时约 70.3125 μs / 候选全局步 15；该终时不得作为 accepted 时间。
顺序固定 original → inverse，只各一个完整根，原策略允许的一次同父同初猜 Newton 回退保留。
不自动换 dt、阈值、父态或追加根。

## 唯一搜索变更和所有原资格

`visibility_inverse_equilibrium.py` 保留显式 A 与 H 的既有 Broyden 更新。
方向先算 H rhs，逐 rhs 计算实际 A d − rhs，仍要求无穷范数 ≤1e−9 max(1,‖rhs‖∞)。
只可尝试一次严格改善实际余项的校正，仍失败则原 `np.linalg.solve` 一次回退并再次检查。
原 condition 估计上限 1e12、逆更新分母/探针守卫、graph 变化丢弃割线、完整残差 Armijo、
24 更新 / 100 evaluate / 至多两次 refresh、全部真实 fresh 接受规则均保留。
H 只用于搜索；不改力、位移约束、图、材料、温度或物理接受余项。

新 `inverse_direction_ab_candidate.py` 是冻结 `expanded_hot_continuation_candidate.py` 的静态复制，
只有显式局部 solver 参数和 detached diagnostic metadata 更新。preflight 以 AST 严格核对此边界，
不运行时修改旧 globals。旧 fresh、完整累计能量、材料、scalar/batch 完整映射、
有限几何、SHOULDER+TOP 域、continuous shell、retained-support whole-path 资格全执行。
每根自己的完整原门槛必须先通过；旧几何和路径查询预算均不变。

所有 checkpoint 均名 `NOT_RESTARTABLE_candidate.npz`，`accepted=false`、`restart_allowed=false`、
`diagnostic_only=true`，直接父 SHA 指向上述真实父，独立 branch origin 不变；不写 accepted.npz。
采用 `expanded_hot_continuation_acceptance.accepted_copy` 只是 fresh 数值资格，不是过程提交。

## 预登记 A/B 数值比较

所有持久 ID/顺序、immutable 字段、图 i/j 和 raw i/j、图 fingerprint、接触布尔/整数分支及
strict crossing IDs 必须一致。两条候选前 14 条历史必须逐项等于真父历史。
以下为相同 dt 下的算法误差筛选，比过程时间筛选严格，不是实验精度：

| 量 | 门槛 |
|---|---:|
| 每粒子位置差的最大 Euclidean 范数 | 1 nm |
| 最大温差 | 1e−4 K |
| 每粒子应力张量 Frobenius 差 / max(原根张量范数,1 MPa) | 0.001，即 0.1% |
| 最大累计塑性应变绝对差 | 1e−7 |
| 最大密度相对差，分母 max(原根密度,1 kg/m³) | 1e−6 |
| 接触历史向量最大差 | 1 nm |
| 每个顶层能量 `_J` 标量差 | 1e−12 J + 0.001 max(双方绝对值) |

完整粒子 properties、完整 full-map 树（含 contact/material/thermal/graph 数组及其非数组字段）
逐叶记录实际差值，observer 分别保存；不能只输出筛选量，也不能声明 bit 等价。
位置/应力/ep/温度每粒子差保存 NPZ。其他字段的实际差及单独原门槛需随结果评估，
PASS 仅表示以上事前数值筛选通过，不表示所有输出零差或可忽略所有未设门槛字段。

## 初始真实切线及成本

两根前只做一次独立真父初猜完整 map 与切线，保存同一 A/H/rhs、free mask、ID 和初猜。
对完全同一只读 A/H/rhs 分别测原 dense+余项与 verified H+condition+余项时间；
两者各自必须通过原线性余项门槛，记录方向差，不宣称机器位等价。
这份矩阵不被任一根复用，根各自从真父独立建立切线。
初始切线不是非线性平衡解，也不作为新接受态。

已测完整候选最大 28.01225140000315 s（原 campaign fine 第 1 步，含一次 Newton 回退）。
固定预测 `30 + 1.5 × 2 × (28.01225140000315 + 1)` =117.03675420000945 s <120 s <450 s。
每根增加 1 s，30 s 覆盖额外初始切线 probe、输出和前后来源检查；这是有界经验预测，非保证。
两根总内层截止 525 s、watchdog 540 s、外部监督 600 s；每根启动前再次检剩余成本。
初始 probe 和 A/B 用单 BLAS/OpenMP 线程。每个方向的校正/回退、每根迭代/evaluate/build、
完整 map/初始切线/整体候选成本都保存。单个顺序 A/B 比率不等同统计稳健或生产过程加速。

入口：

```powershell
.\.envs\pysph\python.exe -B scripts/run_inverse_direction_ab.py --output outputs/inverse_direction_ab_20261005/v1 --preflight-only
# 下行仅在 root 排期后执行：
.\.envs\pysph\python.exe -B scripts/run_inverse_direction_ab.py --output outputs/inverse_direction_ab_20261005/v1
```

15 项小矩阵/制造残差资格的 receipt 固定为
`outputs/verified_inverse_direction_20261005/algebra_v1/receipt.json`，SHA
`9912b64a09c5ce7432ac07d2c40558ad9033291a4189a33c1b060191a597232f`。
这些测试没有真实 AFSD evaluate 或 NNPS；最大 192×192 测速不可替代本 1961 DOF A/B。
源、参数、父态、旧失败、资格证据全在实际运行前写 protocol 并做运行后 hash 复查。
任何资格失败均留完整诊断并停止，整个目录始终不提供可重启过程输出。
