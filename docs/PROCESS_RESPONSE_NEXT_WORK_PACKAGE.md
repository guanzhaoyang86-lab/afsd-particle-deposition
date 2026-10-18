# 下一工作包：真实342粒子完整候选响应切片

2026-10-03。**DESIGN_ONLY / 全部下述新接口测试与四次评估 NOT_RUN。**
本次只落盘本文件，不实现新驱动、不运行物理、不改现有模型/容限。
不再扩展27粒子wall基准范围；不声称该切片能完成过程尺度输运。

## 1. 已有进展与本包边界

本轮 `outputs/process_response_20261003/saved_trial005_initial/report.json`
已经实际得到 `PASS_PRIVATE_NUMERICAL_ROOT_NOT_EVENT_OR_PROCESS`：4次评估，
内部3.1268056s、外部3.5592058s，accepted=false、restart_allowed=false。
这是27粒子冻结trial005的PRIVATE数值求根结果，不是认证停滑、接受新时序、
真实342粒子接口接通或AFSD进展。其来源、初FAIL和拒绝/PRIVATE语义保持不变。

下一包只增加一个真实构形候选响应驱动及最小回归，复用现有
`ProcessThermalTrial.evaluate`。不复制材料/接触求解器，不把新HHT固定法向
Jacobian直接嵌入不同的事件平均残差，不修改默认真实求解入口。
新增价值是把正塑性、多面接触、执行器和BE热的**同一完整候选响应**导出，
为后续实际342点过程求解提供可测对象，而非只再跑一个已有基点。

本包没有新接受步、没有自由粒子轨迹驱动、没有插粒或几何位移域扩展。
最新真实有效状态仍为850.393632910904ns，crossings=0；旧公开历史不回退。

## 2. 固定输入与选择依据

所有路径相对项目根；未来运行必须核对各自原manifest、accepted_history及hash：

| 对象 | 确切现存路径 / 用途 |
|---|---|
| 原真实342点接受parent | `outputs/coupled_afsd_20261002/short_07us/plastic_accepted.npz`，step62，t0=6.503936329109041e-7s；非零S/ep |
| parent来源与参数 | 同目录`manifest.json`及`resolved_config.json`；经`scripts/probe_updated_stage.py::inputs()`现有校验/恢复 |
| 固定完整candidate C | `outputs/nonlinear_continuation_20261003/continuation/stage_001/accepted_trial.json.gz`中的`candidate`五字段；同时保留其`evaluation`仅作历史对照 |
| C曾经获接受的凭据 | `outputs/nonlinear_continuation_20261003/continuation/manifest.json`、`accepted_history.json`第一项及`stage_001/accepted.npz`、`stage_001/acceptance.json` |
| 最新状态不回退凭据 | 同一`accepted_history.json`第四项及`stage_004/accepted.npz`；本包不得修改该history或将重放结果写入它 |
| 既有完整评估成本 | `outputs/nonlinear_continuation_20261003/cost_probe/manifest.json`，complete_trial_wall_seconds=4.04492670000036s |

固定dt=5e-8s，重放范围650.393633→700.393633ns，只是历史支路诊断，
不是从最新850ns继续，也不把这一段再次累计为物理时间。
原342点几何、质量/模量、FKS、摩擦、轴向控制、温度和全部安全域保持原值。

选择stage_001是因为它已有真实的三机制同时活跃证据，不依靠人为加载得到：
PID269新增塑性功6.484598586830537e-9J、实际率1.0478688727805792/s；
104接触、7581接受接触子区间；raw摩擦5.241676700179344e-9J；
轴向执行器功-5.595200466721313e-10J。其余后续三个接受阶段无新增塑性。
这些是选择历史切片的证据，**不是本次BE重放必须被调到相同的独立物理解**。

## 3. 复用函数与同一history接口

- `probe_updated_stage.inputs()`：复用受证据约束的原parent加载；
  `output.restore()`、`updated_trial.particle_hash()`：身份/完整史快照。
- `process_thermal_trial.ProcessThermalTrial(pa,p,t0)`及其`evaluate(C,dt)`：
  已包含当前构形NNPS/corrected伴随算子、同S0/ep0材料、真实多面event-aware
  接触/轴向约束、塑性/摩擦源和BE导热，不重新实现。
- `axial_actuator.control()`：取得实际中点free_axes及prescribed_velocity；
  `stage_iteration.FIELDS/residual_norms()`：沿用五场定义和原残差尺度。
- `probe_updated_stage.compressed()`、`run_coupled_afsd.clean/write_json()`：
  复用原始候选证据和严格JSON；新增报告不得借用旧accepted标签。

输入C字段：`position_m, velocity_m_s, density_kg_m3, temperature_K,
mean_velocity_m_s`。每次评估都从同一parent的x0/v0/rho0/T0/S0/ep0及原ID开始。
不得在两次评估间把map回写parent或累计material/thermal工作数组。

输出沿用原evaluation的`map/residual/material/internal_force_N/boundary/
thermal/midpoint/spatial/feature_guards/normal_before/normal_after`，额外包裹
本次source/config/parent/candidate hash、调用序号、完整成本和检测状态。
所有候选强制`accepted=false,restart_allowed=false,can_publish=false`。

## 4. 固定四次评估及唯一方向

一次进程最多4次完整`evaluate`，顺序 **C、C+、C-、C**。无全局求根、
线搜索、epsilon细化或追加方向。C使用JSON中的原值并核五场形状/有限性。

运行前冻结方向构造，不能按力、塑化结果或标签选择：

```text
X = parent当前LAB位置，按persistent ID升序确定归约顺序
xc = sum(mi Xi)/sum(mi)
ell = max_a(max_i Xi,a - min_i Xi,a)，必须>0
A = [[1,1/4,1/8], [1/4,-1/2,1/16], [1/8,1/16,-1/2]]
P_i = 实际中点axial_actuator.control(...).free_axes[i]
D_raw,i = P_i elementwise (A @ (Xi-xc)/ell)
D_i = D_raw,i / max_j norm(D_raw,j)，分母必须>0
epsilon = 1e-7 m/s

C±.mean_velocity = C.mean_velocity ± epsilon D
C±.velocity      = C.velocity      ± 2 epsilon D
C±.position      = C.position      ± dt epsilon D
C±.density       = C.density
C±.temperature   = C.temperature
```

该方向只用于未接受候选的数值响应测量，不是自由材料的边界速度场；约束轴
扰动为零，tracer/ASRS不进入D。参数/工具程序/时刻固定，原字段不裁剪。
报告 `(R(C+)-R(C-))/(2epsilon)`及每面/塑性/热账对应差分，保留单位和原分子。
这是一个完整残差方向响应测量，不是全Jacobian或其准确性认证；不覆盖所有
温度/密度方向，也不能以有限差分自身证明新解析导数正确。

记录实际contact_id/feature_id、`result.events`、各`interval.kind`、
`material.ep>parent.ep`、实际应变率与FKS floor分支及feature guards。
若发生分支/域变化，保存原始结果并标`NONSMOOTH_NOT_CERTIFIED`或
`REJECTED_DOMAIN`；不得靠改epsilon/调参数把它改判光滑导数PASS。

## 5. 从原字段独立重建的账目

检查器只用数组代数和声明的SI常数，不调用被测contact/material函数生成
参考答案。基点C的历史evaluation用于定位/回归，不冒称独立连续体参考。
对未收敛候选，以下账目使用**真实map/子区间结果**，不能假称C本身已平衡。

| 对账 | 原字段与独立计算 |
|---|---|
| 总冲量 | `boundary.traces[].result.velocity/normal_impulse/tangential_impulse/constraint_impulse`，核m(vmap-v0)=dt*internal_force_N+sumJn+sumJt+Jconstraint |
| 执行器 | `trace.free_axes/prescribed_velocity`与`result.constraint_impulse`；自由轴反力应为舍入零，规定轴mapped velocity等于命令；逐interval独立计算约束冲量与中点速度内积，对`constraint_work_J` |
| 原始摩擦热 | `trace.contact_rows[].normal/wall_velocity`、`result.intervals[].v_before/v_after/tangential_impulse`；先独立重建相对中点切向速度，再核逐face Q=-Jt·u_mid及原`friction_heat_raw_J`/`friction_rounding_bound_J`，不能用总和掩盖负值 |
| 同史塑性 | `material.ep/deviatoric_stress_Pa/plastic_work_J`、parent.ep及`midpoint.density_kg_m3`；独立q_new=sqrt(1.5*S:S)，P=(m/rhomid)*q_new*(ep1-ep0)；检查单调ep和正P，return_dissipation单列 |
| 热源与BE热 | `map.temperature_K`、`thermal.particle_energy_J`、`thermal.pair_conductance_W_K`、`midpoint.i/j`、`thermal_boundary_certificate.wall_row_W_K/wall_weighted_temperature_W`；重建pair两端相反热、dt*(HTwall-H*Tmap)、chi*P及eta*Q，核mcpΔT逐粒子/总量 |
| 重入与身份 | 第1/4调用map、material、boundary、thermal数值数组与ID逐位一致（时钟字段除外）；parent全字段/状态hash前后完全相同；不累计ep/热/冲量 |

沿用现有已注册界，不建立更宽界：逐候选动量每分量m_i*1e-9m/s+1e-17Ns；
自由轴反力1e-17Ns、规定速度1e-12m/s；接触功1e-15J+1e-9*相应绝对功尺度；
rawQ身份1e-25J+1e-12相对，负Q仅使用原逐事件舍入证书；BE逐粒子
1e-12J+1e-10*sum(abs(热源))、总热1e-12J+1e-9*同尺度，pair净热
1e-15J+1e-10*sum(abs(pair热))。塑性代数身份用1e-18J+1e-10相对，
是SI浮点公式重建界而非材料准确性界；不得修改材料返回结果以过该界。
完整字段重入/ID/hash要求精确相同。新检查器实际测试全部NOT_RUN。

候选全局残差单独输出原五场norm，不要求通过求根门，也不得因此发表物理态。
不得把`certify_process_thermal(C,...)`对候选T的检查替换成map检查后仍沿用其
“候选通过”标签；本包仅宣告映射内部账目及未提交响应，绝不是整体接受。

## 6. 时间预算与事务约束

旧完整试探4.0449267s，仅作为初始成本依据，不是本次BE接口实测。
预测 `1.5*4*4.0449267+5 = 29.26956s`，登记为约30s。
一次调用内部50s（含加载/归档/检查）、外部60s。每次调用后重算：

```text
elapsed + 1.5*remaining_calls*max(4.0449267, 已测完整单次最大秒) + 5 <= 50
```

否则剩余调用NOT_RUN并正常收束；任何单次卡住由内部watchdog和外部timeout终止。
不得以4次评估为名循环全局根，或拆分多个60s进程规避总预算。
保存每次完整输入/输出后再开始下一次；异常/超时未完成态明确REJECTED。

运行前与结束后必须核同一source/config/父状态/原archive清单；源清单包括
实际加载闭包、PySPH pin、未来驱动/检查器/本合同及所用现有材料/接触/热/空间
文件，不仅顶层脚本。`inputs()`返回的旧sources须与原凭据核对，再显式补入
`process_thermal_trial.py`、`implicit_thermal.py`等本次实际调用文件。
结束sourcehash漂移、父数组变化、原归档变化均FAIL，不能接受“开始时一致”。

输出必须新目录、原子写入临时文件后发表诊断report；不得覆盖旧
accepted_history/checkpoint/manifest，不调用`continuation_state.commit`或
HHT migration/physical checkpoint发布。新NPZ如需保存只能PRIVATE并拒绝restart。
最新850ns的接受史只读绑定，不因650ns重放改变latest指针。

## 7. 未来实现范围与完成判定

建议新增 `scripts/probe_afsd_process_bridge.py`、必要的纯数组读回检查器、
`tests/review/test_afsd_process_bridge.py`；均为计划，当前NOT_IMPLEMENTED/NOT_RUN。
先对输入来源、固定方向/约束轴、计数/预算、同史不累加、拒绝态和hash漂移做
小接口回归，再执行唯一四评估进程。这里的小接口测试不是新27粒子物理wall基准。

本包交付：真实342点的可重跑完整候选接口、四次或明确部分完成的响应/成本、
独立原始账目与严格非提交证据；允许结论为PASS_RESPONSE_MEASUREMENT、
NONSMOOTH_NOT_CERTIFIED、BLOCKED_BUDGET或明确FAIL，不追求全部变绿。
若没有接触、正塑性或执行器实际活动，不得声称“三机制已覆盖”。

当前工程近似授权允许这一同物理纯候选测量，无需新的物理选择。但它不授权
删除惯性、清零850ns速度、改变接触势/材料参数、取消有限feature保护或50um域界，
不授权GPU、缩放、裁剪、长跑/付费或三工况生产。真实旋转有限边界的整路径
表示、一般3D粘着/执行器事件、宏步准确性/成本和输运仍须后续实际证明。

**状态：只完成本设计文档；本包新代码、回归、四次评估和输运全部NOT_RUN。**
