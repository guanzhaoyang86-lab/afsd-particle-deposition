# AFSD output schema v0.1.0

## 当前：率相关材料数值响应与保存journal作用（2026-10-04）

新增MATERIAL_RESPONSE组件输出不是物理更新：d_ep/d_equivalent/d_stress为一个
方向参数的导数，SI单位分别无量纲/Pa/Pa；FKS rate偏导Pa*s、T偏导Pa/K。
C_alg=dS1/dL单位Pa*s（已含材料步dt），不能再视为弹性模量重复乘dt。
数值P矩阵kg、rhs与algebraic_reaction_impulse_Ns为Ns；反力仅搜索代数，不进功热。

`outputs/material_response_20261004/saved_action_initial/report.json`状态为
COMPLETE_SAVED_ACTION_NOT_ROOT_OR_PROCESS；new_evaluations/new_afsd_steps/new_transport_events=0，
accepted/restart_allowed/can_publish/production_authorized=false。源和组件测试receipt均绑定。
`action_arrays.npz`是同父候选差dC/dR/dcorrection/E及持久ID，不是checkpoint或新时间帧。
`UNEVALUATED_BASE_PROPOSAL.json.gz`保存真实基点的数值proposal/correction/diagnostics，
明确UNEVALUATED，不得用于补帧、重启、真实输运或用其字段宣称温度/应力已求解。

比较表保存全部五场和分组，温度未改善必须保留；零输入方向不计算虚构消减比。
新P仍省略几何/体积/B、接触及BE反馈导数，NONSMOOTH原范围不取消。
原始物理输出/参数均不因后处理或此数值补丁改变，真实公开末态仍850.393633ns。

## 当前：真实342点同史响应，仅诊断（2026-10-03）

`outputs/process_bridge_20261003/real342_initial`为
`REAL_AFSD342_SAME_HISTORY_RESPONSE_MEASUREMENT_NOT_PROCESS`。
四次C/C+/C-/C完整评估已运行；status=NONSMOOTH_NOT_CERTIFIED是数值分支变化，
不是接受末态。所有accepted/restart_allowed/can_publish=false，无新checkpoint或VTK。
公开最新仍850.393632910904ns/cross0；650ns父→700ns诊断不得串接为新时序。

- `parent_state.json.gz`：完整原PA属性、持久ID和原solver state；只读历史来源。
- `configuration.json`及report：参数/方向/epsilon/dt、源码和输入hash、时间预算与实测成本。
  source hash不等于所有安装二进制的完整闭包认证。
- `direction.json`：按原ID/几何/质量/自由轴构造D；约束轴为零；保留实际浮点C+-C-。
- `base/plus/minus/base_repeat_candidate.json.gz`：五场候选及其实际evaluation，
  含材料、内力、逐face冲量、约束反力、raw摩擦热与BE热账；不累加成四个物理步。
- 各`*_audit.json`：独立map账目、负rawQ逐事件和原舍入界；未用abs/max裁负值。
  `*_branches.json`：face身份、塑性/约束/事件及自适应区间签名，不仅比较总数。
- `directional_response.json`：完整残差R=C-map的单epsilon方向商及每face响应，
  分母2epsilon（m/s）；位置商单位s，速度/均速商无量纲，密度商kg*s/m^4，温度商K*s/m。
  不是已认证Jacobian、实验数据或新材料属性。差商为负不等于原始热量为负。
- `last_accepted_reference.json`：分别绑定最新850ns公开checkpoint与650ns历史父，
  history_rewound=false。PRIVATE候选不得成为重启输入、补帧或真实输运统计。

逐位重入相同、独立账目PASS只证明本切片的确定性与数值记账；基点BE温度根残差
仍超原门，不能用map温度替换C温度后宣布接受。7个粒子的算法分支变化原样保留，
不得为可视化把NONSMOOTH改成已验证。完整范围见PROCESS_BRIDGE_REPORT及两份独立读回。

后续只读预条件作用若保存NPZ，必须标ALGEBRA_ONLY、无物理评估/接受/重启，
不得同真实时间序列合并。它评估数值搜索方向，不定义粒子速度或物理力。
本次已有action_arrays.npz保存两P的dC/dR/dcorrection/E，使用字段自身SI单位，
不是新时刻的x/v/rho/T。因果分解保存的应力差为Pa、力差为N，属于同父候选差；
正负分量不得当作新功热。各分项相消时禁止用绝对值占比拼成“能量百分比”。

## 当前：固定节点响应的PRIVATE数值诊断（2026-10-03）

`outputs/process_response_20261003/saved_trial005_initial`分类为
`FIXED_QUADRATURE_DISCRETE_CONTACT_PRECONDITIONER_EXPERIMENT_V1`。
report的PASS_PRIVATE_NUMERICAL_ROOT_NOT_EVENT_OR_PROCESS仅表示同保存区间
满足旧原始根门限。new_afsd_steps/new_published_local_steps均0；没有新accepted_history。

`last_accepted_state.json`只指向旧局部6us公开checkpoint及hash，不是新增末态。
`configuration.json`绑定canonical inputs、实际参数、原dt与原contact_id事件方向。
`PRIVATE_NONRESTARTABLE.npz`保存40属性和完整私有数值史；accepted=false、
restart_allowed=false，项目restore已实际拒绝。通用NPZ读回仅用于审计。
`PRIVATE_DIAGNOSTIC.json.gz`保存原完整残差/功热/节点及新增响应记录。
若求解失败，REJECTED_CANDIDATE.json保存失败上下文，不能转为末态。

contact_response_history的矩阵单位分别是Kx[N/m]、Kv[N*s/m]、
actual_impulse_acceleration_response_kg_s[kg*s]、endpoint_force_acceleration_response_kg[kg]、
root_contact_response_kg[kg]；matrix_norm_kg/matrix_asymmetry_norm_kg为数值搜索矩阵诊断，
numerical_matrix_condition及merit无量纲。它们不代表新物理力/热源。
固定节点响应不含adaptive节点选择导数；continuous_cross_zero_jacobian_certified=false。

即使数值根PASS和总Q为正，完整多项式内部零点及whole_path_passive=false仍禁止发布。
有限节点无负Q也不能推导整路径被动。连续尖点和私有跨零尾段须保留解释，
见PROCESS_RESPONSE_LIMITS、REPORT、OUTPUT_REVIEW。后处理不得改权限、补帧或重启。
真实AFSD有效终时仍850.393633ns、cross0；局部公开终时仍6us。

## 当前：停滑事件开发与私有候选（2026-10-03）

`HHT_LOCATED_STOP_ELASTIC_BE_STATE_V1`与
`HHT_ACTUAL_CONTACT_IMPULSE_WITH_EXPLICIT_ZERO_EVENT_V1`是显式新数值定义，
不与旧checkpoint自动混用。当前实际输出在outputs/process_stop_20261003/scout_initial。
accepted_history仅两记录：5us同物理字段的定义迁移、5→6us新接受区间；只有后一条推进时间。
无已发布的耦合零滑移/右极限记录。真实AFSD没有新增状态。

accepted/restart_allowed均为true且完整semantic/particle/parameter哈希及发布索引一致，
才是可继续的局部态。PRIVATE探根候选这两字段均false；其source、父态、duration、
残差、接触积分及有符号节点功热全部保留，不能因joint root暂时收敛便提升为有效物理态。
REJECTED JSON含pending candidate/每次回溯及失败原因，不是末态。

`event_context`按真实contact_id而非tracer索引；模式为private_continuation、left_zero、
right_zero，incoming_direction意义不变；right_zero另需力学可行性给出的outgoing_direction。
私有延拓输出physical_contact_path_certified=false，负节点Q保留；绝不abs/max。
严格端点证书保存完整滑移多项式、单侧导数、原判零界和端点舍入别名，不改任何速度。

候选保存last_interval实际Jn/Jt、材料/墙功、raw Q、delta impulse/work和算法/root功；
严格提升同时绑定父态、前后时间/dt、诊断candidate hash以及保存/重算的全部积分数组。
诊断算法功不进入热。temperature仍由真实BE热账推进。

将来的零时长event_right记录必须有独立record_kind与父hash，物理时间/全部PA属性连续，
累计物理冲量/功热不重复累计。保存左右实际力、旧canonical A、新Fplus/m和有符号
acceleration_restart_storage_change_J（仅数值储能诊断，不是材料热/完整自由能）。
输出acceleration向量以state中的canonical值为准；事件处保留的PA历史属性不被伪装成新解。
当前这一3D事件发布路径为NOT_RUN/未通过，不可从组件PASS继承认证。

独立参考的event_state NPZ包含每个测得事件前后状态与ledger；两者逐位相同。
reference_history只认证至6.5us（105共同样本），参考是独立离散ODE基准而非论文实验数据。
随后的可行静摩擦事件候选明确REJECTED，两个不同终时不得用于宣称收敛。
完整结果和单位见PROCESS_STOP_REPORT、REFERENCE_RESULTS与OUTPUT_REVIEW。

## 当前：冲量一致曲线接触积分输出（2026-10-03）

`outputs/process_curved_20261003/S0_initial`只有13发布态到5us，分类
LOCAL_CURVED_CONTACT_IMPULSE_HHT_SLOW_STARTUP_NOT_AFSD；真实AFSD新增0。
pilot与saved_red_initial是独立成本/私有回归，不得拼接到S0时序。

- 新schema=`HHT_CURVED_CONTACT_IMPULSE_ELASTIC_BE_STATE_V1`，integration_definition=
  `HHT_DENSE_PATH_ACTUAL_CONTACT_IMPULSE_CORRECTION_V1`。不能按旧线性接触schema解释。
- `force_N`是真实端点物理力，`acceleration_m_s2`是数值HHT加速度，不等同force/m。
  `last_interval`绑定实际接触冲量/材料功、DeltaJ/DeltaW、原始Gauss节点/权重、
  壁速/滑移/法向及切向冲量、signed摩擦热、积分误差估计和动能/数值功账，全部纳入state_hash。
- `physical_force_work_J`是修正后的总力功；`baseline_linear_force_work_J`仅作基线。
  `curved_algorithmic_filter_work_J=vbar·DeltaJ-DeltaW`是signed数值功，绝不能加入温升。
  原始Q按原热分配进入BE一次；法向势/功与墙功分开记录。
- 累积面冲量/墙功使用实际积分行，不再假设每个节点权重为dt/2；当前热库面归属仅限不共享粒子的两面基准。
- `REJECTED_UNREFERENCED.npz`的外层accepted=false/restart_allowed=false为权限事实；
  嵌套candidate_state只表示数值根返回过，不能从中提取restart_allowed来绕过拒绝。
  history无stage014，6us不是有效末态。VTK只含发布态、已解字段，不补帧。
- `saved_red_initial/PRIVATE_NOT_RESTARTABLE.npz`是显式保留物理数组的数值方法迁移回归，
  未发布；其PASS只说明原失败区间通过新检查。旧失败仍是有效旧方法证据。
- `logs/process_curved_20261003/replay_initial`是已发布4→5us一次重算校验，无新接受态；
  完整属性/state/last_interval/累积账bitwise PASS，不代表跨事件重启或时间收敛。
- 接触whole_path_passive来自明确Coulomb公式和source-bound几何/零点证书；有限Gauss差值仅为积分误差估计。
  新精确Coulomb支路要求回调与公式舍入内一致，不能悄悄覆盖其他摩擦律。

独立参考仍只认证5.9375us前缀；粗时序误差及延迟停滑见PROCESS_CURVED_CAUSAL_REVIEW。
后处理不得改变冻结物理、私有/拒绝权限或给未验证的时刻补上模拟字段。

## 当前：慢载局部启动与私有诊断（2026-10-03）

`outputs/process_slow_20261003`分类LOCAL_SLOW_LOAD_FULL_STARTUP_FEASIBILITY_NOT_AFSD。
pilot有一个250ns成本步，S0有8步到2us；两者不是连续拼接，new_afsd_steps=0。
NPZ/VTK已解字段及HHT历史定义沿用下节；5Hz壁位姿/完整矢量滑移记录在boundary_snapshot，
不再用旧.05m/s或z符号判定停滑。新制造输入/几何/边界及参数单位来源在PROCESS_SLOW_CONTRACT。

- `accepted_history.json`是唯一发布权威；每个接受态先完整读回，再原子登记。
- `REJECTED.json`保存父hash/最后接受时刻/失败区间/实际诊断；`rejected_boundary_evaluations.json.gz`
  是候选边界调用，不是已接受时序。数值候选序号不是物理事件时间排序。
- S0首失败为face0/PID3，预测值方向误差；未发生停滑/负热，不能用独立参考的PID4事件替换它。
- `reference_repair1/reference_history.npz`保存两档在相同time_s的states/coarse_states；
  仅0–5.9375us共同前缀有误差界。各自最后前事件时刻/REJECTED候选另存，不能作20us末态。
- `predictor_diagnostic_initial/PRIVATE_CANDIDATE_NOT_RESTARTABLE.npz`仅为私有数值根，
  report.accepted=false/restart_allowed=false，虽原始残差收敛仍最终方向FAIL；不是检查点，
  不能添加到任何accepted_history。候选探索暂缓检查不等于最终门被放宽。
- `runner_mock_publication_*`均MOCK_NOT_PHYSICS，测试事务而非力学，不计接受物理步数。
- 参考source/config哈希、pilot/preflight/history绑定及每个输出artifact哈希必须匹配；
  初次参考元数据修复前日志保留，不能用其旧source授权新执行器。
- 所有独立读回/raw误差见PROCESS_SLOW_OUTPUT_REVIEW；它不证明20us收敛、AFSD输运或生产可行。

禁止用私有/拒绝态补齐时序，禁止把接触插值误差转成温升，冻结物理不由可视化修改。

## 当前：局部HHT真实三维基准，不是AFSD输运（2026-10-03）

`outputs/process_hht_20261003/sweep_repair1`分类
LOCAL_BENCHMARK_HHT_CURRENT_3D_TWO_FACE_THERMAL_NOT_AFSD。
8/16/32档各到2us共56接受态；new_afsd_steps=0。最终FAIL_LOCAL_TEMPORAL_ACCURACY：
仅细档速度绝对误差超界；不能把内部stage accepted解释成全过程精度PASS。

- report.json：真实命令阶段、source/config/reference/preflight/artifact SHA256、墙钟/内存、
  同终时误差与容限、new_afsd_steps=0、production_authorized=false。
- 每档accepted_history.json是唯一发布权威。候选NPZ先保存、全字段读回和source核查后才原子发布；
  未在history中的文件不允许公开重启。失败REJECTED保留最后发布时刻/评估候选，不能冒充末态。
- accepted.npz使用PySPH设施；solver_data schema=HHT_CURRENT_GEOMETRY_ELASTIC_BE_STATE_V1。
  含全部粒子数组/ID/dtype、canonical acceleration[m/s²]、真实endpoint force[N]、边界完整历史、
  parent/state/particle/parameter hashes、time_s、dt、显式壁速跳变左右极限与数值储能变动。
  事件处canonical acceleration以solver_data为准，不能把保留的PA导数scratch当右极限。
- accepted.vtk只写已解位置、速度、canonical加速度、温度、密度、压力、六应力、等效应力、
  塑性应变、质量、持久ID和被动标签；不写未认证自由能/未解e。时序由history和frames.json绑定。
- diagnostic.json.gz：根残差、候选评估/线搜索、实际空间功、whole-path摩擦/滑移、
  raw Gauss Q、BE热分项及实际保存T的余额、法向功/势差、力过滤算法功与根残差功各自独立。
  算法功可正可负，不是物理热；不以总能量余量补热。低率屈服裕度只作弹性域诊断。
- control_switch保留两侧边界/加速度/力及physical hash；物理场不清零。
  full_algorithmic_energy_jump_certified=false。全接触Coulomb方向只做注册采样误差检查，
  全路径被动性/零滑移另作多项式检查，不能混为精确连续摩擦律证明。
- integrals.json记录两独立面的累计冲量[Ns]、壁功[J]、原始摩擦热[J]、热库输入[J]，
  signed filter/root功[J]。本基准同一粒子不跨两面；不得把此面热累计直接用于共享粒子角接触。
- reference_history.npz是独立连续半离散参考，不是实测数据或AFSD生产输出。
  comparison.json按各档共同采样时刻最大误差比较，并保留1us左右力极限；不是不同失败终时对比。
- restart_hht_probe：只允许同源PASS三档后的fine stage016→032逐字段/积分bitwise重算。
  本次精度门FAIL，所以实际重算NOT_RUN；7项接口测试不能替代运行证据。

详细边界/缺项见PROCESS_HHT_REPORT、FIXTURE_CONTRACT、RUN_CONTRACT与RESTART_CONTRACT。
正塑性功、HHT共享多面/执行器/停滑事件、真实方孔/入口与AFSD输运均不在本次通过范围内。

## 当前：事件响应候选诊断与拒绝门（2026-10-03）

`outputs/process_event_20261003/response_probe`是未推进物理的数值探针；accepted_history=[]，
total_evaluations=0、accepted_steps=0，正式父态仍stage004。REJECTED不能作为重启点。
source/config/parent与REJECTED数值线性化锚点分别hash绑定，不能混用两种历史。

- particle_000001.json.gz：available=false、accepted=false、B/C=null，单位s/kg；
  每probe记录force[N]、端点/实际均速[m/s]、停滑事件、锥状态、signed Q[J]、
  原舍入界、完整kick结果hash、h/h2中央及两侧导数。不是实际外力/物理时序输出。
- 有界包装器结果另含attempts、refinement_certificates、refinement_halvings<=2；
  kick_calls/kick_seconds为全部attempt总和，probes顶层只属最后一次，不可两者混数。
  原实现和所有失败仍独立保留；smoothness证书只决定是否值得细化数值差分。
- Schur的B/C和matrix为候选数值量；不能把代数反馈或反力写入真实接触/执行器功账。
  physical trial、源项和raw接受规则未换；当前没有实际完整cloud矩阵或accepted.npz。
- preflight_stable为138 PASS/1 retained FAIL的数值限定门，不是all-tests-PASS；
  preflight_refined为144 PASS/2 FAIL，FAIL_PREFLIGHT。独立wiring1 PASS不覆盖失败门。
  independent_readback.json只证明来源/产物/正式状态一致；process_accuracy未认证。
- refined cloud、cost、fullroot及同终时时间/空间敏感性均NOT_RUN。
  PROCESS_EVENT_NEXT_WORK_PACKAGE是设计，不输出虚构HHT加速度或算法能量字段。

详细限制/命令见PROCESS_EVENT_REPORT与logs/process_event_20261003/COMMANDS.md。

## 当前：具名联合候选与原始残差控制（2026-10-03）

`process_global_20261003`分类STRUCTURAL_ROOT_BE_THERMAL_CONTINUATION_APPROXIMATE_NOT_PROCESS_VALIDATED。
对应trial为PROCESS_GLOBAL_EVENT_AVERAGE_BE_THERMAL_APPROXIMATE；全局SI，未缩放物性。
preflight103项PASS不增加物理时间；一次500ns/root_attempt最终BLOCKED_NO_NEW_ACCEPTED_STATE，
accepted_history=[]，accepted_steps=0；REJECTED_candidate不可恢复成有效物理末态。
source/config/父checkpoint/hash、成本探针和32个candidate、最新完整trial均在目录中。

- 五场raw残差与maximum_scaled_residual是接受依据。maximum_scaled_proposal仅数值候选量。
  numerical_proposal_only包含完整原缺陷、Schur rhs[N s]、density回代[kg/m³]、
  代数反力[N s]与线性backward_error；不能冒充实际执行器反力或接触冲量。
- BE热字段使用pair_heat_into_i_J、pair_conductance_W_K、particle_energy_J、ledger；
  former_explicit_dt_bound_s仅旧显式界诊断，不伪造旧pair_power_W字段。
  boundary热记录统一T1，旧Tmid另存diagnostic；逐事件Q保持原signed值和舍入证书。
- 新具名checkpoint/VTK接口已通过制造零平衡存储测试，但本次没有生成真实accepted.npz。
  restore是内部全字段读回，不授权孤立NPZ公开恢复；重算重启等价NOT_RUN。
- process_accuracy=NOT_VALIDATED；旧法向/总能量诊断和Jaumann风险不因求根残差下降消失。
  未解e/scratch不当真实能量输出；无新材料输运/沉积结果。

细节/命令见PROCESS_GLOBAL_REPORT及logs/process_global_20261003/COMMANDS.md。

## 当前：过程组件与离线重放，不是新增物理时序（2026-10-03）

`scripts/check_process_readiness.py --output 新目录`产生逐组件.log和report.json：
classification=PROCESS_COMPONENT_PREFLIGHT_NOT_COUPLED_ADVANCE，project_complete=false，
accepted_progress=false、new_coupled_steps=0、production_authorized=false。
PASS_COMPONENT_PREFLIGHT_PROCESS_STILL_BLOCKED只说明独立组件测试/旧证据完整性通过；
source_sha256、previous_freeze_sha256、actual command/exit/墙钟/测试数均记录；remaining_gates不得隐去。

- finite_tool_geometry.query：signed_distance_m、closest_points_m、boundary_features、normal或NONSMOOTH、
  application_point_m和两位置壁速；normal_penalty为单实体势能[J]/法向力[N]。NOT_AFSD，非旧逐面罚力等价。
  不输出未解摩擦/热/粒子运动。非光滑活跃拒绝不是成功通过孔口。
- inlet_protocol：候选PA+registry+event，accepted=false/NOT_INTEGRATED；new IDs单调分配，完整字段unit/source/status，
  UNRESOLVED scratch不能当已解。source_increment质量[kg]/动量[kg m/s]/显热[J]/动能[J]/体积[m3]。
  显热可负，累计cp/Tref固定；pressure_flow_work/elastic_energy/mechanical_reaction_work NOT_COMPUTED。
  registry全属性按ID绑定+hash，只证明存储身份；实际推进后的授权rebind尚缺，不能由孤立候选恢复为物理状态。
- implicit_thermal：COMPONENT_SOLVED_NOT_INTEGRATED；temperature、pair_conductance_W_K、pair_heat_into_i_J、
  particle_energy_J、ledger、linear_system_residual_J；former_explicit_dt_bound_s是比较诊断，不是BE验收上限。
  摩擦raw和舍入证书由已验证力学调用者给出；塑性/壁/摩擦分账，不补残差热。
- outputs/process_completion_20261003/thermal_replay.json：PASS_OFFLINE_IMPLICIT_HEAT_REPLAY_NOT_COUPLED_ADVANCE。
  分别读取4个旧stage的各自初态、原dt/真实源和中点邻域，仅重新解热；新T不驱动下个旧stage。
  温度数组仅离线热解，restart_allowed=false；旧显式温度差是诊断，不是新耦合接受态或论文误差。
- structural_mobility：FROZEN_CONFIGURATION_ELASTIC_PRECONDITIONER_NOT_PROCESS_SOLVER，K[N/m]、A[kg]、
  delta_velocity_m_s、reaction_impulse_Ns、端点/中点contact_response[1/kg]。
  仅未预应力线弹性冻结构形；full_tangent=false；零刚模/条件数原样报告，不做pinv或偷偷约束刚体轴。

各组件实际结果和未集成限制见PROCESS_COMPLETION_REPORT及对应COMPONENT/INDEPENDENT_REVIEW文档。
当前真实三维接受轨迹仍上一节4个50ns阶段；本节任何PASS都不增加物理时间。

## 当前：有界整体求根与接受态连续推进（2026-10-03）

分类UPDATED_STAGE_ANDERSON_CONTINUATION_APPROXIMATE；全局SI，production_authorized=false。
目录outputs/nonlinear_continuation_20261003；不是论文验证或旧显式求解器的同源恢复。

- cost_probe/manifest：COMPLETE_UNCOMMITTED_COST_PROBE，未增加物理时间；完整耗时和source/config用于选总evaluate预算。
- continuation/manifest：COMPLETE_REGISTERED_FOUR_STAGES_NOT_DEPOSITION；342点，61evaluate，4接受步，终时850.393633ns。
  source_time_s=650.393633ns，parent为旧step62；settings总上限72，budget实际选69，非69个物理步。
- accepted_history.json为正式发布列表；孤立NPZ不是已发布接受态。每阶段完整验收/输出/严格重载后才原子发布。
- stage_NNN/accepted.npz保存全PA、原始标签/坐标、S/ep/T/rho/v/x、真实累积塑性功/return[J每粒子]、source/config/全属性hash。
  state.valid_properties列实解/持久属性；e/ae/ax/ay/az为scratch_not_solved，不冒充内能。
  restart_allowed=true仅对同变体严格内部loader；外部安全恢复尚无CLI，需正式history授权；两轨迹重求解等价NOT_RUN。
- accepted.vtk真实位置/速度/温度/密度/压力/六应力/等效应力/ep、ID与标签；不输出伪造自由能。
- tool_pose.json实时间对应工具位姿；acceptance.json保留原容限、候选-map功差/Es诊断/法向时间残差/角动量项。
- accepted_trial.json.gz是本接受candidate实际evaluate；candidate_NNN.json.gz和iterations.json保存非线性候选/分支诊断。
  非线性候选不是物理时序；未接受candidate与失败时REJECTED不能取代正式末态。
- transport.json按ID对齐相邻接受态；固定lab z=layer=1.5mm，上/下穿越中心事件质量与唯一ID质量分开。
  累计事件质量可含重复穿越，不是净沉积质量或Fig.8二维面积比。全棒诊断与局部标签分开，不影响力学。
  参考系只随工具平移；原始肩下75点不是新穿越。NO_CENTER_CROSSING_NOT_DEPOSITION本轮各步均成立。
- checkpoint ledger_scope=CONTINUATION_ONLY，PA累计塑性/return另包含parent史；禁止二者直接比较同起点。
  ledger.boundary只累加真实功/冲量/事件次数；boundary_state_latest保存当前中点法向势能、接触数与当前极值，
  全史最大穿透/最小raw Q另取极值，不把势能或最大值逐步相加。热ledger均J增量，法向端点账在acceptance另存。
- independent_readback.json：PASS_ACCEPTED_CONTINUATION_READBACK_NOT_DEPOSITION，43源/97产物/23411子区间；
  不调用被测物理生成参考，确认4态/真实账/位姿/VTK/穿越。未验证自由能、同终时收敛或工艺尺度。
- process_feasibility_corrected.json才是有效成本/几何解释；process_feasibility.json及首日志有负成本错误，INVALID保留。
  名义距离/执行器速度不是材料旅行时间；固定步线性成本外推不授权长时运行。详见CORRECTION。

下面旧段中的“条件接受NOT_RUN”仅指旧Picard当轮；不覆盖本段新四步结果，也不解除旧自由能风险。

## 同史更新构形trial / 唯一全局步（2026-10-03）

分类`UPDATED_CONFIGURATION_EVENT_AVERAGE_STAGE_IMPLICIT_APPROXIMATE`，全局SI；目录outputs/updated_trial_20261003。
不是原HostLeapfrog普通重启。step62来源保留原非零v/S/ep/T，不清零历史；新配置hash包含命名算法与dt/预算。
material_scope真实运行必须UNCHANGED_FKS_APPROXIMATE_PARAMETERS；制造law仅独立单测，不可晋级真实接受态。

- `cost_probe/manifest.json`：COMPLETE_COST_PROBE_UNCOMMITTED、accepted_progress=false，完整含I/O耗时3.81249s。
  `uncommitted_trial.json.gz`保存candidate+evaluation；不是已推进一步的物理末态。
- `global_attempt/manifest.json`：本次BLOCKED_REJECTED_GLOBAL_STAGE、32evaluate、accepted_progress=false。
  `iterations.json`存每次五场map残差（rho按初态rho逐点归一）、最大缩放残差和evaluate墙钟。
  同一个物理dt的32次迭代不是32个物理时间步；温升、ep均不能累加32次。
- `latest_uncommitted_trial.json.gz`：最后一个完整evaluate及其实际candidate；阶段accepted=false/restart_allowed=false。
- `REJECTED_candidate.json.gz`：明确classification=REJECTED、accepted=false、restart_allowed=false；error、candidate、
  last_completed_evaluation、iterations和failure_context。早期异常时last_completed可能属上一候选，
  必须读`last_completed_evaluation_may_refer_to_previous_candidate`，不可假定同一配对。此次迭代上限失败配对是同次最后计算。
- `last_accepted_reference.json`引用未改变的旧plastic_accepted.npz的项目相对路径/hash/time/step；
  不是另造新的接受checkpoint。本次没有accepted_snapshot.npz或accepted_stage.vtk。
- `resolved_config.json`保存原SI参数及新数值settings；逐项物理来源仍由已冻结旧resolved_config/合同索引提供，
  不把本轮预算/迭代参数假称论文参数。
- source_sha256、input_sha256、parent_configuration_sha256、configuration_sha256和artifact_sha256区分新算法/旧物理/文件；
  禁止将重新命名或转换格式视为同源验证。

candidate五字段：position_m[N,3]、velocity_m_s[N,3]、mean_velocity_m_s[N,3]、density_kg_m3[N]、temperature_K[N]。
mean_velocity是独立时间区间平均未知量；map使用真实接触子区间Σh(vbefore+vafter)/2dt，不是端点均速。
evaluation存map/residual、连续方程残差[kg/m³/s]、实际particle_id、同一中点x/rho/T/vbar/L/sigma、i/j/g/B和质量。
material存S1[Pa]、ep、区间塑性功[J/粒子]、return耗散[J/粒子]和等效偏应变率[1/s]，都从同一原史起算。
boundary保留全部particle/face/subinterval/冲量/反力/原始Q与舍入证书、法向势能和conductance；
thermal保留逐粒子热源、Tmap、成对导热、真实热行和dt界。mechanical return/残差不是热。
before/mid/after有限feature守卫结果、实际空间条件数、功/净力/净力矩诊断均独立字段，不以渲染证明三维。

`independent_readback.json`为最后完整候选的只读诊断：本次status=PASS_BLOCKED_EVIDENCE_ONLY_NOT_ACCEPTED，
accepted_physical_stage=false。dominant_residuals按五字段保存PID/轴/构件/实际接触feature、候选值/map值与残差；
还有独立Q/冲量/空间/热账、候选与map功差和角动量项。该PASS只说明检查内容一致，不能覆盖global_attempt的REJECTED。
没有读回全部32次完整候选档案，也没有授予自由能认证或新物理步。

条件接受接口（本次真实路径NOT_RUN）：只有五场残差、空间、热账、动量、约束、特征域及来源复核完成才保存。
保存**最后实际evaluate的candidate**，不得悄悄snap到未经evaluate的map。接受账记录候选减map动能/内力功差、
总能量诊断残差和法向时间残差，不宣称Es为自由能；不加残差热。
端点角动量额外项ΣmΔx×(v0+v1)/2单列；冻结中点外矩与真实轨迹外矩差异仍APPROXIMATE，不补外矩。
若未来输出accepted_snapshot，元数据restart_allowed=false、valid_properties明确实解字段；继承scratch/e/旧累计功不是新解。
真实阶段续算校验NOT_RUN，不能通过旧restore绕过。本轮VTK/checkpoint不输出拒绝态为有效末态。

## 过程时间残差与只读路线评估（2026-10-02）

outputs/process_route_20261002/assessment_complete.json 为完整只读结果；没有新增真实轨迹/VTK/checkpoint。
分类UNSCALED_INERTIAL_RESIDUAL_INTERFACE_ONLY，status=READ_ONLY_ROUTE_ASSESSMENT_COMPLETE_NOT_NEW_PHYSICAL_SOLUTION，
accepted_progress=false、production_authorized=false。前assessment.json/assessment_final.json是保留的入口FAIL，不能当有效结果。

纯process_residual.midpoint_residual输入m[kg]、x[m]、v[m/s]、dt[s]、各力[N]、free_axes和规定端速度。
输出均N×3：kinematic_residual_m、dynamic_residual_N、quasistatic_free_residual_N、
constraint_velocity_residual_m_s、free_reaction_force_N、inertial_force_N、physical_force_N。
它不算接触力、不解全局根、不更新物理时间、不求热、不读取标签；要求外部从同一候选提供真实力。

assessment每个all_steps记录全部原67步的dt/time、Jn/Jt/R逐face/子区间合计（内部核验）、
inertial_free_force_norm_N、quasistatic_free_residual_norm_N、dynamic_residual_norm_N、
internal/normal/tangent_free_force_norm_N及force_denominator_N。
quasistatic_force_ratio=||Rq||/(||P Fint||+||P Fn||+||P Ft||)，不是全局合力比或实验误差；分母0输出null+UNDEFINED_NO_FORCE。
没有设置该比例的生产通过阈值。没有把小名义Mach当作瞬态已静平衡证据。
summary列步骤覆盖、全域比例/残差极值与非静平衡步清单；first/first_plastic/last仅提取实际步，不构造新参考轨迹。

actual_contact_work_J来自原完整子区间已审计账；endpoint_contact_work_surrogate_J只为总Jc·端点平均v。
endpoint_minus_actual_contact_work_J只诊断，surrogate_status=NOT_FRICTION_HEAT_NOT_PHYSICAL_SUBINTERVAL_WORK，
不得将任何差加热或据全步动量残差替代滑移事件处理。
nominal_scales保存原物性推导声速/工具速度/Mach；不更改质量/时间/模量。

input_sha256锚定原freeze、已冻结run manifest、全部33来源及产物；上游5源通过已冻结manifest传递核对实际hash，
项目28源另外匹配直接freeze索引。source_sha256为新reader/残差接口/事前合同，不伪称原物理source变更。
decision只记录当前checkpoint不能直接去惯性；新全局非线性步/过程时间AFSD均NOT_RUN。
独立弹簧/振子测试日志标制造解析基准，不能与原342点SPH输出混用。

## 四关统一热力联动适配器（2026-10-02）

分类ENGINEERING_APPROXIMATE_COUPLED_AFSD_LOCAL_DEBUG；目录outputs/coupled_afsd_20261002。
全局SI；原始小组件/历史真实AFSD输出保留，不将不同几何数据拼接。
入口run_coupled_afsd及只改停止预算的run_coupled_short_window，物理均IntegratedSolver。
schema沿用已有字段；这里T、应力、ep均实际求解，不伪造完整自由能e。
关卡/合同/近似见COUPLED_AFSD_REPORT、CONTRACT、SHORT_WINDOW、CODE_REVIEW。

- manifest.json：classification、production_authorized=false、mode、status、requested/accepted时间/步数、
  accepted_dt_s、step_wall_seconds、RSS/PA字节、source/config/artifact SHA256、final_record和错误。
  原成本mode=COST_PROBE；短窗mode=REGISTERED_INTEGRATED_TRIAL，requested_end_s=7e-7。
  COMPLETE_REGISTERED_INTERVAL_NOT_SUSTAINED_DEPOSITION只指该终时，不是原2us或宏观沉积。
- resolved_config.json：全部SI参数/来源/状态；短窗另存original_explicit_configuration及
  short_window_budget_override，唯一变化为观察终时；保留原2/4us NOT_RUN预算。
- history.json/jsonl：每个接受态（含初态）的particle_id、component_id、material_id、原始tracer/wholebar/ASRS、
  actuator观察标签、position_m[N,3]、velocity_m_s[N,3]、temperature_K[N]、density_kg_m3[N]、
  pressure_Pa[N]、deviatoric_stress_Pa/stress_Pa[N,3,3]、equivalent_stress_Pa、ep、strain_rate_per_s、
  particle_mass_kg、逐粒子累计plastic_work_J/return_dissipation_J、tool_origin_m/tool_angle_rad。
  正压为压缩；stress=S-pI。工具位姿只驱动真实边界，不写自由材料路径。
- steps.jsonl.gz：每个接受步完整journal；含x_before/mid/after、v_before/after/work_velocity、
  S/ep/T/rho前后、mass、邻域i/j/raw径向g/B、L_before/L_midpoint、force、真实contact traces、
  free_axes/原始z0控制、normal_before/midpoint/after及feature transitions、thermal每粒子源与累计账增量。
  子区间记录真正Jn/Jt/约束冲量、相对中点速度/功，不用全步总速度代替摩擦子区间功。
- series.json和accepted_*.vtk：实际接受时序，每10步及首末输出；本窗8帧，非逐步VTK。
  VTK含真实3D位置/速度/温度/密度/压力/完整应力/等效应力/ep、ID/初始位置/观察标签。
- last_accepted.npz：最后接受态；nonzero_accepted.npz为step2；plastic_accepted.npz为第一次新增塑性（本窗62）。
  必须accepted=true/restart_allowed=true及源/config/参数hash匹配；保存初始K/P/L/U/能量和完整账。
  真正重启执行另存short_restart.json，不将仅检查checkpoint数组称为执行过重启。
- REJECTED_candidate.json仅在失败生成，含来源、候选数组/journal、原异常contact/feature上下文和last_accepted。
  不能作为有效物理末态或restart；整步材料/热失败回滚全部数组、账、time/count。

控制范围：全部板棒可变形，只有原始顶层0.5mm棒料的z规定=-feed；xy自由。
标签actuator/tracer/ASRS均不可选受力；真正约束由原始z0与显式程序决定。
基板没有固定底层。有限棒/封闭质量系统：external_mass_kg=external_feed_enthalpy_J=0是明确无连续入口，非默认缺值。
feed_plane_crossed_mass_kg是既有材料几何穿越量，不能冒充外部质量注入。
接触ID由持久particle_id:feature构成；0底座/1肩面/2..5旋转孔壁；本窗1未活动。
法向势能按实际有限feature集合记录，不能继续用消失旧面的切向或U。

热学/功量定义（J，正号为入材料，除显式out）：

```text
P = V*qnew*delta_ep                     # 实际材料塑性功
Qfr_raw = sum_subinterval(-Jt dot u_mid) # 原符号，不裁负
Qplastic = chi*P; Qfr_material = partition*Qfr_raw
Qcontact = dt*hA*(Twall-T); Qinternal = dt*成对导热功率
deltaEth = sum(m*cp*(Tnew-Told))
Qwall_out = (1-partition)*Qfr_raw - Qcontact
Rn = Uafter-Ubefore+Wnormal_relative
Ediagnostic = K+Es+U+Eth+(1-chi)*P
R = Ediagnostic-E0-Wwall-Wact+Qwall_out
```

boundary.contact_power_W只含接触导热；摩擦仅thermal分配一次。return_dissipation、Rn、R都不作热源。
cumulative保存raw摩擦、两侧摩擦热、塑性/非热化塑性、接触热、内部导热净和、热账残差、
工具/底座/执行器冲量与角冲量、normal_material/wall/relative_work、normal_temporal_residual。
material_old_work_J与midpoint_stress_work_J分别是材料kick前梯度与真实区间速度梯度的功；不得合并或删差。
Es/total_J仅诊断；e为上游scratch，不是认证内能；material_energy_certified=false。
energy_status的WITHIN_SCREEN不能当作功势一致或物理精度通过。本窗全局R不是热残差，二者量级差很大。
thermal_dt_bound_s是真实内部+壁行和界；超界拒绝回滚，不是自动热子步求解。

独立short_readback_initial.json统计全部68态/67步/189382子区间，PASS_SHORT_WINDOW_NOT_2US_NOT_DEPOSITION。
mechanisms字段分ACTIVE、ACTIVE_FINITE_WINDOW和NOT_DEMONSTRATED；只有2步塑性不叫持续流动。
三工况/收敛/定量论文比对NOT_RUN；读取器PASS只证独立身份式/已存数组，不是另一个连续体解。

## 局部轴向驱动与释放组件（2026-10-02）

分类ENGINEERING_APPROXIMATE_LOCAL_AXIAL_FEED_NOT_AFSD；outputs/axial_feed_20261002。
合同/范围见AXIAL_FEED_CONTRACT/REPORT。全局SI；不覆盖旧固定底层/旋转肩面变体。
同一PySPH粒子/NNPS/积分/输出，history、steps、VTK、accepted/rejected语义保持。
单点reference单独标MANUFACTURED_AXIAL_ACTUATOR_REFERENCE_NOT_AFSD，不冒充SPH固体。

本组件全部64点可变形，没有材料底部夹持。boundary_role=2、pa.actuator仅是原始顶层观察元数据；
不能用其决定受力。真正控制从原始z0与冻结几何/程序求得，不使用tracer/ASRS/当前位置。
control包含eligible、active、free_axes[N,3]、prescribed_velocity[N,3]、phase及commanded_displacement_m。
eligible为原始顶0.5mm材料层；active只在t<400ns；free_axes的1是自由、0是约束。
受控只z=-0.00212m/s，xy自由；RELEASED所有轴自由，prescribed_velocity=0不表示把实际速度置零。
commanded_displacement_m只记录边界命令，不是自由材料轨迹或沉积长度。
control_before/midpoint/after明确时间层，释放右侧状态t=400ns即RELEASED。

初始量显式输出initial_kinetic_J、initial_linear_momentum_kg_m_s、
initial_angular_momentum_kg_m2_s、initial_normal_potential_J；均按persistent ID规范求和。
当前linear/angular_momentum及其residual扣除初值和累计外冲量/矩，不假设初始静止。
history提供cumulative_base/actuator_impulse_Ns及angular_impulse_Nms向量；cumulative账还保存
base/actuator_impulse_{x,y,z}_Ns、base/actuator_angular_impulse_{x,y,z}_Nms标量，便于checkpoint。
traces中actuator_force_on_material_N=constraint_impulse/dt，正号为对材料作用；对执行器反力为负。
actual contact_rows均BASE；gap基于真实底座z=-0.01m，法向+z，mu=0隔离本轮法向账。
normal_before/midpoint/after对全64点保存gap、signed penetration、active、Fn和弹簧U。
单边penetration的正部是接触几何定义，不是负功/热裁剪。

功定义：Wnormal=Jn·v_material_mid，Wnormal_wall=Jn·v_wall=0，
Wnormal_rel=Wnormal-Wnormal_wall；Rn=deltaU+Wnormal_rel保留原符号，不补热。
kinetic_contact_residual=K-K0+Q-Wwall-Wact-Winternal-Wnormal_rel。
energy_diagnostic_residual=K-K0+Es+U-U0+Pplastic+Q-Wwall-Wact；Es未获自由能认证。
energy_status区分ENERGY_WARNING与DIAGNOSTIC_WITHIN_SCREEN_NOT_CERTIFIED_FREE_ENERGY。
boundary_mass_flux_kg_s=0、system_mass_status=CLOSED_NO_INFINITE_INFLOW；不得误读为连续供料。
plastic/return耗散独立存，不能用法向或时间残差替换。未解T/e不在VTK伪造。

VTK增加actuator_eligible/actuator_active并保留真实x/v/rho/p/S/等效应力/ep、ID/原始位置/标签与法向字段。
actuator/boundary_role即使被修改也不改变力学；后处理不得把它们解释成额外固定点。
last_accepted.npz与release_accepted.npz均accepted/restart_allowed=true，后者400ns。
恢复校验classification、source/config hash、初始K/P/L/U；REJECTED不可restart。
*_restart.json记录实际释放checkpoint续两步与保存轨迹差，不是新的独立物理解。

manifest mode区分TWO_STEP_COST_PROBE与FULL_REGISTERED_AXIAL_COMPONENT；BLOCKED不能称有效末态。
comparison.json保存482帧独立读回、161共同帧最大绝对分量差及逐ID接触活动变化区间。
PASS_AXIAL_FEED_READBACK_NOT_AFSD只证保存账/算术/VTK一致；SENSITIVE_WARNING是整体dt比较状态。
contact事件区间为接受步采样，不是独立精确波到达时间。reader没有独立解整个非线性材料轨迹。
production_authorized始终false；source/artifact哈希和freeze只保证证据完整性。

## 法向开离组件与解析基准（2026-10-02）

新接口不覆盖旧固定z变体。详NORMAL_CONTACT_REPORT/CONTRACT；全局SI，显示单位另转。
单点分类MANUFACTURED_NORMAL_IMPACT_NOT_AFSD；3D分类ENGINEERING_APPROXIMATE_NORMAL_RELEASE_SPH_NOT_AFSD。
目录outputs/normal_contact_20261002；manifest源/配置/全部artifact SHA256，production_authorized=false。
mode与status区分2步成本探针、完整注册600ns3D、BLOCKED；单点有独立PROBE_COMPLETE状态。
最初未加重启guard的单点输出及源码副本保留；正式粗静壁使用static_dt005_guarded。

3D保留原history/steps/VTK/checkpoint结构，但boundary_role=2现在仅标识初始顶面，
不是z夹具；只有role1底层xyz固定，role0/2全部xyz自由。标签不定义接触或施力。
normal_contact按初始顶层persistent particle_id记录gap_m、signed_penetration_m、
penetration_m、active、normal_force_N、potential_J以及tool位姿/速度/offset。
signed_penetration=dx/2-gap，正部仅定义单边几何接触，不是热/耗散裁剪。
normal_before/midpoint/after分别是半漂移前/接触几何/接受末态，不能互换时间层。
实际接触traces逐feature记录normal_force_N、Jn、Jt、约束反力、各子区间原始Q。
壁速z包含刚工具抬升，底座静止。相同点位置用于速度与力功共轭。

法向账：normal_work_J=Jn·v_material_mid；normal_wall_work_J=Jn·v_wall；
normal_relative_work_J=两者之差；normal_temporal_residual_J=deltaU+normal_relative_work。
U0显式初始接触储能；global diagnostic=K+Es+Pplastic+Q+U-U0-Wwall-Wconstraint。
Rn不是热，Es不是获认证自由能。负壁功可代表执行器取回初始储能。
raw负Q不裁；reader negative_raw_intervals是独立重算计数，不等于保存trace的负数计数，
接近舍入尺度的符号可能因求和顺序不同。本轮细档保存4项，重算1项，日志单独列出。

3D VTK新增signed_gap_m、normal_signed_penetration_m、normal_penetration_m、
normal_force_N、normal_potential_J、contact_active，及旧已求解3D力学字段。
这些平面法向诊断只对冻结远离孔/外缘的块成立，不是任意有限工具/多面接触输出定义。
temperature是规定输入，不输出假温度解；checkpoint的e仍scratch，不是热能。
单点VTK只输出ID、tracer、AS/RS、mass及x/v；不伪造stress/rho/T。

last_accepted.npz与prelift_accepted.npz均accepted/restart_allowed=true，后者在200ns，
仅允许同source/config组件恢复。REJECTED另存，不可作为物理末态或重启。
两个实际checkpoint续2步的readback报告与主运行分目录存，不延长主物理终时。
opening_brackets_by_id_s为接受步前后活动性变化的时间区间，不是精确事件根。
block_comparison.json含121共同帧六字段差与开离区间；SENSITIVE_WARNING不是PASS_ALL。
所有新状态仅组件证据，不授权AFSD生产；旧schema和失败日志保留。

## 旋转肩面—可变形三维SPH小块（2026-10-02）

独立companion为ENGINEERING_APPROXIMATE_ROTATING_DEFORMABLE_BLOCK_NOT_AFSD，
不是上一节平面四质点，也不是三工况沉积。接口/参数合同ROTATING_DEFORMABLE_CONTRACT。
outputs/rotating_deformable_20261002/分cost_probe、dt5e9、dt2p5e9新目录；
manifest的mode区分两步探针和完整200ns组件。生产授权始终false。

history每帧存particle_id、boundary_role(0自由/1底层xyz夹具/2顶层z夹具)、
tracer_type/as_rs_tag、x/v、真实连续密度、压力、3x3偏应力/全应力、ep与速率。
rate_field_status初帧NOT_EVALUATED、以后PRECONTACT_TIME_LAYER；不是最终速度梯度。
工具origin/angle、质量、K、Es_diagnostic、累计rawQ/Wwall/Wconstraint/Winternal、
material_old_work/midpoint_stress_work/plastic/return_dissipation都单独保存。
Es含旧剪切平方诊断与EOS体积能，不是获有限应变认证的自由能；全局残差不补成热。
293.15K是规定输入，thermal NOT_SOLVED，不在VTK输出假温度或假总内能。

steps.jsonl每步含半漂移x、before/after/work-average速度、持久ID与角色、
rho_before/rho_mid/rho_after、drho、实际i/j/g/B、condition、L_before/L_midpoint、
sigma_mid、S_before/after、internal_force、塑性及算法耗散、全部逐粒子kick子区间。
work_velocity为各子区间dt加权中点，不预设等于整个kick端点中点。
traces的肩面feature/contact_id稳定，Fn规定值和normal impulse、constraint reaction
分存；无接触的内部点也记录其内力冲量。法向载荷被顶层z夹具平衡，不是块体预压解。
自由层可有真实z响应；不能把初始四层点云本身当三维物理证据。

VTK导出已求x/v/rho/p/六个对称应力分量/等效应力/ep/速率、mass/h、ID和标签。
checkpoint使用原PySPH原子NPZ，仅固定组件restart_allowed=true；source/config哈希
校验且拒绝REJECTED。checkpoint内e仍为未解scratch，不是可用于热重启的物理能量。
失败候选另存rejected_candidate.json，包括last accepted和candidate；非有限诊断值
用显式nonfinite标记保存，不剪裁为0。已接受history/checkpoint不含失败末态。

独立reader不导入求解器，校验来源/全输出哈希、字面核导数和完整三维邻域、
伴随功/净力/净矩、子区间Coulomb及signedQ、实际dP/dL、EOS及材料时间层、VTK。
读回PASS只代表该组件的这些检查，不是连续体精度/本构能量/论文验证。
两档比较要求同参数/几何/源代码与共同200ns，报告绝对差和事前工程筛查，不宣称收敛阶。

## 旋转刚壁接口与肩面点质量组件（2026-10-02）

rigid_contact_kinematics.contact_rows为SI数组：particle_index/particle_id、
feature_id=0..5（BASE,SHOULDER,X-,X+,Y-,Y+），contact_id="ID:feature"；
normal/gap/penetration、application_point=x、contact_point=c=x-gap*n，
wall_velocity在x，wall_velocity_at_face在c，body_origin/translation_velocity/omega。
c只为有限面投影，不是完整刚体最近点；角部同粒子可有两个面，标签不作为输入。
relative_kinematics未给L时face material velocity为None/NOT_SOLVED，不补零；
给L后的重构也仅诊断，尚无转置力分配。contact_power保留raw signed摩擦功率W，
normal_relative_power_W独立；偏心torque_Nm/power_W不转成热。

outputs/rotating_contact_20261002/dt2e5和dt1e5是平面自由点质量接触动力组件，
不是3D SPH固体/AFSD。history.json含每步t、x/v、持久ID/标签、工具pose、K与累计
signedQ/壁功/正常/执行器功/残差。contact_intervals.jsonl每步四个真实kick及其
子区间、Jn/Jt/Jconstraint、before/after/midpoint速度、roundoff bound与几何row。
温度、应力、塑性、密度、NNPS内力NOT_SOLVED，不导为真实场。
VTK只含位置、速度、ID/被动标签/规定mass。必要PySPH rho/e是checkpoint scratch，
h为接触检测元数据；不要在可视化中改标物理求解字段。

last_accepted.npz用原原子checkpoint，restart_allowed=true只用于本制造组件，
production_authorized=false；同时校验config/sourcehash。失败candidate另存REJECTED，
先回滚accepted粒子；不能从REJECTED继续。manifest含来源/规范化config/全文件hash。
reference_dt1e7/5e8.json独立RK4保存literal参数、来源和参考误差；仅持续滑动有效。
readback.json的PASS_COMPONENT_READBACK_NOT_AFSD同时核验原signed功、Coulomb法则、
转矩冲量/角动量、Wwall=omega*deltaLz、同终时参考差、全部VTK/ID，非生产证书。

## PredictorOnset同终时4/8敏感性输出（2026-10-02）

outputs/process_predictor_sensitivity_20261002/n3_i8复用冻结schema，
56accepted/57帧到.2s；旧n3_i4与旧i8/NOT_RUN不覆盖。
comparison.json是后处理而非物理态：status=COMPARISON_COMPLETE，
screening_status=WITHIN_SCREEN或SENSITIVE_WARNING，不是生产/能量证书。
要求两档完整且实际reader PASS；不对partial比较不同末态，禁止插值。
common_time_s共17时刻，matched_indices标对应accepted REQUESTED帧；
metrics有名义尺度/实际分母、绝对与归一化全路径/末态/逐时刻差、原有符号末态。
力采用contact_force_N合力，非近零support+contact系统reaction_N。
events按phase:particle_id:feature_id:occurrence匹配，存分组/顺序/缺失与时间差；
只指CONE_LIMIT新到边界候选，不是全部黏滑切换。
cases有成本/步数/原材料告警/独立signed功账，不改负数或补热。
provenance含事前合同/比较脚本、23项来源、输入report和readback哈希。
5%为事前工程提示，不是paper误差、收敛率或求解门槛；两档材料能量未认证。

## PredictorOnset集成输出（2026-10-02）

新命名变体Q7_PREDICTOR_ONSET_V1；原Q7的history/journal/checkpoint/VTK物理字段
不变，原Onset的候选/事件/接受时间链也不变。新增predictor_provenance.json包含
新wrapper/config/helper共3来源、onset_source_sha256共4来源、baseline_source_sha256
共16来源，并必须进入report.artifact_sha256。不能只看旧provenance的16项。
resolved_config.predictor_variant记录新入口；180s是外部运行墙钟限制，不改变物理dt。
旧父配置60s保留为来源历史，不是当前新变体的外部墙钟上限。

attempt的iterations分为INTERIOR_PREDICTOR（corrections_used、interior_force_relative、
normal_force_N）与CONE_ROOT（iteration、root_relative、normal_force_N）。前者可保存
原始负N，不是求解出的接触态；两阶段修正合计<=14。记录不能重置计数，
预测行不能冒充收敛根。失败/预测候选都不得变成accepted。
新inspect_process_predictor_onset.py先验证三组来源/精确配置/共享预算与事件时间链，
再实际调用原Q7物理数据审计。PASS_PARTIAL_READBACK退出2仍代表整段BLOCKED；
PASS_COMPLETE_READBACK仅是此制造基准的数据完整/局部物理检查，不代表AFSD认证。

## 可行初猜诊断（2026-10-02；不产生新的accepted）

outputs/process_feasible_predictor_20261002/target_0125仅存纯分支对照。
inputs.json记录旧输入/config/source与新helper/诊断/prereg哈希；
old_comparison/predictor_comparison.json及diagnostic.json存完整状态前后指纹、
原始trial/修正阶段、共享预算、原接触评估。branch_status=ROOT_RETURNED
不是接受；本轮status=ROOT_RETURNED_EVENT_SPLIT_REQUIRED_NOT_ACCEPTED、exit2。
REJECTED_<old|predictor>_raw_NNNN.npz保存每次真实trial及ID、signed N/tau、gap，
均classification=REJECTED_DIAGNOSTIC_NOT_RESTART、restart_allowed=false。
最终分支根也标REJECTED；没有新的checkpoint/history/VTK有效末态。
raw机械残差不是cone残差；候选功与预测能量差不能加入accepted或热源。
所有raw的base_step=16、target=.125；accepted仍.1125。未求解热/惯性字段不补零伪造。

## Onset v1补充（2026-10-02）

字段仍为Q7实际accepted物理数据。新增onset_provenance.json，补充4项来源
（onset入口/配置＋恢复入口/配置），关联父16来源，进入report.artifact_sha256。
resolved_config.onset_variant显式标Q7_BOUNDARY_ELIGIBILITY_V1；原恢复14预算
及所有物性/验收限制不变。不得只用父16来源作为完整入口清单。
新reader inspect_process_plastic_friction_onset.py先验来源/精确配置和attempt
时间链，再执行原Q7数据/功账审查。onset_audit只表示对应检查，非连续体认证。

每个请求单独attempt_<base_step>_<target>.json，class=ATTEMPT_TRACE_NOT_RESTART，
restart_allowed=false；stage为BOUNDARY_PROBE/CLASSIFICATION/TARGET/
INTERIOR_BISECTION/BRACKET/ACCEPTED。ROOT_RETURNED_NOT_ACCEPTED仍是诊断候选；
只有ACCEPTED与history/journal/report一致的步才是有效物理数据。
新失败可能在state构建前抛出，原REJECTED_manifest会标candidate_available=false，
不补假数组。诊断另存的REJECTED_DIAGNOSTIC_NOT_RESTART原raw trial不能进history。
negative_normal_diagnostic同时保存迭代调用位置、signed N、约束gap和机械残差，
机械残差不是未定义域内的cone root。accepted末态仍.1125s而不是目标.125s。

## Q7R 恢复变体补充（2026-10-02）

基础Q7字段/accepted与REJECTED分离规则不变。
新增recovery_provenance.json，存实际恢复入口与配置的source_sha256、父实现
baseline_source_sha256、variant=Q7R_BOUNDED_NEWTON14、production_authorized=false。
它必须进入report.artifact_sha256。旧provenance/report.source_sha256仅为父实现16项，
不能把它单独当完整运行来源清单。resolved_config.recovery_variant记录实际覆盖，
plastic_parent.numerical.maximum_iterations=14，其他物理/容限不变。

独立检查必须用inspect_process_plastic_friction_recovery.py；返回recovery_audit
补充变体、预算、来源和配置差异验证。PASS_PARTIAL_READBACK/exit2仍表示路径BLOCKED；
不能从元数据PASS推导物理整段PASS。新增读取器的合成测试不是力学解。
本轮accepted时间仍.05s。若事件定位失败，REJECTED候选可能保存target而非最后二分点；
必须看专用诊断trace，不得猜测candidate时间代表失败二分时间。
专用事件重放不允许commit，classification=PURE_EVENT_REPLAY_REJECTED_DIAGNOSTIC_NOT_RESTART，
status=DIAGNOSTIC_COMPLETE_NOT_PATH_COMPLETE、restart_allowed=false；不能用于生产重启。
diagnostic.json包含37次branch及可行性检查、独立向量判据、接受态前后指纹和源hash；
inferred_bisection是从观察trace推算的边界，不是原异常保存的局部变量。

## Q7 existing J2 + vector friction companion (2026-10-01)

outputs/process_plastic_friction_q7_20261001/n3_i4 is a real PARTIAL coupled
benchmark, status=BLOCKED at accepted.05s/step6, not the requested.2s endpoint.
n3_i8/NOT_RUN.json is explicit: no common-end increment comparison was performed.
Classification ENGINEERING_APPROXIMATE_QUASISTATIC_J2_TWO_TANGENT_FRICTION_BENCHMARK,
production_authorized=false; no FKS/AA6061/thermal/AFSD certification.

Same Q6 fixed-normal geometry/IDs/two-vector contact history and independent
contact work definitions. Q7 replaces Q6 static finalep[N] with material histories
[F,N] (F=accepted_steps+1, N27): equivalent_plastic_strain, plastic_work_J,
return_dissipation_diagnostic_J, interval_plastic_work_J, interval_return_diagnostic_J,
equivalent_strain_rate_per_s. No rejected candidate is added to these frames.
The initial values0 mean actual zero initial state, not fabricated missing fields.
Trial/Newton/probe/bisection increments never enter accumulated P/R until commit.

P is actual constitutive q_new*dep*V. For the declared constant-Y90MPa case,
independent P_i=Y*V_i*ep_i (initialep0); interval return R_i=1.5*G*V_i*(dep_i)^2,
G35GPa. R is algorithmic return diagnostic, NOT extra physical dissipation or heat.
Friction Q is indexed by contact [F,M9], not by all material particles. Do not add
R or a material work defect to P/Q to force an energy balance.

records add interval_plastic_work_J/interval_return_diagnostic_J and
material_work_warning. The latter is abs(Wstress-Es-P)>1e-6J, a preregistered
diagnostic, NOT certification or posthoc force/temperature modification.
report.material_work_warning_steps enumerates all flagged accepted frames.
Actual n3_i4 has P=.2396461244J,Q=.00528460208J and material defect-.00428535J
at.05s. Partial readback PASS verifies faithful accounting, not zero energy defect.
Temperature293.15K remains prescribed input; no temperature/dynamic velocity field
is exported to VTK/history. Quasi-static interval rates are not inertial velocities.

Before solve, provenance.json stores config/source hashes. Every accepted state
also gets unique accepted_####.npz: same frame arrays under the driver's short
names, persistent particle_id/tracer_type,step,physical_time_s,record_json,
classification=ACCEPTED_DIAGNOSTIC_NOT_RESTART,restart_allowed=false.
These preserve prior states on an external timeout; an incomplete final file is
invalid and cannot supersede an earlier verified file. Complete final reports hash
every journal and the reader compares them against history/records/IDs exactly.
checkpoint.npz also remains diagnostic-only, not production restart support.

REJECTED_manifest records attempted_time_s=.0625,lastaccepted=.05/step6,
candidate time/base step,exception,Newton iteration trace; REJECTED_candidate.npz
contains candidate arrays separately. final.vtk/checkpoint contain ACCEPTED .05s,
never the attempted.0625s displacement/stress/ep. Reader returns
PASS_PARTIAL_READBACK with exit2 for this honest BLOCKED run; it cannot promote it
to COMPLETE. Finer case NOT_RUN is not a zero-field result.

branch_diagnostic_12_14_layout/ is a separate PURE_ACCEPTED_STATE_BRANCH_DIAGNOSTIC_NOT_RESTART.
Both REJECTED_diagnostic_budget12/14.npz remain REJECTED_DIAGNOSTIC even though
budget14 returned a branch root on correction13. The forced-stick ID21 contact
still violates its circle, so event_split_required=true; accepted time stays.05s.
No diagnostic candidate can be concatenated onto history or displayed as an
accepted physical frame. Original blocked case source/config and output remain intact.

## Q6 vector tangent circular cone companion (2026-10-01)

outputs/process_friction_cone_q6_20261001/{oblique_i4,oblique_i8,basis37_i4}:
ENGINEERING_APPROXIMATE_QUASISTATIC_TWO_TANGENT_CONE_BENCHMARK, NOT_AFSD,
production_authorized=false. TWO_TANGENT_CONE_PATH_COMPLETE only means the
declared manufactured program reached.2s. Fixednormal-X,9X+wall contacts,
27deformable points,54freeDOFs,clampedX0. BOTH Y/Z have circular Coulomb friction.
The prescribed wall tangent program has global direction(.6,.8); material
displacements and traction directions are solved, not set to this wall direction.
37degree tangent basis is a coordinate representation change, not a rotated tool.

Keep Q5 identity/reference/position/stress/mass/support/normalforce/work fields.
Replace scalar tangent fields with these GLOBAL YZ two-vectors [F,M,2]:

- tangent_force_global_yz_N;
- interval_relative_slip_global_yz_m;
- tangent_impulse_global_yz_N_s;
- relative_kinematic_rate_global_yz_m_s;
- contact_law_residual_global_yz_N.

interval_modes[F,M] now0=zero slip within1e-11m,1=nonzero slip over ENDING interval;
it is neither slip sign nor the Newton release mask. Forces satisfy a disk,
not independent rectangular component bounds. released_rows[M] bool in events
means natural-map rows were solved; a released row may end up sticking.
events: kind=CONE_LIMIT,step,time,bracket_s,ending_interval_modes,released_rows,
cone_margin_N=mu*normal-norm(tau). Ordered events match accepted ledger frames.

Q_i=-dot(tau_mid_i,du_yz_i-dwall_yz)=-dot(Jt_i,relative_interval_rate_i).
Raw negative roundoff tails remain stored, without abs/clip/cancellation.
Rates are quasi-static interval kinematics, not solved dynamic velocity.
contact_law_residual uses disk projection only to evaluate a root equation;
saved physical tau comes from material equilibrium, never from posthoc force clipping.
At commit/readback, independent cone magnitude, natural-map residual and
dimensionally propagated vector-flow error are checked; see Q6CONTRACT addendum.
Nearzero slip directions cannot be interpreted more accurately than slip tolerance.

records add maximum_contact_law_residual_N,minimum_cone_margin_N; inherit explicit
wall/material/support/stress work and elastic diagnostic. Final ep is0 because
this benchmark actually remained elastic; plasticity+friction is NOT_SOLVED.
Final compressed/friction-supported displacement is not stress-free plastic set.
Manufactured input temperature293.15K is not a solved thermal field.

Same NPZ/VTK accepted-state and hash rules as Q5. checkpoint is diagnostic-only,
restart_allowed=false; REJECTED_candidate/manifest never becomes a physical endpoint.
History stores27persistent IDs and observer tags; output basis rotation does not
change them. No extra solved temperature, dynamic velocity or inferred deposition.

Each *_readback.json outside the run folder is independently regenerated from
saved data, not solver status alone. Synthetic reader test fixtures explicitly
say SCHEMA_FIXTURE_NOT_PHYSICAL; they cannot serve as mechanics evidence.
comparison.json compares only COMPLETE paths at the same .2s and common requested
times; no interpolation across event states or comparison of unequal failed endpoints.
4/8increment differences are sensitivity, not a demonstrated convergence order.
selected_array_bytes excludes Python/PySPH process overhead and is NOT peak RSS.

## Q5 directional-Y friction companion (2026-10-01)

outputs/process_friction_q5_20261001/n3_i4 and n3_i8 are manufactured NOT_AFSD.
Classification ENGINEERING_APPROXIMATE_QUASISTATIC_DIRECTIONAL_FRICTION_BENCHMARK.
production_authorized=false; DIRECTIONAL_FRICTION_PATH_COMPLETE means the stated
quasi-static program completed, not isotropic Coulomb, paper or production validation.
Only X=0 face clamped all3 components; no inherited Y/Z rollers. X+ normal stays
compressed; unsupported normal detachment must BLOCK. Friction acts only inY, Z free.

history.npz: particle_id/tracer_type[N], reference_position_m[N,3], mass_kg/reference_volume_m3[N],
support_components[N,3], contact_particle_id[M]/contact_feature_id[M]=1 (M9,N27),
physical_time_s[F], position_m/displacement_m[F,N,3], stress_Pa[F,N,3,3],
normal_force_N/tangent_force_N[F,M], interval_modes[F,M], contact_force_N/support_force_N[F,N,3],
friction_dissipation_J/interval_friction_J[F,M], interval_relative_slip_m[F,M],
tangent_impulse_N_s/relative_kinematic_rate_m_s[F,M], kinematic_interval_rate_m_s[F,N,3].
equivalent_plastic_strain[N] is an elastic-only final diagnostic (verified0), not a
missing plastic history disguised as zeros. Plasticity-with-friction NOT_SOLVED.

lambda=Fint_x positivecompression, tangent_force=tau onmaterial, gap=wallX-uX.
mode0 denotes sticking over the ENDING accepted interval; ±1 denotes sign of
material-relative-wall slip inY, not wall motion or force sign. A wall reversal
may produce nonzero frictionforce and zero slip during elastic sticking/unloading.
Normal point is (L+wallX,current_y,current_z), tangent axis+Y, normal-X. No arbitrary
rotating normal or two-component friction cone is implied by 3D particle positions.

Per interval: slip=duY-dwallY; Jt=tau_mid*dt; relative_kinematic_rate=slip/dt;
rawQ=-Jt*relative_kinematic_rate. These are quasi-static interval integrals, not
inertial kick midpoint velocities. Raw small negative Q values remain stored;
initial interval arrays0 are initialization, not clipping. Dissipation is NOT a
solved heat/temperature field. No thermal equation or heat partition is run.

records distinguish wall_work_J, contact_material_work_J, support_work_J, stress_work_J,
friction_dissipation_J, interval_friction_J, minimum_raw_contact_dissipation_J,
friction_work_identity_J=Wwall-Wmaterial-sumQ, Es and stress-work/Es residual.
wall_displacement_m[3] describes prescribed rigid boundary only; wall_force_on_material_N[3]
is not confused with total reaction_N (wall+support). mode_counts and raw minnormal
are diagnostic. No residual energy is added toQ. Final normal compression remains,
so nonzero storedenergy is not a closed unloadedcycle error. Residual tangential
displacement with ep0 is a friction-supported elastic state, not plastic deformation.

events: step/time, kind=FRICTION_LIMIT, bracket_s, ending_interval_modes[M], cone_margin_N[M].
Ordered event steps match ledger FRICTION_LIMIT frames one-to-one, no missing/duplicates.
Program corners .05/.1/.15 and final.2 are requested states. Extra events add frames;
F=accepted_steps+1, not always4*increments_per_phase+1. Bracket width is internal
search resolution, not an analytic event error bound. Candidate probes are not frames.

final.vtk contains accepted IDs/tags/position/displacement/stress/ep, no temperature
or dynamic velocity. checkpoint.npz is diagnostic-only restart_allowed=false;
PySPH scratch T/u/v/w/e cannot become solved output. REJECTED states stay separate
with manifest lastaccepted/candidate/attempted time; never promote them to endpoints.
Source/config/artifact hashes and evidence freeze preserve old failures. Inspector
independently reconstructs saved contact work/Q and Es, but does not independently
rebuild the body gradient or pointwise force. See PROCESS_FRICTION_Q5_REPORT.md.

## Q4 isolated orthogonal two-wall benchmark (2026-10-01)

outputs/process_corner_q4_20261001 is NOT_AFSD, production_authorized=false.
Classification ENGINEERING_APPROXIMATE_QUASISTATIC_SMALL_STRAIN_TWO_WALL_BENCHMARK.
Four actual paths finish0.2s; CORNER_PATH_COMPLETE is program completion, NOT all
scientific tests PASS. A strict unloading ep-constancy FAIL remains in the core suite.
Do not replace it with COMPLETE or the independent readback PASS status.

history.npz retains Q3 material geometry/stress/ep/energy/ID/time arrays, but normal
contact is indexed by M=18 rows: contact_particle_id[M], contact_feature_id[M] (1=X,
2=Y), contact_axis[M] (0/1), contact_normal[M,3] (-ex/-ey), normal_force_N[F,M],
active_contacts[F,M], gap_m[F,M]. There are27 material IDs,15 distinct contacting
particle IDs and18 unique (particle_id,feature_id) pairs;3 edge particles each have
two features. Repeated contact_particle_id is intentional, not duplicated mass.
contact_force_by_wall_N[F,2,N,3] replaces Q3 contact_force_N. Sum over walls for the
particle force; do not overwrite shared-edge contributions or double-count full vectors.

positive lambda means compression; force=-lambda*e_axis; gap=d_wall-u_axis.
Wall-contact point is the material point projected to x=L+dX or y=L+dY; normals
are fixed. Candidate set is initial boundary only, not general collision detection.
X/Y walls use symmetric or staggered programs recorded in resolved_config/report;
wall_displacement_m is never a free-material trajectory. Z+ is a free material surface.
Only X=0/Y=0/Z=0 normal rollers are permanent. displacement/position include all3 DOF.

records wall_force_on_material_N[2,3], wall_displacement_m[2], wall_work_J[2],
contact_material_work_J[2], interval_wall_work_J[2], interval_contact_material_work_J[2],
wall_minus_material_work_J[2], active_contacts_by_wall[2], minimum_gap_by_wall_m[2],
maximum_gap_by_wall_m[2], minimum_normal_force_by_wall_N[2] retain per-wall values.
support_work_J/stress_work_J/plastic_work_J/return_dissipation_diagnostic_J and Es
are distinct global ledgers. A wall and material work discrepancy cannot cancel
against the other wall or become heat. reaction_N remains summed wall+support force,
not either wall alone. No ambiguous cap_displacement/cap_reaction/actuator_work fields.

events carry step, physical_time_s, kind=RELEASE_LIMIT, bracket_s, near_zero_contact_rows,
near_zero_feature_ids. The ordered event steps must match ALL RELEASE_LIMIT ledger
frames one-to-one, without omission/duplicates/reordering. Bracket width is a search
tolerance, not total analytic-time error. Every accepted internal event is an extra
history frame; F=accepted_steps+1, not necessarily2*requested_increments+1.

final.vtk and checkpoint.npz contain accepted states only; diagnostic snapshot
restart_allowed=false. Upstream scratch velocity/internal energy/temperature must
not be published as solved physics. kinematic_interval_rate_m_s=delta_u/physical_dt
is not inertial momentum velocity. Finite-strain density, friction, heat, FKS contact,
inertia, general rotating tool, AFSD and production restart are NOT_SOLVED.
Raw small negative force residues or small ep increments are retained, not clipped.

REJECTED_manifest/candidate policy follows Q3, with outer attempted time distinct
from candidate time/dt/base_step and last accepted state; no rejected state is a
physical endpoint. Constructor failures remain outside catch. Source/config/artifact
hashes are saved, and the new freeze preserves all previous failures. The read-only
inspector rederives per-wall work and material diagnostics independently; this checks
the stored benchmark, not paper validation. Details: PROCESS_CORNER_Q4_REPORT.md.

## Q3 isolated unilateral platen benchmark (2026-10-01)

outputs/process_contact_q3_20261001 is NOT_AFSD. Classification:
ENGINEERING_APPROXIMATE_QUASISTATIC_SMALL_STRAIN_UNILATERAL_CONTACT_BENCHMARK.
elastic_n3_i8/j2_n3_i8/j2_n3_i16 reach the manufactured0.2s endpoint. An internal
RELEASE_LIMIT adds an accepted frame; accepted_steps is not always2*increments_per_phase.
Time is physical loading-program time, not an AFSD result or inertial integration.

history.npz includes persistent particle_id/tracer_type[N], reference_position_m[N,3],
reference_volume_m3/mass_kg[N], support_components[N,3] (permanent rollers only),
contact_particle_id[M], contact_feature_id[M] (single plane=1), physical_time_s[F],
position_m/displacement_m[F,N,3], stress_Pa[F,N,3,3], equivalent_plastic_strain,
plastic_work_J/return_dissipation_diagnostic_J[F,N], contact_force_N/support_force_N[F,N,3],
normal_force_N/active_contacts/gap_m[F,M], kinematic_interval_rate_m_s[F,N,3].
That rate is delta_u/physical_dt, NOT momentum velocity. Temperature/heat/current
density evolution are NOT_SOLVED. Reference density and prescribed T remain inputs.

gap=d_wall-u_x; normal_force_N=lambda positive compression; material contact=-lambda*e_x.
Active/open at exactlyzero gap/zero force is a degenerate numerical classification,
not evidence of sticking/sliding. Raw signed force/gap values retain iteration or
roundoff residues; no clipping. Initial boundary nodes only: NOT general collision detection.
Contact point can be reconstructed on wall as (L+d, current y, current z); normal is
constant(-1,0,0). Tangential motion is unrestrained and carries no contact force.

report.records carry step/time, wall_displacement_m, wall_force_on_material_x_N,
wall_work_J, contact_material_work_J, support_work_J, stress_work_J, raw plastic/return
and small-strain elastic diagnostics, separate interval_wall_work_J and
interval_contact_material_work_J, cumulative wall_minus_material_work_J, original
internal work identity, gap/normal force extrema, complementarity magnitude[J],
active count, Newton/active-set trace and event_kind. Parent cap_displacement_m is
the WALL program, not material cap displacement; actual material displacement is in history.
Parent reaction_N is summed wall+support force, not either boundary alone.
Parent accepted_rate_floor_evaluations=0 is inherited scratch diagnostic: Q3 never uses FKS.
Torque is computed with reference X only, not finite-rotation angular momentum verification.

report.events: accepted step/time, kind=RELEASE_LIMIT, bracket_s[lo,hi]. Time bracket
width is not a bound on total error relative to an analytic event; equilibrium tolerance
also contributes. Saved events and full history are authoritative, not uniform animation time.
UNILATERAL_PATH_COMPLETE means only the stated boundary program completed; production=false.
BLOCKED keeps last accepted state plus REJECTED_manifest and available candidate arrays.
Manifest distinguishes outer attempted_time_s from candidate_time_s/candidate_dt_s/
candidate_base_step; candidate_available=false if no candidate formed. Constructor exceptions
are outside the runner catch; no universal forensic guarantee. Never restart rejected data.

final.vtk: actual accepted position, reference position/displacement, IDs/labels,
six stress components, ep and plastic work. No fabricated thermal/dynamic field.
checkpoint.npz is upstream diagnostic-only restart_allowed=false; scratch u/v/w/e/T
must not be published as solved dynamics/internal-energy/temperature. Contact state is
authoritative in history, not a complete restart interface inside this diagnostic NPZ.

report.source_sha256/artifact_sha256 plus PROCESS_CONTACT_Q3_FREEZE preserve evidence.
Reader independently recalculates wall program, geometric gap, unilateral conditions,
actual force-displacement work, J2 plastic/return formulas, elastic diagnostic, IDs/labels
and record/frame correspondence. PASS_COMPLETE_READBACK is data consistency only;
readable BLOCKED data stays PASS_PARTIAL_READBACK. Selected retained array bytes exclude
temporary buffers/RSS. Full process timing lives in logs, separate from advance_wall_s.

## Q2 isolated plastic compression-return benchmark (2026-10-01)

outputs/process_plastic_q2_20261001 is NOT_AFSD. Classification:
ENGINEERING_APPROXIMATE_QUASISTATIC_SMALL_STRAIN_PLASTIC_BENCHMARK.
Four cases j2_n3_i8/j2_n3_i16/j2_n5_i8/fks_n3_i8 reach a manufactured0.2s program,
not AFSD process time. Return of an axial actuator is NOT a free unloaded shape.
Production_authorized=false; no solved thermal/contact/inertia/density-evolution field.

history.npz: particle_id[N], passive tracer_type[N], reference_position_m[N,3],
reference_volume_m3[N], mass_kg[N], constrained_components[N,3], cap_particles[N];
physical_time_s[F], position_m/displacement_m/reaction_N[F,N,3], stress_Pa[F,N,3,3],
equivalent_plastic_strain/plastic_work_J/return_dissipation_diagnostic_J[F,N].
plastic_work_J is accumulated accepted V*q_new*dep, not solved heat.
return_dissipation_diagnostic_J=accumulated V*1.5G*dep^2, not a residual-heat source.
kinematic_interval_rate_m_s[F,N,3]=delta_u/physical_dt is not momentum velocity;
deviatoric_strain_rate_per_s[F,N] is the actual material input rate before FKS floor.
Fixed initial IDs/labels are not recomputed from deformed positions.

report.records[step] stores physical_time_s, cap_displacement_m/cap_reaction_x_N,
all support reaction_N, actuator_work_J/stress_work_J, raw plastic/return/elastic
diagnostic energies[J], work_identity_residual_J and work_minus_elastic_minus_plastic_J,
ep/rate extrema, accepted_rate_floor_evaluations, free-force residual/Newton trace,
mass_kg, normalized internal-force/torque defects. Torque is about reference X,
NOT a current-configuration finite-rotation angular momentum certificate.
Es=V*(S:S/(4G)+K*(trace grad u)^2/2) is a small-strain diagnostic only.

report.status=DISCRETE_PLASTIC_PATH_COMPLETE refers only to accepted discrete
equilibrium at the prescribed endpoint; errors retain BLOCKED and last accepted
history. REJECTED_manifest/candidate arrays are debugging evidence, not final physics.
Only exceptions carrying a candidate create REJECTED_candidate.npz; early ValueError
may have a manifest only. This is not a full restart/forensic guarantee for all faults.
resolved_config.json includes input sources/units/status, law, prescribed temperature,
unchanged E2 geometry/elastic inputs and actual FKS parameters. FKS rate_floor1/s and
dimensionless b1 remain explicitly DEBUG; no temperature, modulus or mass scaling.

final.vtk publishes accepted geometry, IDs/labels, displacement/reference position,
six stress components, ep, plastic/return ledgers and reference_density_kg_m3 only.
checkpoint.npz is an upstream diagnostic dump with restart_allowed=false: its
ParticleArray scratch u/v/w/e and prescribed temperature/reference rho must NOT be
rendered as solved dynamics/internal energy/thermal/density evolution. Use history
and schema-declared VTK fields, not arbitrary scratch properties.

setup_wall_s/advance_wall_s are separate; wall_before_report_s includes output writes
before final report. Full interpreter/import/serialization time is in command logs.
selected_array_bytes excludes transient buffers and RSS, not peak process memory.
source_sha256 identifies direct inputs; artifact_sha256 hashes all pre-report files.
PROCESS_PLASTIC_Q2_FREEZE retains transitive parent dependencies and current evidence.
*_readback_v2.json independently checks hashes, frame/ledger correspondence, actual
reaction work, raw monotonic accepted plastic ledgers, independent J2/FKS work and
return formulas, IDs and output identity. PASS_COMPLETE_READBACK is not AFSD
validation; readable BLOCKED history stays PASS_PARTIAL_READBACK (CLI exit2).
Old *_readback.json and failing regression logs remain historical, not overwritten.

## Q1 isolated quasi-static benchmark (2026-10-01)

outputs/process_time_q1_20261001 is NOT_AFSD. Classification is
ENGINEERING_APPROXIMATE_QUASISTATIC_SMALL_STRAIN_BENCHMARK; no production authorization.
Each of n3_i4/n5_i4/n5_i8/n7_i4 reaches the manufactured0.1s ramp, not0.1s of AFSD.
history.npz stores persistent particle_id/tracer_type, reference_position_m,
reference_volume_m3, mass_kg, constrained_components[N,3], physical_time_s[frames],
actual position_m/displacement_m[frames,N,3], stress_Pa[frames,N,3,3],
external_force_N/support_reaction_N[frames,N,3] and actual ep(final,N).
kinematic_interval_rate_m_s is displacement increment/physical interval; it is
NOT velocity obtained from momentum integration. Thermal/inertia/contact/FKS are NOT_SOLVED.

final.vtk contains actual coordinates, displacement, reference position/density,
ID/passive label, six solved stress components and ep; no temperature or dynamic velocity.
checkpoint.npz reuses upstream diagnostic dump, restart_allowed=false: scratch zero
u/v/w and prescribed temperature/reference rho inside ParticleArray MUST NOT be
published as solved dynamic/thermal/density fields. Production restart is unavailable.

report.status=DISCRETE_EQUILIBRIUM_COMPLETE means only the discrete bounded load
problem reached its endpoint. continuum_accuracy is reported, not a precision gate;
non-monotone stress refinement does NOT become PASS. Records separate external_work_J,
stress_work_J, constraint_work_J and elastic_small_strain_diagnostic_J. No residual heat.
wall_seconds covers setup+advance before files are written; full subprocess timing
is separately recorded in n7_i4.log. array_bytes is selected arrays, NOT peak RSS.
comparison.json is READBACK_ONLY: hashes, common endpoint, independent load/actual
force-displacement integral, 4/8 increment differences. It is not AFSD convergence.
Source hashes in each run are a direct subset; PROCESS_TIME_Q1_FREEZE retains
complete previous dependencies and all current artifacts. Old schemas/FAILs unchanged.

## E3-4 accepted chain and read-only mechanisms (2026-10-01)

Physical VTK/checkpoint schema is unchanged. The two accepted segments at
outputs/short_window_e3_20261001/dx05_cfl01_05us and
outputs/continuation_e3_20261001/to_2us reach absolute2e-6s/step191.
Derived chain_report.json contains192unique frames, one exactly checked step48
seam removed, original real time grid and hashes. No interpolation or solver advance.
reached_target requires both manifest and runner COMPLETE in addition to accepted
time. A readable FAILED/INCOMPLETE chain is not relabelled as completed physics.
observed_segment displacement is from that segment's first frame;
observed_origin_to_end is only available if actual step0/time0 exists, otherwise null.
Cumulative ledgers are not summed across restarts: segment deltas subtract the
parent endpoint. Supplied preflight_receipt.json verifies before-run parent hashes.
Unsaved child initial internal/scratch checkpoint equivalence is NOT established.

mechanism_report_v2.json schema_version1.1 is ACCEPTED_STATE_MECHANISM_OBSERVATION.
It is authoritative over the retained first mechanism_report.json, whose confusing
net_feed_origin_below_layer_mass_kg name actually meant initially ABOVE the layer.
The corrected field is net_initially_above_layer_currently_below_mass_kg [kg]:
component1,z0>=layer,current z<layer. It is a net position count, not gross mass flux,
and is separate from persistent final_ever_feed_plane_crossed_mass_kg.
outside_bore_layer_candidate_* selects original bar/component1 AND actuator0,
particle centres0<z<=layer outside the currently translated/rotated square bore.
It does not select all substrate or prove newly fed/deposited material.
volume_m3=sum(m/rho), never a2D area ratio or Fig8 point-count surrogate.
free_plastic_* records both raw ep>0 and observation-only ep>1e-10 counts;
IDs are sorted for readback only. No force/heat depends on tracer/AS/RS labels.
Plastic work [J] is read from solved cumulative per-particle plastic_work;
negative raw counts remain visible, never clipped. Not-solved fields stay absent.

cost_extrapolation uses measured diagnostic dt median and totalwall/newacceptedsteps,
with SI prescribed rates; values are constant-current-cost estimates, NOT promises
or long-run authorization. Initial/plastic/deposition/production claims are separate:
this chain has3free plastic bar particles but no observed feed crossing or candidate
outside-bore layer volume. paper_comparison/sustained_deposition_validated and
production_authorized remain false. New methods in PROCESS_TIME_DECISION are NOT_RUN.

## E3-3 explicit short-window acceptance (2026-10-01)

outputs/short_window_e3_20261001 contains one preregistered .5us implementation
window, NOT the original2us goal or deposition validation. Existing schema and
solved fields are unchanged; runner/manifest COMPLETE refers only to5e-7s.
The derived integration_report.json records48accepted steps/49frames, no rejected
state, and keeps paper/deposition/mixing/production claims false.
tests/review/test_short_window_e3.py checks every frame, source/config identity,
saved fields/checkpoint, SI tool pose, passive IDs and the existing energy ledger.
The .5us final shortened step is explicit; performance.dt_min_s must not be
mistaken for stable CFL dt or a valid paper-scale cost basis.
Any future continuation must retain checkpoint_parent/source/config fingerprints
and the actual accepted time grid; no historical REJECTED state is migrated.

## E3-1 bounded integration evidence (2026-10-01)

ENGINEERING_INTEGRATION_E3.md authorizes one engineering-approximate short probe,
not production or quantitative validation. Existing VTK/checkpoint fields and
all physical parameters remain unchanged. No source change occurred in this run.

scripts/inspect_engineering_integration.py reads the actual series/NPZ/VTK,
validates source/config/frame hashes and persistent ID/labels. It selects the
accepted endpoint from series.json/checkpoint.npz, never from the rejected
candidate or performance.completed_steps alone. An alarm can increment the
runner count before rejection; callers must not mistake that count for acceptance.

outputs/engineering_e3_20261001/integration_report.json is a derived READBACK:
status=INCOMPLETE_COUPLED_IMPLEMENTATION, readback_status=PASS are independent.
It records accepted_step/time/frames, requested_end, real free-particle movement,
plastic count, temperature, actual tool angle/travel and prescribed cap travel,
plus rejection identity/features and certificate metrics. Units are SI.
The geometric cap travel is NOT an inferred free-material displacement.
Zero plastic strain/feed-plane crossing is reported as observed, not filled with
an expected deposition pattern. No mixing or deposition claim is authorized.

This attempt accepted21steps/22frames; candidate22/particle476 is REJECTED.
The candidate's auxiliary trace can belong to a previous callback (existing
trace_not_guaranteed_to_be_current_candidate=true); use failed_local_context
for candidate identity, not the stale auxiliary trace. interval_context is not
propagated by the existing wrapper; if absent, it is UNAVAILABLE, not reconstructed.
Local high-precision cone reference is an offline numerical diagnostic, not a
replacement accepted state and never a restart file or reconstructed trajectory.

## E2 engineering block benchmark (2026-10-01)

ENGINEERING_BLOCK_E2.md defines a separate small-strain, fixed-reference-geometry
benchmark; it does NOT replace the production fields or material parameters below.
classification=ENGINEERING_APPROXIMATE_SMALL_STRAIN_REFERENCE_CONFIGURATION,
production_authorized=false. This is force-driven three-component motion on a
three-dimensional NNPS graph, not an AFSD deposition run or a paper comparison.

history.npz stores persistent particle_id, reference_position_m [N,3], volume_m3,
mass_kg, time_s [frames], actual position_m/velocity_m_s [frames,N,3], full
stress_Pa [frames,N,3,3], final equivalent_plastic_strain and body_force_N [N,3].
reference_final_displacement_m/reference_final_stress_Pa are analytical references,
not solved values. Tracer tags are passive; only initial-coordinate boundary mask
selects the roller's normal DOF. No solved temperature/heat/friction is published.

report.json records actual time/dt/steps, setup/advance wall seconds, array byte
counts (NOT process peak memory), source/config/artifact SHA256, all errors/checks,
and per-step work/energy/force/momentum. end_state_reaction_x_N is the instantaneous
force needed at the accepted endpoint; last_kick_midpoint_reaction_N is the actual
last kick reaction. constraint_impulse_N_s is accumulated from the actual kicks;
body_impulse_N_s and linear_momentum_kg_m_s allow an independent momentum check.
constraint_work_J is from the actual impulse/midpoint velocity, not a reference
force substituted into the ledger. Elastic energy is a limited small-strain
diagnostic; current production Jaumann variable-volume issues remain unresolved.

final.vtk contains actual points, velocities, IDs, stress_xx and ep only plus
reference positions/density. reference_density_kg_m3 is a fixed integration
weight property; it is not a solved finite-strain density field. No e=0 scratch
field from the ParticleArray may be interpreted as physical internal energy.

checkpoint.npz uses existing upstream storage but restart_allowed=false: it is
a diagnostic snapshot, not a complete reference-operator/stress-staging restart.
The existing restore() rejects it, by design. Restart capability for this new
benchmark is NOT_IMPLEMENTED, not silently borrowed from production checkpoints.
PILOT_NOT_ACCEPTANCE is timing only. FAIL can be a completed run exceeding
preregistered physical accuracy; raw fields must retain that status via report.json.

## Independent review override — 2026-09-24

Current scientific authority: SCIENTIFIC_REVIEW.md and VALIDATION_CONTRACT.md.
The historical HOST_FKS_DEBUG profile below is preserved, not a current acceptance.
Current variant is REVIEW_APPROX_CORRECTED_WORK_CONJUGATE_CONTACT_MIDPOINT,
classified APPROXIMATE, never production-authorized. Four sensitivity attempts
under outputs/review/ are FAILED; COMPLETE elsewhere only means a requested
computational interval finished, not quantitative validation.

Manifests add validation_contract(version/path/sha256/production_authorized=false)
and model_fidelity. VTK filename/title uses APPROXIMATE; all solved fields and
signs remain as below. performance.json adds requested_end_time_s; an unreached
target is INCOMPLETE or FAILED. In the four review attempts checkpoint/VTK only
contain step0; accepted later-step scalars exist in diagnostics.jsonl, not as
invented particle time series. A failed partial kick is never exported as valid.

Important: internal_energy_J_per_kg / elastic_shear_J are derived diagnostic
quantities, NOT a validated thermodynamic state function at finite volume change.
The independent prestressed-compression counterexample fails (review N08).
Do not cancel this defect with a bookkeeping source or label residual as heat.
Friction heat now uses the actual kick midpoint; negative dissipation fails.

Postprocessing measurement.py outputs R_count/R_volume only, preserving the
complete explicit sampling spec, its SHA256 and selected IDs. experimental_area_ratio
is null/BLOCKED. The function checks sampling, not input-file provenance: callers
must first verify SCIENTIFIC_FREEZE.json, run manifest/config/source and VTK hashes,
and match sampling time/rotation/traverse to that run. No valid experimental ROI
is supplied in this release. See VALIDATION_CONTRACT for all operator edges.

outputs/review/sensitivity_summary.json is a read-only derivative with input hashes,
geometry/BC equality checks, diagnostics extrema and hardware cost projections.
Different terminal times force N09/N10 comparison values to null/BLOCKED.
No new numerical fields or experiment reference CSVs are synthesized.

## Scope

This contract separates particle identity from mutable state and separates
physics runs from geometry/debug artifacts. The machine-readable counterparts
are `schemas/run-manifest.schema.json` and
`schemas/particle-record.schema.json`. The dependency-free VTK writer now serves
both the independent geometry fixture and the tested PySPH HOST_FKS_DEBUG
adapter. These remain distinct output profiles/classifications.

All stored calculation fields use SI. A viewer may convert metres to
millimetres, radians per second to rpm, pascals to MPa, or kelvin to degrees
Celsius, but display conversions must not overwrite stored values.

## Run directory and manifest

Each run directory has one primary manifest: the fixture profile uses
`<run_id>.manifest.json`; the HOST_FKS_DEBUG profile uses `manifest.json`.
The host runner saves an initial checkpoint/state and then an INCOMPLETE manifest
before entering the time loop. JSON metadata is replaced atomically and finalized
as COMPLETE, INCOMPLETE (wall cap, exit 2), or FAILED (alarm, exit 1).
Configuration/pin failures occur before creating a run directory and fail nonzero.
Abrupt process death may leave INCOMPLETE; it is never proof of completion.

Required manifest data:

| Group | Required content |
|---|---|
| Identity | `schema_version`, unique `run_id`, `created_utc` |
| Classification | `run_classification`, `physics_claims_allowed`, `status` |
| Backend | repository/name, exact commit, dirty flag, compiler/build mode |
| Configuration | case ID, resolved-config path and SHA-256; unresolved override list if any |
| Units | SI system and base field units |
| Coordinates | stored frame, axis definition, tool translation/rotation convention, AS/RS rule and paper-figure mapping status |
| Particles | file/series, format, count, SHA-256, canonical field names |
| Provenance | command, log paths, checkpoint parent, code/data origins and known approximations |

The backend commit may be `null` only for non-backend fixtures. A physics run
must record a pinned commit and configuration digest. Restarted runs add parent
run/checkpoint IDs; they do not reuse another run's ID.

## Run classifications

| Value | Meaning | Permitted claim |
|---|---|---|
| `GEOMETRY_FIXTURE` | Synthetic/static geometry, label or I/O test | No physics or reproduction claim; `physics_claims_allowed=false` |
| `DEBUG` | Short solver diagnosis, possibly coarse or incomplete | Only the explicitly tested code path |
| `APPROXIMATE` | Physics run with named approximations/calibration | Qualitative result within documented assumptions |
| `PAPER_COMPARISON` | Configuration cleared through all paper-comparison gates | Comparison only; never automatically means successful reproduction |

Filename text is not authoritative: the manifest classification is required.
However fixture/debug filenames must also include their classification so that a
detached image cannot be mistaken for a result.

## Particle identity and immutable fields

The following logical fields are mandatory in every state file or must be
losslessly joinable from a manifest-declared immutable table:

| Field | Type / shape | Unit | Rule |
|---|---|---|---|
| `particle_id` | unsigned integer scalar | `1` | Globally unique, persistent through sorting, restart and inflow; never reuse an emitted ID |
| `initial_position_m` | float vector[3] | `m` | Position at particle creation; immutable |
| `material_id` | integer scalar | `1` | Mechanical/thermal material identity |
| `tracer_type` | integer scalar | `1` | Immutable initial tracer identity, never recomputed from current position |
| `component_id` | integer scalar | `1` | Substrate, feedstock, tool, base or later explicit component code |

Oxide and copper proxy tags in this project are passive: tagged feedstock points
retain the same AA6061-T6 `material_id` and material state as untagged feedstock.
Force, stress, thermal, contact and constitutive code must not branch on
`tracer_type`.

Neighbor sorting must apply one permutation to all particle fields. The test
suite must compare the `(particle_id, initial_position_m, material_id,
tracer_type)` mapping before and after sorting and after checkpoint restart.

## Time-dependent physics fields

Physics states are expected to expose these canonical names. A field may be
omitted from `GEOMETRY_FIXTURE`; a backend capability/run manifest must explain
why it is absent from a physics run.

| Field | Shape | Unit / ordering |
|---|---|---|
| `position_m` | vector[3] | `m`, coordinates in the manifest's stored frame |
| `velocity_m_per_s` | vector[3] | `m/s` |
| `mass_kg` | scalar | `kg` |
| `density_kg_per_m3` | scalar | `kg/m^3` |
| `temperature_K` | scalar | `K` |
| `pressure_Pa` | scalar | `Pa`; pressure sign convention in manifest |
| `stress_Pa` | symmetric tensor[6] | `Pa`, order `xx, yy, zz, xy, yz, xz` |
| `strain_rate_per_s` | symmetric tensor[6] | `1/s`, same order |
| `equivalent_plastic_strain` | scalar | `1`, non-negative |
| `equivalent_plastic_strain_rate_per_s` | scalar | `1/s` |
| `internal_energy_J_per_kg` | scalar | `J/kg` |
| `smoothing_length_m` | scalar | `m` |

Per-state metadata includes `step`, `time_s`, `dt_s`, live particle count,
tool pose/velocity, total mass, total linear/angular momentum, energy accounting,
minimum/maximum temperature, neighbor statistics and wall-clock step rate. A
field's absence is represented by absence plus a reason in metadata, not by a
zero array or magic number. NaN/Inf causes a failed-state record and non-zero
run exit unless an explicitly named diagnostic mode says otherwise.

## Tracer metrics and trajectories

Trajectory tables key on `particle_id` and include time, position, velocity,
frame and selection-protocol ID. Representative Fig.11/12 selections must
store their initial-selection query; they cannot be chosen after observing the
outcome.

Output `R_count` and `R_volume` as distinct proxy metrics. A paper-area metric
must additionally record section plane/thickness, projection, ROI, AS/RS rule,
weights and pixel/geometry resolution. If the denominator is zero, emit
`value=null`, numerator and denominator counts, and `status=UNDEFINED`; never
add an unreported epsilon.

## Current verified fixture

`outputs/geometry_fixture/AFSD_GEOMETRY_FIXTURE.vtk` is a 118-point synthetic
ASCII legacy VTK data set with all five mandatory identity fields. Its manifest
is `GEOMETRY_FIXTURE`, states `physics_claims_allowed=false`, and records the VTK
SHA-256. `tests/test_geometry_fixture.py` reads the file back and verifies exact
persistent IDs, both tracer labels, host material IDs and initial coordinates.
It does not exercise a solver, neighborhood, stress update, contact, heat
transfer or integration path.

## Upstream PySPH evidence outputs

`outputs/upstream/pysph/` contains NPZ/info/log/profile files created by the
pinned upstream smoke runs. `logs/upstream/pysph/validate_outputs.py` read
their final states back and verified time and particle-array counts. These
files are retained as baseline evidence only: they do not yet contain the
complete immutable identity set or a project run manifest and therefore do
not conform to this schema. The adapter in
`src/afsd_backend/pysph/output.py` performs explicit conversion for new project
runs, not retroactive conversion of these upstream baseline archives.

## Implemented HOST_FKS_DEBUG profile (2026-09-24)

Canonical evidence: `outputs/debug/20260924_physics_verified/`.
Each case has `manifest.json`, `resolved_config.json`, `consumed_parameters.json`,
`source_hashes.json`, `series.json`, `DEBUG_state_NNNNNN.vtk`, `checkpoint.npz`,
`run.log`, `diagnostics.jsonl`, and `performance.json`. There is no synthesized
trajectory field or frontend animation. VTK positions are solver positions.

`series.json` holds step/time, available field names, each VTK SHA256, rigid-tool
translation in metres and rotation in radians about +z. Tool translation is not
applied to free particle coordinates. Rate fields represent the last midpoint
evaluation; `rate_evaluation_time_s` gives that time. The plastic strain rate is
the accepted increment divided by dt, associated with that interval. At t=0,
unevaluated rates are **absent**, with an explicit reason; zero initial stress,
zero plastic strain, imposed temperature and initial velocity are actual declared
initial conditions. The unused upstream checkpoint scratch arrays `e/ae/ax/ay/az`
are not canonical solved-energy/acceleration fields.

Stored field mappings:

| Logical field | VTK representation / meaning |
|---|---|
| position_m | POINTS, laboratory SI coordinates |
| velocity_m_per_s / initial_position_m | named vector[3] |
| stress_Pa | six named scalars stress_{xx,yy,zz,xy,yz,xz}_Pa; Cauchy stress, tension positive |
| strain_rate_per_s | six named scalars strain_rate_{xx,yy,zz,xy,yz,xz}_per_s |
| pressure_Pa | bulk EOS, compression positive; trace(sigma)=-3p |
| equivalent_stress_Pa | sqrt(3/2 S:S), including all three shear components |
| equivalent_deviatoric_strain_rate_per_s | sqrt(2/3 Ddev:Ddev), raw rate before documented FKS floor |
| equivalent_plastic_strain / equivalent_plastic_strain_rate_per_s | accumulated J2 strain / accepted delta per dt |
| internal_energy_J_per_kg | derived Cp*(T-Troom)+linear-EOS volume energy+S:S/(4G*rho); no fictitious heat solution |
| mass_kg, density_kg_per_m3, temperature_K, smoothing_length_m | named SI scalars |
| material_id | 1=AA6061-T6 for every deformable particle |
| component_id | 0=deformable substrate, 1=deformable feedstock; rigid bodies described separately by config/pose |
| tracer_type | 0=untagged, 1=initial exterior proxy, 2=initial centre proxy, only initial 3–4mm band |
| wholebar_tag | same tag codes but separate whole-bar diagnostic selection |
| actuator | 1=top-cap imposed x/z DOFs, 0=free material; not a tracer |

Persistent IDs are positive integer values; this Windows PySPH LongArray is
**int32**, not uint64. The implemented adapter checks assignment for overflow and
refuses reuse; current capacity is 2^31-1 IDs. A future wider-ID migration must be
explicit. All properties, including stride-9 stress/rate tensors, follow PySPH
particle extraction/reordering and are saved in the checkpoint. The material's
initial coordinates and passive labels never get recomputed from current position.

Checkpoint solver data includes physical t, count, next_id, cumulative heat/work
ledger, initial mass/energy, last evaluation metadata, rate-floor counter,
configuration fingerprint and source hashes. `feed_crossed` is a persistent
first-crossing flag. Actual infinite inlet creation is NOT implemented. Source or
configuration mismatch blocks restart. Upstream NPZ loading uses pickle: **never
load untrusted checkpoints**. On a failed step the last valid saved checkpoint
is retained; manifest/checkpoint step can be earlier than the attempted step.

Energy ledger signs (all joules unless labelled power in watts):

- wall_work_J and actuator_work_J: positive external mechanical work into material/contact.
- wall_heat_out_J: positive heat into rigid reservoirs, including their friction share;
  negative means reservoir supplies heat.
- plastic_heat_J, friction_heat_J, contact_heat_J: separate source contributions;
  contact_heat_J may be negative. They are not added again to total energy.
- total_J = kinetic + deviatoric elastic + EOS volumetric + thermal increment from
  prescribed initial T + contact penalty potential. Residual = total -
  (initial total + external work - outgoing wall heat).
- numerical_return_dissipation_J: reported backward-return loss, not silently
  converted to heat or used to cancel residual.
- external_mass_kg/external_feed_enthalpy_J: exactly zero for the finite closed
  particle collection, not missing values filled with zero. feed_plane_crossed_mass_kg
  is internal first-passage mass across the deposition-height plane, not external mass influx.

Diagnostics also include density deviation, penetration, extrema, momentum,
pair count, net conductive power, NaN/Inf alarms and every accepted dt. Thresholds
and normalization are in PARAMETER_SOURCES D10. Thermal boundary reservoirs and
the preloaded contact potential are included; no residual is hidden in the large
absolute initial thermal energy. Returned energies are small-strain/hypoelastic
DEBUG measures; finite-deformation thermodynamic consistency is not yet validated.

`scripts/validate_physics_outputs.py` independently checks manifest/artifact hashes,
all frame identities/finite values, stress trace/equivalent invariants, strict time
ordering and final VTK/checkpoint equality. Machine schemas remain version0.1.0,
with DEBUG/APPROXIMATE config classification and SI joule unit added; the host
series uses explicit named components, not JSON particle-record rows.

## 2026-09-30 coupling diagnostic extension (no solved field invention)

See COUPLING_FAILURE_REVIEW.md and COUPLING_DIAGNOSTIC_CONTRACT.md.
The legacy VTK field internal_energy_J_per_kg is a DERIVED_DEBUG_DIAGNOSTIC,
not a thermodynamically validated finite-strain internal-energy state.
The closed elastic cycle demonstrates why Es=S:S*V/(4G) is not a valid global
free energy for this updater. pa.e remains NOT_SOLVED scratch; no residual heat
or new solved-energy array is introduced. Keep this limitation in downstream plots.

Phase1 provenance.json stores source hashes, source/config fingerprints and the
resolved SI configuration. contacts.jsonl records ALL active faces per kick;
last_two_kicks.json retains a preceding full accepted kick and the rejected one.
Identity: step/stage, persistent particle_id, contact_id=ID:FEATURE, feature/face.
Geometry: midpoint centre, projected contact point, force application point,
normal, gap/penetration [m], before/mid/after activation and feature lists.
After geometry is explicitly an end-drift CANDIDATE prediction, never accepted.
Velocities [m/s]: before/after virtual internal/contact/constraint impulse sums;
the physical upstream integrator performs one combined velocity kick.
Wall translation [m/s], omega [rad/s], actual collocated wall/material velocity,
and reconstructed geometric-contact velocities are separate fields.
The latter are DIAGNOSTIC_RECONSTRUCTION_ONLY, not additional solved DOFs.
Impulses [N s]: internal, normal, friction, every other face, actuator reaction.
Work [J]: raw friction dissipation Q=-Jt dot relative_midpoint, lab friction/normal
work, internal/constraint work, delta kinetic energy and local kinetic identity.
local_energy_residual is ONLY that kinetic identity; full local thermal/free-energy
residual is unavailable before the rejected candidate's thermal evaluation.
No absent thermal candidate is filled as a solved zero or copied final solution.

ACCEPTED_DEBUG checkpoint: last_accepted.npz in diagnostic capture.
REJECTED_partial_integrator_state.npz: actual partially advanced arrays, NOT an end state.
REJECTED_candidate_fields.npz: explicitly predicted fields, temperature_not_advanced.
REJECTED_manifest.json: restart_allowed=false; these files are debug only.
Project restore rejects classification=REJECTED or restart_allowed=false.
Generic upstream pickle loading is not an authorization to use a rejected state.
Existing source/config/variant mismatch also blocks ordinary application restart;
old diagnostic checkpoints need an explicit forensic replay, not silent migration.

event_aware_v1 is recorded in solver state and manifest provenance; legacy remains
the default. Each accepted local interval stores dt, velocities, normal/tangent
impulses, constrained reaction, actual relative midpoint, raw heat per face,
wall/other/constraint work and kinetic residual, plus shared contact mobility.
The complete local trace is attached to rejected candidate diagnostics when available.
If the callback aborts partway, last_completed_contact_callback_trace may be stale;
failed_local_context identifies the actual failing particle/time/load/features
and completed local solves in that REJECTED candidate. Never merge them as one state.
New iteration-only diagnostic context includes projected-mobility cone data and
last four iterations, stored by test_new_candidate_red for the current blocker.

For actual FAILED runs, series/checkpoint.npz refer only to accepted output frames;
REJECTED_candidate.npz and REJECTED_contact_debug.json are separate.
Checkpoints are written to *.pending.npz and atomically replaced after complete dump;
pending files after a hard stop are not accepted output/checkpoint endpoints.
HARD_TIMEOUT.json (if produced) means INCOMPLETE with in-progress state rejected.
Per-run source_hashes.json remains that run's source revision, not automatically
updated by later diagnostic-only patches. New freeze integrity is separate from
scientific acceptance and never grants production permission.

## Contact454 numerical-coordinate extension (2026-09-30)

See CONTACT454_REPAIR_REPORT.md. Local kick traces now include tangent_bases
with shape [contact, Cartesian_xyz, tangent_component_2], dimensionless. These
are orthonormal plane coordinates, NOT a 2D mechanical/neighbor model.
The particle velocity, force, stress and neighborhood remain genuinely 3D.
Each interval adds tangent_coordinates_midpoint [contact,2] in m/s; the existing
tangential_relative_midpoint [contact,3] remains available. Heat is the actual
minus impulse dot relative-midpoint velocity in this orthonormal representation.
Failed cone context stores base_tangent [contact,2], tangent_bases, projected
mobility, coefficient, limits and the last four force iterates in tangent N.
Old failed-cone records with 3D projectors/base are immutable historical schema;
do not silently interpret their dimensions as the new schema.

ZERO_SLIP_WITHIN_ROUNDOFF events include time_s, contacts, residual_slip_speed
[m/s], zero_tolerance [m/s] and following=STICK_FEASIBILITY; VELOCITY_NOT_OVERWRITTEN.
They record resolution at the existing floating-point threshold; no velocity is
snapped to zero, and raw signed heat remains stored. A PASS artifact/readback
never overrides the run status FAILED or authorizes use of a rejected state.
The latest real run has 26 accepted steps/27 frames; candidate27 is REJECTED.

## Method-class numerical diagnostics (2026-09-30)

CONTACT_CLASS_CONTRACT/METHOD_CLASS_REVIEW define the newer numerical revision.
No new real AFSD series was generated in this stage; the previous26-step output
and REJECTED439 remain historical artifacts, never relabelled as successful runs.
Local kick velocities/impulses remain in LAB. Each interval additionally records
reference_frame_velocity [m/s] (a constant numerical velocity translation for
that kick, not a rotating observation frame or prescribed material trajectory),
relative_speed_scale_m_s and friction_rounding_bound_J. Raw signed friction heat
is never replaced by this bound. Cartesian/tangent midpoint representations
remain power-conjugate. Roundoff-level negatives are retained in raw output.

Cone failure schema now records scaled_residual, signed_gap, target, Newton rank/
line-search history and final tangent force; old change_N/threshold_N fields
belong only to the historical projected iteration. A KKT gap is an optimization
certificate, NOT a heat or energy field. base_tangent is a velocity for midpoint
solves and an acceleration for static-feasibility solves; interpret its units
with the caller (coefficient time vs dimensionless), never mix the two contexts.

Class-matrix JSON stores each case's independent reference certificate, measured
stationarity/gap/disk feasibility/objective and scaled endpoint error. Values are
dimensionless after recorded test construction/scales. Multiple valid face-force
splits are allowed; this does not authorize changing friction capacity or heating.
The extreme normal-boost test certifies tangent reversal/heat only, not accurate
subtraction of arbitrarily huge laboratory energy totals.
# E3-2局部接触数值修复补充（2026-10-01）

contact_cone.solve接口/返回三维力未变；失败cone上下文可增加
disk_projected_steps：无量纲归一化QP中每个圆盘的正投影步长，
不是物理时间步。旧上下文没有该字段时保持缺失，不填零。
contact variant仍为event_aware_v1，精确实现由source_hashes区分；
旧E3输出与REJECTED状态不被新局部成功改写，也不可跳过source检查恢复。
logs/contact476_20261001/final_local_artifacts为LOCAL_BENCHMARK，
不是新的AFSD时序；其中reference force来自冻结80位独立参考，
完整kick报告只验证冻结局部区间。新真实2us为NOT_RUN。
结果说明见docs/CONTACT476_SCALE_REPORT.md。
