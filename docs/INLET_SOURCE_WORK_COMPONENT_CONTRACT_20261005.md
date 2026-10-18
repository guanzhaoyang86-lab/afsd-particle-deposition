# 无应力出生端的入口源功组件合同

本轮只实现读取已有 PA 与耦合结果的来源计算和能量恒等检查；不求根、不做 NNPS、不自动签发任何 force/heat/energy PASS。旧求解器、入口 V2 和来源账不修改。

## 最小可执行基准及适用域

采用独立三维热短棒：底部固定、上部按明确世界 z 控制区给定送料，其他方向自由。补入粒子来自外部库存，初始 `S=0`、`rho=rho_reference`、`ep/plastic_work/return_dissipation=0`，温度与速度有逐项来源；两次补入之间真实求解力学与热学并重启。新生位置无重叠、三维邻域有秩、新生初始接触储能严格为零。外部出生端明确为零牵引；z 执行器作为单独的粒子控制力，其实际反力功只计一次。旧/新粒子之间的内部截面不另计“压力流功”。这属于离散库存出生模型，不能直接称作任意压力驱动、应力加载的开放流入边界。

后续真实基准必须与同域预置较长棒的参考比较。新生后完整图、修正梯度和切线重建；控制层释放事件、粒子密度与步长影响另作真实门禁。无联合基准时只能封存 V2 DIAGNOSTIC 链，不能升级旧状态为过程接受。

## 与现有模型一致的 SI 能量

现有材料储能诊断：

`U = Σ(m/rho) (S:S)/(4G) + Σ m (K/rho_ref)[ln(rho/rho_ref)+rho_ref/rho−1]`，单位 J。

它复用原储能定义，不是新的有限应变自由能声明。显热 `E_T=Σ m cp(T−Tref)`，cp 恒定，单位 J。出生弹性源从真实新生初值计算；受限无应力初值使其精确为零，而不是缺项补零。出生显热、质量、动量、动能与体积重新计算并核对旧 source_event。

`Tref` 可与 `parameters['room']` 不同，因为所有 B/P/A 库存和出生源使用同一显热零点，质量增量引起的参考项严格成对抵消；报告明确保存此温度及 cp。累计入口历史必须沿用同一参考和 cp（V1 registry 也强制该约定），禁止混合绝对热基准。

压力流功仅在外部出生边界显式零牵引时计算为零。执行器功来自同一实际耦合结果 `reaction = −force*(~free)` 与真实位移 `Σ reaction_z*displacement_z`。所有非入口固定轴必须零位移；否则本组件不支持额外运动边界功的分配，直接拒绝。该执行器功已包含于原 `prescribed_boundary_work_J`，合账时不得再加第二遍。

在出生前 B、出生后机械步开始 P、实际末态 A 上计算：

`R_open = W_actuator + W_tool + U_birth + Q_birth − Δ(U + E_T + C) − W_plastic_unthermalized − Q_wall_out`。

`C` 为实际接触法向与切向势能。要求 P 的旧粒子接触储能与 B 提供的真实账一致，新增粒子初值零；不掩盖出生引入的接触势能。`Q_wall_out = friction_wall_heat − contact_in`。

恒等式为 `R_open = R_mechanical_raw − R_thermal_external_raw`。后者按真实 `mcΔT` 与塑性热、材料分得摩擦热、接触输入计算；粒子间导热的全局舍入和仍保留为诊断。BE/Jaumann/返回算法的原始误差绝不扣除、转成热或裁零。参考温度改变所带来的出生显热与变质量库存偏移应抵消。

QS 没有求解惯性；入口动量/动能和实际动能变化单列观察，不进入上述 QS 总账，不补成温升。出生热库存可能远大于机械增量，因此不以总输入能量重定义或放宽原机械接受门槛。

读取原机械账时，`np.sum` 与本组件 `math.fsum` 的一致性舍入余量对储能差分采用被减两端的绝对储能量级；法向/切向势能差亦同。逐字段记录余量和来源，原报告 residual 与重算值分别保存。这只避免“大旧储能、小增量”的读取假 FAIL，不是新的开放系统接受容限，不改原 raw/前缀门槛。

## 接口与证据

`evaluate_source_work(parent_pa, prepared_pa, actual_pa, source_event, *, parameters, free, result, mechanical_energy, parent_contact_storage_J, external_port)` 是纯只读计算。结果必须与准备 PA 行顺序一致，实际 PA 可排序；父旧行除 actuator 外不能在出生时改变。`parent_contact_storage_J` 按父 PA 行排列，每行为法向/切向储能两列（J）。`external_port` 明确出生端零牵引、来源和“无额外内部截面功”约定。

返回值记录实际 PA/结果绑定、三个 `COMPUTED` 源项、完整开系统原始残差、恒等式差和限制，不自动设置过程 PASS。保存该原始报告后，`source_work_record(report, artifact_path)` 读回验证文件内容并计算 SHA，提供 V2 三项源功格式；外部 receipt verifier 仍须审核真实力热能门禁。

组件测试使用独立写出的质量/能量数值和制造数组，包括加热、弹性储能、负执行器功、非零原始缺陷、参考温度变换、QS 动能隔离、排序和缺项拒绝。它们不是 AFSD 或真正两次入口联合验证。
