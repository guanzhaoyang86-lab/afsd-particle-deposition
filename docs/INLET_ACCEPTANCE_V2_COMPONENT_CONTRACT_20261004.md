# 持续入口接受事务 V2：组件限定合同

本模块仅为 `COMPONENT_ONLY` 的事务基础层，不运行 AFSD、NNPS、材料、接触或热求解，不证明持续入口已经接入。旧 `inlet_protocol.py`、控制器和已有证据不修改。

开发隔离位置为 `scripts/inlet_acceptance_v2.py`；不在递归扫描的冻结 `src` 树新增模块。正式过程接入前不得把本组件加入旧模型来源集合。

## API 与单位

`initialize_acceptance_state(pa, registry, process_state, control_plane, provenance)` 只从尚无插入事件的 V1 registry 初始化 V2；`process_state` 提供按 PA 当前行排列的 `(N,8,3)` 接触位移历史（m）、`(N,3)` 布尔 `free`、`time_s`（s）及已有真实 `bore_observation`。旧初始资格和穿孔记录原样保留。

`propose_insertion(pa, state, incoming, specification, tool_pose)` 复用旧协议分配 ID 和来源账，生成脱离父态的 PA 和 PENDING 状态。事件时刻必须等于当前接受时刻。接触历史旧行按 ID 保留，新行零值明确是入口初始条件；新自由度默认自由，再应用入口控制区域。旧源账只能在本步骤累计一次。

`propose_mechanical_step(pa, state)` 在固定库存上准备下一次接受，允许按几何重新选择入口控制集合。`process_state_for(pa, state)` 将历史重排回给定 PA 的行顺序。

入口控制面是有来源的固定世界 `z` 平面，`component_id == 1 && z >= plane_z_m` 的材料只有 z 受控；平面以下释放原入口 z 约束，其余原固定轴不动。初始化要求现有 actuator/free 与该面一致，避免悄悄改约束。选择只发生在事务边界，边界间控制切换精度和反力功仍需后续真实事件验证；这里不规定内部轨迹。每次候选记录激活/释放 ID，反力功初始为 `NOT_COMPUTED`。

`seal_mechanical_acceptance(actual_pa, pending, contact_history, free, time_s, receipt, source_work, verify_receipt, pose_at_time)` 只接受外部合格凭据绑定的实际 PA；按 ID 密封更新的全部字段、接触历史和原 allocator/source 历史。不回写候选或父态。密封不计算材料/热/力，也不原子写文件；实际 runner 必须把 PA 和返回状态作为一份 checkpoint 原子发布。

## 封存条件

- receipt 必须绑定 PENDING state SHA、准备 PA SHA、实际 PA SHA、实际接触/free 摘要、准确起止时间及唯一 step ID。
- force、heat、energy、geometry_overlap、NNPS、visible_support、gradient_rank、contact 八项均需 `PASS` 和 SHA256 证据；外部验证回调必须确认这些现有凭据。组件假凭据仅使用 `MANUFACTURED_COMPONENT_FIXTURE_NOT_AFSD` scope；真实 scope 还需外部 `joint_inlet_benchmark` PASS，模块不自行签发资格。
- 独立真实入口基准可用 `DIAGNOSTIC_COUPLED_INLET_BENCHMARK_NOT_PROCESS_ACCEPTANCE` scope，仍需上述八项凭据、三项源功及外部验证，但不要求事先已有 joint 基准。因此可产生“插入—真实求值—诊断封存—重启—再插入”的基准链，再由独立对照审核其资格。诊断封存明确不可充当过程接受，不构成先验基准的循环依赖。
- 第一次封存锁定 receipt scope；制造、诊断和过程三类历史不能混用或升级。每个真实 PROCESS 步仍须引用已存在的独立 joint 基准凭据。状态永久有 `process_acceptance_by_this_module=False`，诊断绑定另记 `diagnostic_only=True`。
- `pressure_flow_work`、`elastic_convected_energy`、`mechanical_reaction_work` 都必须有 `COMPUTED`、J 单位、有符号有限值及证据 SHA；少项/NOT_COMPUTED 拒绝。压力流功与执行器反力分开且无重复计入的约定必须显式匹配。真实求解未计算时保留缺项，绝不默认补零。
- 质量、材料身份、出生/初始位置和被动标签不变；真实应力、温度、塑性及接触历史允许更新，但必须匹配凭据。原 IDs 不可遗失或替换，旧来源账不重复累计，不能通过重新初始化 allocator 续跑。
- 出生事件以出生位置相对实际肩面的上方关系确定新增料棒资格。出生本身从不新增 crossing；后续只对已经出生的真实两端状态调用原几何穿孔观察器。新增 V2 observer 合同明确区别初始资格与出生资格。
- 所有操作保持 `COMPONENT_ONLY` 分类。未经联合基准的制造测试不得写成实际入口或 AFSD 接受。

## 本轮验证范围

制造 PA：两次插入，中间修改位置、温度、应力和塑性模拟外部物理更新（不是 AFSD），接受后重启并排序。检查历史按 ID 保留、来源账不重复、资格/出生与 crossing 分离、入口释放和重选、凭据绑定、缺功项和粒子遗失/篡改拒绝。只运行新组件测试并记录命令、退出码与来源 SHA。

接口限制：本层保存逐事件来源账与三项有符号功凭据，但不计算完整过程能量残差或累计总过程能量账。入口压力流功与执行器反力的分割、弹性储能对流、其他入口牵引分量及 QS 中动能的解释仍由真实联合基准给出；三项字段齐备不等于总能量自动闭合。
