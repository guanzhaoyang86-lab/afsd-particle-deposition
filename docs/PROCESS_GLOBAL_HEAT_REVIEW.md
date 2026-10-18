# 联合过程步热/功独立公式复核（2026-10-03）

## 结论与范围

**真实 500 ns 尝试没有接受新物理状态；最终候选为 REJECTED。**
本次独立公式复核支持“热与逐事件功账不是此次直接阻塞”，不支持“过程步通过”。
最后有效物理时间仍为 `8.50393632910904e-7 s`；候选所请求的
`1.3503936329109039e-6 s` 不能写成已达到终时、输运结果或重启初态。

复核者也是新热变体实现者；这里“独立”指参考公式和读回脚本不调用被测实现，
不是不同机构/人员的外部同行审查。此前“通过”声明没有作为数值参考答案。

证据：

- 最终运行：`outputs/process_global_20261003/root_attempt/manifest.json`
- 被读 journal：同目录 `latest_trial.json.gz`，而非物理末态 checkpoint
- 机器报告：`outputs/process_global_20261003/independent_heat_review.json`
- 可重跑独立读回：`scripts/assess_process_global_heat.py`
- 实际日志：`logs/process_global_20261003/thermal_variant_independent_review.log`

## 源码审查与时间层

逐段读取真实 `process_thermal_trial.py`、`implicit_thermal.py`、`coupled_boundary.py`、
`contact_kick.py`、`axial_actuator.py`、`rigid_contact_kinematics.py` 及真实 driver。
没有修改运行时已哈希的任一文件，没有再跑力学/热模型。

1. 材料/压力/内部力仍从同一 accepted 历史与共同中点几何、密度、温度计算。
   每次候选不累积旧迭代的塑性或温升；独立 `vbar` 仍进入梯度和事件平均位移。
2. 塑性功为 `(m/rhomid)*q_new*(ep_new-ep_old)`，单位 J；只有 chi 份额入热。
   `(1-chi)` 独立记账；材料返回算法耗散与机械残差没有补入热源。
3. 每个实际 contact 子区间先接受原事件冲量与原始带符号 Q；新证书又逐段验证
   `Q=-Jt·u_t_mid`、真实 frame/相对速度尺度及原 256 eps 舍入预算，之后才求和。
   没有用跨事件正负抵消赋予资格，没有 abs(Q) 或 max(Q,0)。
4. 工具速度与接触力都作用于粒子中心；投影 face 只是几何诊断，没有混用投影点速度来算功。
   底座零速度、工具平移加绕 z 旋转，可直接由 xmid 与 time_mid 独立复算。
5. H 为 W/K，B 为 W，均由真实活跃 face 和 base/tool 温库组成。
   `dt*(B-H*T1)` 为 J；热对用相同末温 T1。原 Tmid 通量仅保留为明确 diagnostic。
   BE 热图使用未修正的径向核梯度，没有将机械 B 修正混入导热。
6. 五场原始 residual 与旧机械功、净力/力矩、候选动量、执行器、事件均速等原门槛保留。
   仅不再使用原显式导热稳定 dt 限制；former bound 不是大机械步准确性证书。

本次没有发现需要改写已冻结真实试跑源码的热接线错误。
但新接受函数不是独立认证机构：本报告仍直接重算 journal，未调用其返回 PASS。

## 实际独立复算

脚本只导入 NumPy、标准库及 PySPH 的文件读取设施，不导入 AFSD 的 acceptance、trial、
contact、heat、material、kernel、NNPS 函数。检查运行 source/artifact/parent SHA256，
读取完再复查所有被监视文件未变。用核径向多项式独立重建梯度/热导，
逐事件重算冲量、功和约束，再用节点热容量及源项分别核对 mapped T1 与 candidate T1。
全部容限取既有合同，未依据结果放宽；独立切向投影的 Q 比较额外显式包含同源舍入预算。

| 检查 | 状态 | 真实证据 |
|---|---|---|
| 342点三维库存/ID、8131无向热对 | PASS，候选范围 | 三维秩与父库存，原径向核/对热导重算 |
| 104活跃face、14008接触子区间 | PASS，局部接触范围 | 每区间 Q、normal impulse、相邻速度连续及覆盖完整 |
| 摩擦功/热 | PASS，局部接触范围 | 总原始 Q=`2.2776385892915923e-9 J`；最负单项=`-2.1466732775777463e-24 J`，逐项均在其原舍入预算内且未裁掉 |
| 独立 Q 与记录差 | PASS | 最大 `2.490358642262329e-25 J` |
| 子区间动量 | PASS | 最大 `2.6469779601696886e-23 N s` |
| 接触功恒等式 | PASS | 最大 `1.6543612251060553e-23 J` |
| 刚壁功共轭速度 | PASS | 独立工具/底座速度与记录最大差0 |
| 映射速度整体动量 | PASS，map范围 | 最大 `2.117582368135751e-22 N s` |
| 候选速度整体动量 | **FAIL，不能接受候选** | 最大 `2.4486669854926623e-14 N s`，未达原 `m*1e-9+1e-17` 门槛 |
| 候选与真实事件均速 | **FAIL，不能接受候选** | 最大 `3.6276534772644586e-8 m/s`，原门槛 `1e-9 m/s` |
| 原五场总 residual | **FAIL** | 最大归一化值 `72.55309545164955`；v差 `7.255309545164956e-8 m/s` |
| BE mapped节点热账 | PASS，热子问题范围 | 最大节点残差 `1.7183926145792952e-17 J`，总残差 `1.5325038430944512e-16 J` |
| BE candidate节点热账 | PASS，热子问题范围 | 最大节点残差 `3.2236469862495777e-16 J`，总残差 `-3.0927948404413793e-15 J` |
| 塑性功/分配 | PASS，仅该候选账 | 该候选增量塑性功为0；不能据此验证持续塑性 |
| 温域 | PASS，候选范围 | `723.1499783215213–723.1500768163503 K` |
| 新物理接受态、真实输运、大步准确性 | **NOT_RUN / 未达到** | 新接受步0，不得把局部热账PASS当作全局推进 |

mapped 向材料的壁导热为 `-1.1000005585079798e-7 J`，材料热增量为
`-1.0886123640290182e-7 J`，差额来自材料份额的真实摩擦热，不是 residual heat。
BE 热对净和仅约 `6.56e-27 J`；内部导热没有增加全局热量。

## 仍然存在的风险

- 机械整体候选未收敛：局部 kick 都守恒不代表候选整体守恒。这正是拒绝而不提交的理由。
- 原 Jaumann/Cauchy 本构闭循环自由能问题未修复。Es、体积能及法向罚势的总账依然只是
  diagnostic；没有通过“总自由能闭合”的认证，也没有把该残差转成热。
- 当前所有热导、接触面与法向力在区间中点冻结；BE 只是这一固定系数热图的稳定方法，
  不是移动接触棱边/真实供料/长时传输/大时间步精度的证明。
- 实际有限特征 guard 和50μm位移域依然存在。本次没有改成完整有限工具几何，
  更没有因热通过而授权跨越孔口。
- 该32次预算耗尽的失败保留。不得放大容限、删失败或把1.3504μs候选作为下一段初态。

## 命令与成本

    .envs/pysph/python.exe -B scripts/assess_process_global_heat.py --case outputs/process_global_20261003/root_attempt --output outputs/process_global_20261003/independent_heat_review.json

实际由 `subprocess.run(..., timeout=60)` 包裹，首跑退出0；内1.311574 s，外1.644368 s。
状态固定为 `PASS_REJECTED_CANDIDATE_LOCAL_HEAT_REVIEW_NOT_ACCEPTANCE`、accepted=false、
production_authorized=false。重跑请换新输出文件，脚本拒绝覆盖已有报告。
