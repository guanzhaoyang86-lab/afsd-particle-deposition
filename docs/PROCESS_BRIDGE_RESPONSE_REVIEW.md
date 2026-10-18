# 真实342点完整响应：独立只读解释与下一求根任务

续记（2026-10-04）：第5节是当时的下一包设计；50ns的两种现有P作用检查现已
完成，见PROCESS_BRIDGE_ACTION_REVIEW及ACTION_CAUSAL_REVIEW。实际只用真实
C±/map±的correction差，没有采用原建议的零候选/代数伪map；500ns复算仍NOT_RUN。
原分支分析与NONSMOOTH结论不变，历史建议留存，不可误读为当前全部未运行。

2026-10-03。READ_ONLY；本审查没有新增evaluate、接触kick、轨迹或epsilon扫描。
原结果 `outputs/process_bridge_20261003/real342_initial/report.json` 的
**NONSMOOTH_NOT_CERTIFIED保持不变**，不得重标数学Jacobian认证或物理接受。

## 1. 实际发生了什么

四次真实ProcessThermalTrial调用均完成，四份map内部账PASS，首/末C全部数值
输出hash相同。原父态/来源完整性保持；342点同史切片有104接触、8多面粒子、
25活跃执行器粒子及25非零约束反力粒子；只有PID269本步新增塑性，三次均如此。
真实AFSD最新接受态仍850.393632910904ns，crossings=0；这里是650ns历史支路诊断。

本次独立读回脚本：

    .envs/pysph/python.exe -B logs/process_bridge_response_review_20261003/review_saved_response.py

运行耗时0.8295082s，exit0，new_physical_evaluations=0。它只使用gzip/JSON/NumPy
数组代数，不导入物理模块；逐个验证原report artifact及source摘要，结束再次核对。
读回结果：`logs/process_bridge_response_review_20261003/readback.json`，SHA256
`e223794d565146b3b0c1ec44f1205d59f5d71a31b9de151f3324b4a30935f664`。
同目录 `reconstructed_response.npz` 保存完整342行五场δC、δR、map差商以及
力/反力/塑性/热等数组；SHA在readback.full_arrays中。脚本拒绝覆盖已存在读回证据。
原始完整各face冲量/Q差商与所有字段差商，经独立减法与原directional_response
逐位核对一致；这只验证记录/运算，不验证导数极限。

## 2. 不同的是算法分支，不是已发现的接触拓扑/塑性切换

**CONFIRMED：** base/plus/minus的全部face ID集合及其顺序、free_axes、塑性active、
FKS rate-floor active、before/after normal active IDs均相同；所有粒子的停滑
事件数、顺序和涉及contact ID也相同。没有发现换面、约束区变更或塑性活跃集变更。

真实差异集中下表7个粒子，都是单面接触；feature定义来自
`rigid_contact_kinematics.FEATURE_NAMES`。子区间数不同导致原签名不同：

| PID / contact | feature | base / plus / minus 子区间数 | 额外变化 |
|---|---|---:|---|
| 202 / 202:2 | BORE_X_MINUS | 270 / 272 / 268 | 仅自适应划分 |
| 209 / 209:2 | BORE_X_MINUS | 472 / 480 / 517 | minus识别为roundoff-zero |
| 257 / 257:4 | BORE_Y_MINUS | 502 / 488 / 484 | 仅自适应划分 |
| 281 / 281:5 | BORE_Y_PLUS | 509 / 470 / 509 | plus识别为exact-zero |
| 286 / 286:4 | BORE_Y_MINUS | 341 / 341 / 278 | minus识别为exact-zero |
| 329 / 329:3 | BORE_X_PLUS | 517 / 517 / 474 | minus识别为exact-zero |
| 334 / 334:3 | BORE_X_PLUS | 268 / 268 / 270 | 仅自适应划分 |

总接受接触子区间7581/7538/7502，其中ADAPTIVE_MIDPOINT为7302/7258/7222，
EXACT_CONSTANT_DIRECTION为279/280/280。这里的exact-zero指旧kick的
ZERO_SLIP数值识别标签，不是独立连续解认证。

四个换标签事件仍落在**同一face**；其roundoff识别的完整切向速度norm及原界：

| 事件 | slip norm [m/s] | 原zero_tolerance [m/s] |
|---|---:|---:|
| PID209 minus | 1.44097920e-15 | 1.62391549e-15 |
| PID281 base | 1.53462653e-15 | 1.59057014e-15 |
| PID286 base | 3.45247086e-15 | 3.62973119e-15 |
| PID329 base | 1.47604211e-15 | 1.60232268e-15 |

完整原事件时刻/余量见readback.changed_traces.events。base相对kick的事件时刻为
202:12.4018936435ns，209:35.6491357198ns，257:37.2391972503ns，
281:34.3752453423ns，286:15.1191496326ns，329:35.8054106103ns，
334:15.3216529820ns；扰动时刻自然变化，无硬编码同步/速度投影。

**SUPPORTED而非证明：** 这些差异主要反映自适应积分/roundoff事件定位路线对
微扰的敏感，而不是新物理失效；同contact事件出现并不证明连续解对候选可微。
只读解释不能取消原签名、隐藏失败或改事件判别阈值。一般完整Jacobian仍未认证。

## 3. 五场完整差商与耦合响应

定义 `q_R=(R(C+)-R(C-))/(2e-7m/s)`，R=C-map；下表为逐分量绝对最大值。
实数数组和对应原分子均已保存，**不是拿期望答案生成的reference**。

| 五场 | max abs(q_R) | SI单位 | 极值PID/分量 | max abs(q_map) |
|---|---:|---|---|---:|
| position | 4.94634715e-8 | s | 1 / x | 1.82904906e-8 |
| endpoint velocity | 1.97850182 | 1 | 1 / x | 0.731591365 |
| density | 0.140560132 | kg s / m^4 | 198 | 0.140560132 |
| temperature | 0.00163026925 | K s / m | 269 | 0.00163026925 |
| event-mean velocity | 0.989253640 | 1 | 1 / x | 0.365795682 |

没有扰动candidate rho/T；它们非零的map响应来自真正耦合，不是预设温度场。
map速度/均速最大响应均在PID268的z分量。不能因这一个分量数值接近2倍，
就恢复全体event-mean等于端点平均的假设。

其它差商（每单位epsilon速度；极值保留符号）：

| 字段 | 实际值 | 位置/说明 |
|---|---:|---|
| internal force | -5.037447723 N/(m/s) | PID270,z，分量绝对最大 |
| constraint impulse | +2.518723862e-7 Ns/(m/s) | PID270,z，分量绝对最大 |
| total constraint impulse z | +1.038605807e-6 Ns/(m/s) | 25个执行器总反力冲量响应 |
| plastic work | +4.952412659e-7 J/(m/s) | 全部在PID269 |
| ep | -0.005812831869 /(m/s) | 全部在PID269；仍保持同一步正塑性 |
| return dissipation | -1.737416362e-6 J/(m/s) | 诊断，未当热源 |
| total raw friction | -1.427463182e-9 J/(m/s) | signed差商，不是负热 |
| total actuator work | -2.201844543e-9 J/(m/s) | signed外部功响应 |
| total wall work | +1.833578301e-10 J/(m/s) | signed外部功响应 |
| total heat input | +4.945275344e-7 J/(m/s) | 塑性源主导此方向 |
| local internal conduction | -2.082391445e-11 J/(m/s) | PID269；全域净和约舍入零 |

逐face法向/切向冲量与rawQ完整列表见readback.per_contact_directional_response
（104项）。例如rawQ逐粒子最大绝对响应在PID256；负的差商只表示Q随方向减小，
不表示Q为负，也不能拿abs处理。所有真实热事件原始正负值保留。

base旧C对新BE的五场scaled最大残差为178.32917365，来自PID269温度
8.91645868e-10K；其余力学分量仍满足旧单项尺度。这是旧显式接受候选对新BE
map的非零残差，不是新root接受。plus/minus最大scaled为197.61448/210.93456。

单epsilon左右差商并不完全相同：速度左右差的max约7.67345e-5，均速约3.83671e-5，
温度约1.13687e-6 K/(m/s)。这些包含真实非线性、分支离散及差分舍入，不能仅凭
该一组值分离原因，更不能据此给数学Jacobian误差百分比。未扫epsilon，未认证收敛阶。

## 4. 对真实342全局求根可用与不可用的结论

可用：已经得到同一parent完整五场响应以及力/塑性/反力/热分解，可以直接检验
现有数值搜索矩阵在这个真实方向上漏掉哪些响应。无需先证明处处光滑的Jacobian
才可尝试有明确适用域的近似搜索矩阵；每个候选最终仍由原actual residual和完整
物理检查裁决。数值proposal不得替代真实反力/热、不得变成粒子运动边界。

不可用：不能把该50ns/650ns父态的单方向数据直接当成500ns/850ns父态的局部
B/C或完整Jacobian；也不能据此证明500ns大步精度。原500ns失败记录保留：
`outputs/process_global_20261003/root_attempt/latest_trial.json.gz`，32评估后
velocity residual=7.255309545e-8m/s、mean residual=3.627653478e-8m/s，仍未接受。

实际API核查：

- `stage_structural_root.solve_stage`是原raw残差控制的候选迭代外层，不提供完整J。
- `ContinuityStagePreconditioner.apply(candidate,mapping)`可只做代数搜索作用；
  `StageEventPreconditioner.from_response_blocks/apply`也可复用，但要求完整且同ID序B/C。
- **目前没有已完成的真实342事件B/C**：
  `outputs/process_event_20261003/response_probe/manifest.json`是BLOCKED，
  response_progress仅PID1、5kick calls、ONE_SIDED_DIFFERENCES_INCONSISTENT；
  没有response_blocks.npz。不得把它或27点HHT结果描述为已可直接装到实际root。

## 5. 建议的唯一下一包：存档响应对实际搜索矩阵的离线作用检验

先做一个无新evaluate的有限任务，不再wall fixture、不先续50ns刷接受时间：

1. 固定使用本次存档的 `delta_C=C+-C-`、`delta_R=R+-R-`，不重取epsilon。
   读取650ns parent并按**原50ns、同物性/约束**构建已有
   ContinuityStagePreconditioner；构造矩阵不调用ProcessThermalTrial.evaluate/kick。
2. 对线性探针调用existing apply取得 `P delta_R`，检查
   `delta_C-P delta_R`，按五场、free/constraint、plastic、单/多面及零事件粒子拆分。
   可用五场全零**纯代数数组**作为apply候选、delta_R作为其mapping来避免大数
   加减损失；这些不是PA，不送入任何物理evaluate，绝不清零accepted状态。
   P本来就是近似逆，误差不要求为零，不预设用它通过任何物理门。
3. 只读原500ns最后FAILED candidate/evaluation及numerical_proposal_only，使用其
   **850ns父态、500ns dt**重建旧P，核原proposal可复算、ID/符号/连续性回代/约束
   行一致，记录原raw残差到proposal的作用。不要把50ns响应当500ns参考；若早期
   candidate日志缺完整map，不自行补造该步δR或调用evaluate填档。
4. 比较的是已实际测得的响应与明确近似P的作用；据分解结果只选一个搜索项修正
   （或判现有数据不支持因果定位）。StageEvent完整B/C缺失时标NOT_AVAILABLE，
   不调用旧有限差分kick补齐，不为了对照额外扩大作用矩阵。
5. 有了该报告再登记一次面向原500ns失败candidate的PRIVATE求根修正包：原父态、
   原dt、全部raw阈值/guard保持；任何新proposal须完整evaluate。其预算/候选次数
   先由实际成本预估，不在本审查执行；不要求先获得连续Jacobian全局认证。

本下一包矩阵构造/作用验证当前全部NOT_RUN；建议单次只读分析60s上限，超时保留
部分结果，不进入轨迹。停止条件是交付可重算的P作用误差分解与一个有证据的修正
范围，或明确无法定位；不是反复调参到响应“正确”。

真实本次成本：evaluate 3.43–3.91s，但含完整保存/独立检查为8.72–9.27s/次；
整进程内39.4136653s、外39.9551181s、工作集325861376B。执行前合同已显式把
预算从旧设计50/60s修订为110/120s（不改4次上限），supervision证据一致；
不能再把4.04s旧显式成本用于预测完整新BE根迭代。500ns步还需单独成本依据。

本次没有解除有限工具边界、50um位移域、长期稳定/准确性或过程输运阻塞。
本报告只把下一工作从局部wall基准转到真实342点搜索矩阵的可检查任务。
