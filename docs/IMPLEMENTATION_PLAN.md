# 独立审查后优先级（2026-09-24）

## 当前：材料率相关搜索响应已修；转向剩余停滑块（2026-10-04）

读MATERIAL_RESPONSE_REPORT、INTEGRATION_DESIGN、OUTPUT_REVIEW。新增material_response
及material_stage_response，原物理代码/材料参数不变。独立11材料+12Schur回归首次PASS；
不是23个AFSD过程。新数值P保留非对称性、完整约束列、原后向误差门，不pinv/缩放。
同存档342点/650ns父/50ns方向作用外2.370s，速度inf由1.43e-7降至6.16e-9m/s，
L2为旧连续性P的.09137；温度交叉响应未改善。新候选均未evaluate/接受/重启。

不再在已补的材料项上循环。原500ns最后候选全部弹性、停滑点占均速残差平方99.54%，
不能把这里的塑性修复宣称为其根因或修复。下一针对剩余多面停滑的具名近似数值搜索
响应及独立反向/约束基准；不要求数值搜索矩阵先取得整个连续Jacobian光滑认证，
但每个真正候选仍须原实际map/功热/物理门裁决，不静默覆盖自由速度。
现有root刷新wrapper仅设计未实现，journal构造实际跑过，完整root仍NOT_RUN。

用户允许退而求其次：优先小型、真实受力、有明确局限的工程近似过程，不强求论文
定量三工况；不靠动画/预设流线、删除50um域guard、缩放或清零惯性冒充输运。
有限工具孔口与接触势/功仍须接入；最新真实850.393633ns/cross0不变。
自动goal当前usageLimited，进度落盘但未保证睡后自动运行。本节优先于下文历史计划。

## 当前：真实342点四次同史响应已取得（2026-10-03）

先读PROCESS_BRIDGE_REPORT、OUTPUT_REVIEW、RESPONSE_REVIEW；命令在
logs/process_bridge_20261003/COMMANDS.md。273回归PASS，唯一真实进程四次评估
外39.955s，全部独立map账通过、首尾完整输出逐位相同、父态未改变。
它包含正塑性/104接触/8多面/25执行器/BE热，不再仅凭27点局部块推测真实接口。
实际评估预算在运行前修订为110/120s；旧下节50/60s为历史设计，不能混用。

最终NONSMOOTH_NOT_CERTIFIED：7个单面粒子的自适应区间/停滑数值识别变化；
face、塑性active、约束集合不变。原始负Q事件保留并逐项核原舍入界，无裁剪。
旧显式接受C对BE的温度根残差仍178.329倍原门，不冒充新接受态。
真实AFSD仍850.393633ns/cross0；四次诊断没有推进物理时间，目标未完成。

当前工作直接转向真实五场求根：复用stage_continuity_preconditioner.apply，
用保存的C+/map+、C-/map-只读计算delta_C-P*delta_R，拆分自由/约束、塑性、
单/多面及零事件粒子，按实际证据选最小搜索项补丁，不再扩大wall fixture。
这项检查10月4日已完成（外2.639s），原工具序列化FAIL保留；连续性P的速度
L2余项比.06805、inf比.36685，最大余项在无直接接触粒子。独立全账户读回亦已
通过（外10.823s），空数组初次工具FAIL保留；无新实际evaluate。
下一补丁依据ACTION_CAUSAL_REVIEW对材料/几何/应力内力的存档分解，不由接触函数名猜测。
该分解已定位268,z的最大余项紧邻唯一塑性269；这一方向内力差主要为塑性返回，
几何/客观旋转项小。首选下一最小包：已有material.update的同支FKS rate-aware
塑性响应组件与独立材料点基准，再明确接入数值搜索；实际本构/门限不动。
必须涵盖塑性硬化/率/温度隐式导数、非比例增量，rate-floor和屈服边界不假称光滑。
尚未实现该补丁、尚未重跑500ns；不得把存档分解当作修复已通过。
完整342事件B/C矩阵目前NOT_AVAILABLE，不从27点HHT移植。不能把50ns单方向
直接当500ns导数；旧850ns父/500ns失败候选及raw门保持，后续独立注册有界修复。

过程终点仍须跨过有限孔口并真实输运：当前50um有效位移域小于250um肩距，
不能删guard续微步。有限工具组件尚未进入实际势/功与热接触；有限库存首穿越
可以先于持续插粒，但不能先做动画替代。整体自由能/同终时敏感性/成本风险未解除。
本节取代下文已过期的“342接口NOT_RUN”，不改变旧FAIL及原冻结证据。

## 当前：原trial005数值根修复，转向真实342点接口（2026-10-03）

先读PROCESS_RESPONSE_REPORT、LIMITS、OUTPUT_REVIEW、END_TO_END_REVIEW和NEXT_WORK_PACKAGE。
245组件回归与130独立数学检查实际PASS。新增contact_path_response、
contact_response_preconditioner、hht_response_stage；只改具名变体搜索方向，
实际evaluate AST不变。旧trial005同父态/dt/参数/门限，唯一PRIVATE诊断4评估收敛，
速度残差6.75346e-14m/s；旧18评估FAIL保留。整次外3.559s。
但跨零延拓仍不被动、不可发布/重启；连续跨零Jacobian不存在的独立反例保留。
4516旧项除三索引外不改。真实AFSD仍850.393633ns/cross0，局部接受仍6us。

不再把27点scout当作下一默认动作。下一单一包是实际342点五场响应接口：
复用ProcessThermalTrial.evaluate和公开非零塑性650.393633ns→700.393633ns历史切片；
最多基点/正向/反向/基点四次非提交评估，独立对账真实多面/轴向约束/塑性/BE热。
方向仅按几何与自由DOF构造，固定一档epsilon，不调到绿；所有评估从同一材料史开始。
该测量目前NOT_RUN；事前30s预测须第一调用后更新，内部50/外部60s。
保留latest850ns，不把历史切片当新进度；完整字段和预算见NEXT_WORK_PACKAGE。

工程路径仍需：真实残差的可负担过程积分、有限孔口feature功/势连续、多面执行器、
同史正塑性，以及有限库存首次中心穿越。现50um总位移域不允许250um肩距穿越；
不能删guard后冒充有效域扩展。持续插粒可晚于有限库存输运。
完整成本/同终时误差未证明可行；不生产三工况、不GPU、不静默删惯性或更改接触势。
本节取代下文“再做trial005”的历史下一步，但不删除其FAIL与原证据。

## 当前：独立反向参考到6.5us，耦合有效到6us（2026-10-03）

先读PROCESS_STOP_REPORT、IMPLEMENTATION、REFERENCE_RESULTS、CAUSAL_REVIEW和OUTPUT_REVIEW。
233组件回归PASS（189旧+44新）；不代表所有旧科学FAIL消失。旧4365项及失败保留。
新增hht_zero_event/stop_contact/stop_stage/stop_flow及独立reference_hht_event。
解析反向、完整矢量判零、单侧方向与状态/功热绑定已验证；耦合停滑仍未通过。
独立reference_endpoint2通过原9项界到6.5us；后续两网格粘着事件不同终时，不能称收敛。
scout_initial仅一次外4.576s：新解5→6us并发布；第三探根区间因近共线域/全局根失败停止。
不能把私有6.30/6.28us候选作为末态或重启；真实AFSD仍850.393633ns、crossings0。

下一最小工作包：先对完整contact impulse correction的acceleration Jacobian建立
独立解析/有限差分验证，区分物理滑移与PRIVATE反向延拓符号，包含法向gap和全路径积分。
再决定是否一次saved trial5私有回归；不得只改endpoint force导数、放宽原界或投影速度。
路径：hht_stop_stage的preconditioner、hht_stop_contact的实际积分、hht_stop_flow事件根；
保留旧定义，记录新策略来源/hash与所有FAIL。具体证据/限制见PROCESS_STOP_CAUSAL_REVIEW。
耦合顺序事件、3D粘着、正塑性、多面执行器、有限开放工具入口及真实输运仍待完成。
本轮不再跑scout；下一次必须先有独立基准及<=600s完整预测，不开三工况生产/GPU。

## 当前：曲线接触冲量一致修复通过，局部到5us（2026-10-03）

读PROCESS_CURVED_REPORT、DERIVATION、OUTPUT_REVIEW、CAUSAL_REVIEW、REPLAY_REPORT、NEXT_WORK_PACKAGE。
新增hht_curved_contact/identity/coupled及driver，旧4221冻结项除三索引外保留；不是新后端/换本构。
实际路径Coulomb冲量修正进入同一联合根；配套算法功单列、不转热。新schema绑定原始积分历史。
189预检PASS。原500ns RED重现后，新变体12评估PASS；没有改dt/物性/容限。
固定S0现在13步到5us，外9.296s；6us候选因参考只到5.9375us拒绝，不可续算。
独立读回+64点功热重算通过；4→5us一次真实checkpoint重算40属性/完整历史/账目bitwise一致。
但5us速度RMS误差.01325m/s、应力.409MPa，真实参考5.97466us停滑尚未解决；不能仅说缺参考。
下一包独立定位停滑—静摩擦可行性—反向滑动参考及事件积分，另冻结时间准确度/成本判断。
S1/S2/20us收敛/1ms/生产NOT_RUN。局部1ms粗估769–1742s超600；不是AFSD成本预测。
真实AFSD仍850.393633ns、穿越0；正塑性/多面执行器/真实有限孔口入口/输运未接入。
PROCESS_CURVED_FREEZE只证完整性；以下为保留历史范围。

## 当前：慢载启动揭示接触时间表示缺口（2026-10-03）

先读PROCESS_SLOW_REPORT、EVENT_AUDIT、OUTPUT_REVIEW、NEXT_WORK_PACKAGE；命令见logs/process_slow_20261003/COMMANDS.md。
新慢载制造输入/独立参考/有界driver/预检已实现，168回归PASS。原4067冻结项除三个索引外不变。
独立参考共同前缀到5.9375us，真实全向停滑约5.97466us；内部力4.675N>摩擦容量2.739N，支持反向滑动而非粘着。
实际S0在2us后首个500ns预测值触发方向误差；不是该停滑。一次私有同dt根12评估收敛，
但最终方向误差.012526>.01，故仅延后预测检查不足。未放宽门/调dt，私有态不能续算。
下一包是与冲量/做功/热一致的曲线滑移接触时间表示，先独立基准；随后定位停滑并检验反向滑移。
S1/S2、20us同终时敏感性、1ms、实际慢载重启NOT_RUN；不堆旧尺度微步。
真实AFSD接受仍850.393633ns、穿越0；本轮9个局部接受文件含一个成本重复步，不是AFSD进度。
后续正塑性、多面执行器、真实方孔/入口、输运仍须接入。PROCESS_SLOW_FREEZE只证完整性。
下述各节保留历史，以本节的新证据和明确范围为准。

## 当前：HHT局部耦合贯通；同终时速度精度FAIL保留（2026-10-03）

先读PROCESS_HHT_REPORT、OUTPUT_REVIEW、NEXT_WORK_PACKAGE；命令见logs/process_hht_20261003/COMMANDS.md。
新增hht_linear/path/coupled及probe/check/restart脚本；原AFSD/物性/空间/材料路径不替换。
已修正加速度预条件中的密度—压力残差漏项，独立整块代数及原失败step6回归通过。
最终预检145PASS；8/16/32档均到2us，完整三档外6.314s，3D受力/导热/接触输出是真实求解。
细档8/9绝对精度PASS，速度2.0520e-4m/s>2e-4；9个细化比均PASS。整体仍FAIL，不调阈值/追加网格。
旧失败、初次第6步失败与每个接受/拒绝文件保留；独立56状态账目/源码读回通过。
严格半程重算脚本已实现、7回归PASS，但实际重算因精度门失败NOT_RUN。
新真实AFSD接受0；仍850.393633ns、穿越0。本次27粒子制造弹性块不是棒/基板沉积。

下一单一包是慢载、保惯性且包含完整启动瞬态的过程步精度/成本判断；先注册独立参考与预算。
本轮500kHz制造壁运动不代表300rpm过程步经济性；不再追这2.6%超界或堆同尺度微步。
不得省掉启动、清零旧速度或静默切准静态；模型变更需用户方向。
后续正塑性功、共享多面/执行器与停滑事件、有限工具/入口、真实输运仍未接入HHT。
PROCESS_HHT_FREEZE只证完整性，不是生产/论文复现许可。以下记录均为历史范围。

## 当前：事件导数有界实验已收束；改进过程时间积分，不再调h（2026-10-03）

先读PROCESS_EVENT_REPORT、NEXT_WORK_PACKAGE；命令在logs/process_event_20261003/COMMANDS.md。
实际新增event_response(B/C分离)、stage_event_preconditioner和有界driver；原物理trial/kick/材料/空间未改。
147项不同回归145 PASS/2 FAIL：连续ODE精度和真实PID1两次差分细化预算均失败；原RED与原日志保留。
原cloud在PID1停下，外1.103s；修订cloud/成本/500ns根均NOT_RUN。新接受0，真实末态仍850.393633ns、穿越0。
已证PID1早期侧差是平滑O(h)曲率，但另一轴在注册预算内未过；不增加减半次数或放宽门。
独立时序审查确认旧kick返回细解，但只控局部末速度，未控制全过程误差/实际均速积分。

下一ONE工作包为PROCESS_EVENT_NEXT_WORK_PACKAGE的具名保惯性时间积分与算法能量验收，
先导出并验证光滑能量身份及事件路径/功一致性，再小3D同终时精度与成本交集。
HHT候选是DESIGN/NOT_IMPLEMENTED，不可直接把旧kick插入端点Newmark公式。
若精度和成本无交集，报告该路线不可行；不再把若干500ns接受态当成输运进展。
有限孔口、入口接受事务、真实材料输运和重启/敏感性仍是后续未完成项。
PROCESS_EVENT_FREEZE只保文件完整性，不会把FAIL变PASS或授予生产/完整复现。
下面为保留历史，以本节与独立报告的明确范围为准。

## 当前：结构/连续性/事件接触/BE热已接入同史联合候选（2026-10-03）

先读PROCESS_GLOBAL_REPORT、ROOT_CONTRACT、STEP_AUDIT；命令在logs/process_global_20261003/COMMANDS.md。
新增ContinuityStagePreconditioner→stage_structural_root→ProcessThermalTrial/acceptance；
process_global_state保存具名新变体，原物性/接触kick/空间/材料/自由运动未改。
103项预检PASS；唯一500ns尝试32评估、外265.559s后预算FAIL。raw速度7.2553e-8m/s>1e-9；
不扩迭代或换dt重跑，REJECTED完整保留。新接受0步，真实终时仍850.393633ns。
BE已在真实候选连线，不再仅离线组件；但尚没有新BE耦合接受轨迹。

下一联合数值缺项是事件平均接触响应进入结构Schur（先独立反向/多面/执行器响应基准），
不是降低raw门或继续堆50ns。详细残差集中证据见本轮独立读回报告。
过程时间精度仍需要保留惯性的具名积分/算法耗散验证；500ns高频欠解析不能当过程许可。
有限孔口/肩缘事件功、真实穿越、入口接受事务与重算重启仍未解决，最终门不变。
以下保留之前组件/轨迹记录，其中NOT_INTEGRATED以本节明确更新范围为准。

## 当前：面向完成态的过程依赖组件；尚未联合集成（2026-10-03）

用户要求继续推进到项目结束，完成门见PROCESS_COMPLETION_CONTRACT；不能以测试数/短步数替代真实料流。
固定PySPH、旧动力学轨迹/FAIL不变；当前真实终时仍850.393633ns，未启动长跑或生产。
新组件：finite_tool_geometry（有限开放方孔/连续单实体法向势能）、inlet_protocol（持久ID/原子候选与显热源账）、
implicit_thermal（同径向热图的保守BE）、structural_mobility（M+dt²K/4受约束弹性结构响应）。
这些不是新后端；都显式NOT_INTEGRATED，不替换当前UpdatedTrial，不自动放开50um守卫。
入口：scripts/check_process_readiness.py --output 新目录；详细证据见PROCESS_COMPLETION_REPORT。

下一统一集成单元，而非继续堆50ns：

- 构造PROCESS_INERTIAL_EVENT_STEP_APPROXIMATE同史联合残差，复用现NNPS/B/材料/原子提交和输出。
  structural_mobility仅弹性冻结构形预条件；补实FKS率/几何/热耦合方向导数和事件contact响应，不能冒称完整一致切线。
- 隐式热已能独立取消显式稳定界，但旧实际界约0.292ms远大于50ns；它不是当前力学瓶颈的独立解药。
  需从实际每个contact face提供H和sum(H*Twall)，用同一T1热时间层记账，不能直接传旧wall power。
- 有限工具单实体距离罚力不同于逐face弹簧。尖棱/中轴多最近点活跃时仍拒绝，需联合接触事件/法向锥设计与独立反例。
  保留平肩/开放旋转方孔，不暗自圆角、封孔或交换AS/RS；未通过集成前不放开位移域。
- 入口PA候选还需要物理重叠/NNPS/接触检查、接受history授权的rebind+commit、质量/动量/焓含压力流功账和新输运算子。
  当前registry绑定全部状态，物理推进后必须合法重新绑定，禁止重建注册表造成ID复用。
- 做一个联合动态-接触-热同终时基准及实际成本，再决定唯一<=600s过程探针；各组件PASS不授权真实宏步。
  若改热态准静态起点，只能单独设计/初始化并说明省略启动瞬态，不能清零当前checkpoint速度后续接。

当前过程输运、完整焓/自由能、宏步准确性与联合重启仍未完成；项目不标complete，不授予论文完整复现。
以下是旧阶段记录，不覆盖最新组件范围或把旧微秒结果当过程完成。

## 当前：整体阻塞已跨过，4个相连真实接受步；过程输运仍缺（2026-10-03）

先读NONLINEAR_CONTINUATION_REPORT、AUDIT、CONTRACT及CORRECTION；命令见logs/nonlinear_continuation_20261003/COMMANDS.md。
同一PySPH pin、同一五场UpdatedTrial/物性/接触/热/接受容限；仅新增有界Anderson求根与提交/观测driver。
40不同短回归最终PASS，首错/旧FAIL保留；唯一真实推进61evaluate、4×50ns，外237.387s。
新分支从旧step62的650.393633ns推进至850.393633ns，每步NPZ严格重载后用于下一步。
独立读回43源/97产物、23411接触子区间PASS；自由棒最大净位移36.174nm，固定肩面穿越上下均0。
仅第一阶段PID269新增塑性，后3步无新增塑性；这是连续短步接受，不是持续料流/沉积或有效75%。

下一主要单元不再是局部接触微调，也不是堆相同50ns步：

- 过程时间步：保留updated_trial.py/stage_anderson.py基线，先为可用宏步定义独立惯性/事件/热/客观率准确性与成本。
  大步根收敛不等于准确；当前预压瞬态不可直接去惯性，旧制造Q7不可当实际过程后端；不得隐藏缩放。
- 过程物理域：原始方棒75点肩上/75点肩下，最近绝对间隔250um>冻结50um位移域，当前配置不可能中心穿越。
  先补geometry.py/rigid_contact_kinematics.py/coupled_boundary.py的有限孔口/肩缘进出与功测试，不能只放guard。
- 持续供料：particles.py/output.py新增ID/质量/动量/焓账与可重启入口；当前有限棒和observer不支持动态插粒。
  标签只观测，已有肩下材料不是新沉积；测量平面不为结果而移动。
- 重启与同终时对照：continuation_state.restore当前仅内部严格往返；安全外部入口还须检查正式accepted_history。
  重新求解的两轨迹重启等价、dt/空间敏感性均NOT_RUN。禁止比较不同终时宣称收敛。
- 成本正确文件process_feasibility_corrected.json：.117925s仅名义执行器距离时标，当前固定步代价粗外推约1614天，非可信长时预测。
  首版process_feasibility.json负成本逻辑FAIL/INVALID保留。无理由启动长跑/GPU或三工况；先解决路线与范围。

NONLINEAR_CONTINUATION_FREEZE仅证完整性；旧材料自由能/低率正则/热启动/预结合近似和FAIL继续有效。
以下均为保留历史，不得用旧“无新增步”或制造准静态PASS替代本轮真实状态。

## 当前：同史整体trial已接通，真实Picard全局步拒绝（2026-10-03）

先读UPDATED_TRIAL_REPORT/CONTRACT/INDEPENDENT_REVIEW；命令见logs/updated_trial_20261003/COMMANDS.md。
新增updated_trial.py复用当前3D邻域/伴随空间、原Jaumann/J2/FKS、事件多面/轴向约束和热分账。
五场x/v/rho/T/事件平均速度来自同一accepted史；trial不提交、不累计ep/热，不以首末均速代替事件积分。
35不同短回归最终PASS，两条fixture首FAIL保留。一次真实完整候选含证据3.81249s，保守预算213s。
唯一50ns整体尝试从旧塑化step62出发，32evaluate完成、外145.590s，Picard残差停滞后REJECTED；无新增有效步。
全局失败不是负摩擦耗散异常；不能标75%，旧有效轨迹仍342点/67步/700ns，非持续沉积。

下一主要单元是整体求解策略，不再反复补旧接触局部比值：

- 保持UpdatedTrial的五场同史map接口，改进独立整体求根模块；不要同时换材料/核/热/接触。
- 先注册有独立根的惯性、多约束、分段非光滑联合基准与成本，不用制造线性PASS代替真实root。
- 明确非线性尺度/步长与失败回退，不靠更长Picard或降低事前容限把本次REJECTED变接受。
- 复用最后实际evaluate的candidate保存约定；map端点不能未经feature检查直接提交。
  候选/接触map功差、法向时间势能残差及角动量额外项必须原样报，不补热/外矩。
- 真正多步推进可用后再处理持续入口/有限边缘、同终时敏感性和沉积输运；三工况仍不授权。
  自由能材料及现有边界近似另立方案，旧FAIL不删。新snapshot续算、过程尺度、论文定量验证NOT_RUN。

UPDATED_TRIAL_FREEZE只证完整性；以下保留历史计划和结论，不覆盖其失败。

## 当前：过程时间路线已判定，下一实现更新构形纯候选（2026-10-02）

读PROCESS_ROUTE_REPORT/DESIGN/CONTRACT；复跑见logs/process_route_20261002/COMMANDS.md。
本轮不重复微秒试跑；新process_residual.py实现保留真实质量/惯性的Rx/Rm/Rc及Rq诊断，未接生产。
新38不同回归PASS（23+15）；读取器首次FAIL和两次真实入口FAIL保留，最终全量只读67步/189382子区间完成，外13.117s。
新读取器共3修正超过自设2次预算，最后一次仅登记兼容例外已在PROCESS_ROUTE_REPAIR_NOTE明确，不改物理/容限。
实际Rq比例.149–.777，末步.348；动态残差max7.028e-15N，说明当前预压瞬态不能直接删除惯性接准静态。
不是从此否定准静态稳态AFSD；只是旧Q7冻结几何/制造材料不能当当前checkpoint续算器。
反向摩擦解析例证全步中点动量可闭合但热/位移错误；隐式仍须同一子区间事件与功。

下一实际代码单元（保持1–4整体，不再拆旧接触比值）：

- 新纯trial_increment：accepted PA + candidate x1/v1/rho1/T1 + dt/p/rigid_program；构造独立中点PA，
  当前NNPS/corrected B/L、旧材料从同一S0/ep0起算、同构形压力/伴随力、真实工具face/轴向约束。
- 将本轮Rx/Rm/Rc与Rrho连续方程一起输出NOT_ACCEPTED；试探不能累计材料史/热或写accepted checkpoint。
- 精确复用mechanics.py/material.py/rigid_contact_kinematics.py/axial_actuator.py/coupled_thermal.py。
  Q7仅借鉴有界迭代/纯trial/完整受约束切线框架，不继承其制造J2、固定X面、冻结D/B或常rho。
- 联合基准先验更新3D仿射/非零应力刚转/功与密度/反复trial不改接受态；方向导数说明几何项，不能称Q7局部切线完整。
- 此后才测trial/Jacobian成本并决定一次<=600s全局试步；全局接触未知力来自同一candidate，滑移/feature事件必须解析。
  热宏步若越行和界需已验证子循环/隐式热，不把原热函数拒绝当自动子步。
- 真实持续入口、有限边缘、材料自由能、预结合/热启动依据仍缺。三工况、宏步/沉积和当前敏感性NOT_RUN。
  新隐式残差不是新隐式求解器，不能更新总体进度为“已完成过程时间模拟”。

真实状态仍上一轮342点/67步/700ns；未删除旧FAIL，不授权生产/GPU/付费/长时。
PROCESS_ROUTE_FREEZE仅证完整性。以下是保留的历史阶段。

## 当前：四关统一真实三维热力联动已运行（2026-10-02）

用户要求整体推进1–4，以下历史“逐项下一关”不再是当前计划。
先读COUPLED_AFSD_REPORT、CODE_REVIEW、CONTRACT及SHORT_WINDOW；可复跑命令在logs/coupled_afsd_20261002/COMMANDS.md。
同一PySPH pin；新增IntegratedSolver适配器复用现有粒子/NNPS/材料/积分/输出，旧物理和FAIL不改。
开放旋转方孔+原始顶层z-only送料，全192点基板/150点棒可变形，刚工具/底座，真实导热/塑性/摩擦/接触热。
新83不同短回归PASS，旧选择41PASS及raw2FAIL保留；首次solver测试失败/reader舍入失败日志全部保留。
成本原2/4us超预算NOT_RUN；另预先登记同物理700ns短窗，67步实际外132.273s完成，不拆分绕600s。
全量独立68状态/189382接触子区间/8VTK/checkpoint读回PASS；step2和首次塑性step62各续2步逐字段零差。
2个自由棒点塑化，ep_max1.8228662e-5，只有2步新增塑性；基板塑性0，送料穿越0，肩面未接触。
已达四机制同例共同起作用，未达持续塑性输运/沉积。Es非认证自由能；R4.11617e-7J约外功7.26倍，风险不能隐藏。

下一有效单元是整体过程时间尺度路线，而非重复相同接触局部比例：

- 固定本轮为接触/热力联动回归；132s仅推进0.7us，同一显式路径直接秒级推进不现实。不要继续靠微秒反复重跑充进度。
- 评估现有准静态组件扩展更新几何、真实开放工具/轴向约束、材料和热源的工程近似方案；Q7当前不是可用AFSD长时后端。
  先给不忽略关键受力/热源的数学范围、独立功/热/材料点与小3D基准和成本；路线确认后再实现，不静默替换当前显式分支。
- 相关源路径：solver.py/integrator.py/coupled_solver.py时间接口；material.py独立自由能候选；
  geometry.py/rigid_contact_kinematics.py有限孔口/肩缘；particles.py/output.py持续入口ID/质量/焓账。
- 未决物理：有限接触棱边、预结合板棒、热启动/预压程序、FKS参数与低率floor、Jaumann自由能/分裂误差。
  工程近似可继续方法开发，不能把缺项改成已验证或用残差补热。
- NOT_RUN：本轮时间/空间敏感性、原2/4us目标、连续供料、肩面承载沉积、宏观温场/三工况/论文定量对照/GPU。
  不授权三工况生产或长时付费任务；单次仍<=600s，可视化只读冻结物理。

COUPLED_AFSD_FREEZE只证完整性，不证物理通过。以下保留此前阶段，不覆盖旧证据。

## 当前：局部轴向驱动—自由材料—底座—释放十项完成（2026-10-02）

优先读AXIAL_FEED_REPORT/CONTRACT及logs/axial_feed_20261002/COMMANDS.md。
新增轴向边界companion，仍同一PySPH pin/原corrected+adjoint、材料与contact_kick；不改旧生产默认。
当前64点全可变形，没有底层夹具；只对原始顶16点规定z速度，xy自由，400ns后全部xyz自由。
驱动速度2.12mm/s，行程仅0.848nm；系统封闭，没有连续物料入口，没有基板，非塑性沉积。
新40短回归PASS，旧选定31PASS；首次排序checkpoint ERROR保留，已固定初始账按ID求和顺序。
两档160/320步同至800ns，外8.254/14.456s；482帧独立读回PASS，真实释放checkpoint各续两步逐位一致。
最小邻域rank3、有xz/yz/内部自由响应、末态16点接触；ep全程0，thermal NOT_SOLVED。
161共同帧位置/速度/应力/密度screen内，但四个功/能量量超过1e-12J，整体SENSITIVE_WARNING。
粗/细ENERGY_WARNING为141/156帧；旧Jaumann闭循环及空间/上一法向阶段警告不改判。
本轮不因能量告警无限重复同一块，不加残差热、不缩放物理，也不声称dt收敛。

下一有效工作单元：轴向边界与开放旋转方孔耦合，先有限feature与功共轭，再一个局部组件。

- 必读路径：src/afsd_backend/pysph/geometry.py(surfaces)、rigid_contact_kinematics.py(contact_rows)、
  contact_kick.py(domain_kick/kick)、axial_actuator.py(control)；复用而非另建后端。
- 新独立接口基准必须覆盖：孔内无接触力、孔壁/角部多面、共转零滑、自由xy+受控z的反力功、
  有限面切换的gap/势能一致。采用投影点速度重构时要配同一映射的伴随力，不能只改壁速。
- 通过后在新命名companion组合现有边界；从两步成本开始，冻结实际几何/温度/速度/接触尺度，
  只对入口/工具施规定运动，不为自由材料规定轨迹。单次外部<=600s，不静默延长原预算。
- 接口失败只修最小适配路径；不借机替换本构、恢复raw空间路径、修改本轮警告或删除反例。
  工程近似允许带已知精度风险继续独立接口开发，不等于通过宏观塑性/热学验证。
- NOT_RUN：本轮第三dt/空间加密、旋转孔壁动力传递、热、连续供料、棒+基板塑性沉积、三工况。
  真实AFSD历史仍2us/191步无持续沉积；不授权生产。AXIAL_FEED_FREEZE只检查完整性。

以下阶段保留历史证据，不把其中旧“下一任务”覆盖本节。

## 当前：法向承载—卸载—开离十项完成，精度警告保留（2026-10-02）

优先读NORMAL_CONTACT_REPORT/CONTRACT。本轮新增命名companion，不改旧生产物理默认。
新43短回归PASS，旧选定26PASS；四个解析单点、两个3D组件和两个预抬升重启实际完成。
64点底16固定，其余48点xyz自由；Fn由gap求，正常弹簧能/壁功/材料功/时间残差分账。
5/2.5ns到共同600ns，外6.178/11.365s；362帧独立读回PASS，全16顶点约275ns开离。
121共同帧速度差.01648m/s、应力差194151Pa超过事前screen，SENSITIVE_WARNING保留。
其余筛查量及开离区间通过，不等于全轨迹连续体精度/收敛阶通过；没有再调参或跑第三dt。
ep实际0、thermal未求；旧Jaumann能量/空间FAIL和AFSD2us/191步历史不变。

下一工作单元转向局部送料接口，不继续打磨同一个卸载比值：

- 新命名局部执行器companion，复用particles、NNPS、HostLeapfrog、contact_kick。
  冻结入口层厚/目标轴向速度/固定物理几何，只有执行器层规定运动，自由材料不得覆盖速度。
  先用独立小基准验constraint impulse、执行器功、自由点传力、ID/标签不干预。
- 涉及旋转孔壁前，读geometry.py、rigid_contact_kinematics.py、contact_kick.domain_kick：
  面点速度重构必须配伴随力；有限feature切换的势能/间隙要成对。当前平面U诊断不能泛化。
- 不借机改material.py本构/模量或恢复raw空间路径；保持corrected与adjoint成对。
  新送料测试先预估预算+两步成本，单次<=600s，不把本轮screen警告改成精度已验证。
- 多面送料耦合、热、塑性沉积与空间细化仍NOT_RUN；本轮警告作为后续精度门槛保留。
  三工况生产不授权。本轮运行/失败日志及可重跑命令见logs/normal_contact_20261002/COMMANDS.md。

以下各“当前”标题是保留的历史阶段记录。

## 当前：旋转肩面—真实三维可变形小块五步完成（2026-10-02）

读ROTATING_DEFORMABLE_REPORT/CONTRACT。本轮冻结64点/2mm立方体制造问题，
新companion复用原3D NNPS、校正梯度/伴随内力、Jaumann/J2/FKS、contact kick及HostLeapfrog。
新26PASS(14短组件+12独立读回回归)，旧选定26PASS；所有旧FAIL保留。
两档40/80步同到200ns，外2.070/3.714s；独立读回PASS，122个VTK均核对。
41共同时刻比较WITHIN_SCREEN，最大速度差.000466523m/s、全应力差4530.068Pa。
实际内层xyz响应、xz/yz剪应力非零，ep=0；是三维弹性载荷传递，不是塑性沉积。
原材料能量/空间自由表面风险未消除；源代码与参数/几何相同，两档不声明收敛阶。

下一不继续重复固定z夹具：验证正常承载/开离后接局部送料。
- 先读contact.py WallContact、contact_kick.py domain_kick、geometry.py surfaces。
  当前Fn为规定摩擦容量且在顶层被z夹具抵消，不能冒充工具压入/基板预压。
- 冻结独立单点/小块压入—卸载基准：间隙、接触刚度单位、正常势能、壁功、
  动能、反弹/开离事件，严禁吸引力/负功裁剪，明确法向冻结力的时间误差。
- 用命名companion接入已有penalty或已证明的法向更新；不静默改旧生产默认。
  测后才解除本小块顶层z夹具，底部支承保持明确；不得给自由点规定压实轨迹。
- 接口输出contact feature/活动状态/gap/normal impulse/reaction/normal potential，
  正常开离回归通过再考虑局部轴向执行器。先测dt/两步成本，单次<=600s。
旋转孔壁面点速度重构仍需配转置力传递；热/送料/有限应变/FKS校准/空间细化未做。
三工况生产不授权；真实AFSD历史仍2us/191步、无持续沉积。

## 当前：旋转刚壁接口+自由质点动力组件两项完成（2026-10-02）

读ROTATING_CONTACT_REPORT/CONTRACT；以下小块状态为历史，均保留。
新增rigid_contact_kinematics.py，区分中心作用点/有限面投影/两处壁速，
显式诊断偏心功差、不作热；开放孔与角部双面、共转伪滑反例17PASS。
新增probe_rotating_shoulder_component.py复用旧HostLeapfrog/contact_kick，
自由xy响应，z规定Fn/反力记账；这只是平面4点质量/三维接口，不是SPH固体。
两档各一次40/80步到.0008s，外.617/.908s；独立RK4参考和逐子步读回PASS。
新45PASS（17接口+11短组件+17读回），旧contact12PASS；首次严格零断言FAIL日志保留。
不是完整AFSD、不是热/有限应变能量认证，旧材料/空间FAIL不变。

下一工作单元：旋转肩面与可变形SPH小块接入，不继续打磨点质量比例。
- 阅读src/afsd_backend/pysph/solver.py的evaluate与contact_kick.py的domain_kick，
  复用既有PySPH粒子/NNPS、mechanics/material/HostLeapfrog；不得用本轮无内力
  的点质量组件替代金属模型或把全块刚化。
- 新命名companion（建议scripts/probe_rotating_deformable_block.py），先冻结等温、
  小应变/短程物理几何、载荷/法向约束和入口范围；只给刚壁/明确执行器规定运动。
- 以rigid_contact_kinematics.contact_rows联接逐接触ID/feature/作用点/壁速，
  原kick保持共享mobility及真实reaction；逐步记录内力功/接触功/塑性功/存能告警。
- 初次运行前实测稳定dt、2步成本与内存，外限<=600s，单次有界试探。
  温度未解不填假字段。先验力/质量/功/NaN门槛不靠剪裁放宽。
- 旋转孔壁的真实面点重构若用于受力，需与速度重构转置分配配对；当前只是诊断，
  不借平肩特例宣称孔壁已验证。正常开离/轴向送料/热/FKS/有限应变继续列缺项。
三工况生产/空间细化/GPU本轮NOT_RUN；真实AFSD仍2us/191步无持续沉积。

## 当前：同物理4/8检查完成，推进旋转接触接口（2026-10-02）

读PROCESS_PREDICTOR_SENSITIVITY_REPORT/CONTRACT；以下Onset状态为历史。
冻结求解器/23项来源/物性/几何不变，唯一8增量88.626s到.2s/56accepted，
与4档.2s/40accepted同终时比较；各24事件组，30条接触边界条目ID/顺序相同。
新27后处理PASS，旧组合38PASS+1历史RED，实际reader PASS_COMPLETE_READBACK。
10项共同17时刻归一化差均<.807%，不是收敛阶/实验误差/AFSD认证。
Ws-Es-P仍-.002018197J告警，能量未认证，旧闭循环/空间FAIL保留。

下一只做旋转刚壁接触运动学/功接口，不再回头调此小块：
- 读src/afsd_backend/pysph/geometry.py的tool_pose/local_xy/surfaces；
  定义per-contact ID/feature/point/normal/wall velocity与刚体pose接口，
  区分旧surfaces粒子位置壁速与实际力作用位置，先核对功共轭。
- 独立测试v_wall=V+omega cross r、共转零滑、静壁/移动旋转壁功、
  旋转方孔开放与feature变化、标签不影响运动学；先冻结参考与范围。
- scripts/probe_process_friction_cone.py的wall/state/vector_contact_work
  只支持统一平移/固定法向；只做命名companion，不静默改生产默认。
- 接口审查后才考虑有界旋转小组件运行；正常开离、局部送料、热/FKS、
  有限应变仍分开列缺项。空间/第三档/AFSD/生产NOT_RUN，不启动GPU/长跑。

## 当前PredictorOnset：塑性—摩擦完整小块.2s已贯通（2026-10-02）

先读PROCESS_PREDICTOR_ONSET_REPORT/CONTRACT/REFERENCES，以下为历史。
新PredictorOnsetBlock仅branch接旧predictor，advance/commit/材料/空间/物性不变；
预测与cone仍共享14修正，原事件/功/法向检查不放宽。未改旧源码/配置/FAIL。
新38PASS+1旧RED（8独立、5实际、25元数据），旧predictor13+RED/Onset27+RED。
唯一全程4增量76.738s：16请求经24事件形成40accepted步/41帧到共同.2s。
实际reader PASS_COMPLETE_READBACK；不是AFSD沉积或材料能量认证。
Ws-Es-P=-.004347974J仍告警，旧闭循环/空间FAIL保留；8增量等NOT_RUN。

下一：configs/benchmarks/process_predictor_onset.json保持同一物理问题，
用scripts/probe_process_predictor_onset.py跑8增量，输出到新目录；先冻结
同终时反力/位移/ep/P/Q及事件测量与工程容限。预算先估90–150s、<=180s，
不得以不同失败终时算收敛。reader用inspect_process_predictor_onset.py；
不再回头针对已过的.1125初猜循环调参。完成有界对照后讨论旋转壁/局部送料；
正常开离/有限应变/热/FKS/开放旋转方孔仍未在此基准验证，不授权生产。

## 当前可行初猜：分支阻塞已越过，下一接回既有起滑事件（2026-10-02）

先读PROCESS_FEASIBLE_PREDICTOR_REPORT/CONTRACT/REFERENCES；以下Onset为历史。
新纯分支solve_with_feasible_predictor仅预平衡非接触自由DOF，复用原切线/力，
预测与cone共享14修正，不裁N、不改变accepted或材料。不改旧源码/配置/FAIL。
独立13PASS+1旧guard RED；旧Onset27+RED/Q7R39/Q6 31兼容PASS。
一次冻结.125目标对照.675s：旧initial N=-6.972N；新1预测+4cone返回根，
minN=14.599N、root2.924e-10。八个强制黏着面超摩擦圆，故exit2且全候选REJECTED。
accepted仍.1125/16步，完整路径NOT_RUN。输入/状态指纹均未变。
下一项：companion接回Onset现有事件驱动，先验证起滑定位和原子commit，
不直接接受超限分支或整段强制释放，不增预算/换本构/加正常开离。
完整.2s重跑前重估成本；8增量须共同终时，AFSD生产仍不启动。
旧材料能量警告及闭循环/空间FAIL保留。

## 当前Onset v1：已越过起滑阻塞，接受至.1125s；初猜法向域阻塞（2026-10-02）

先读PROCESS_PLASTIC_FRICTION_ONSET_REPORT/CONTRACT/REFERENCES和
FAILURE_DIAGNOSTIC_RESULT。以下Q7R为历史。旧代码/配置/失败未改。
新advance_onset用正容量近边界资格＋整体探测区分即时起滑和明确向内卸载；
资格不是强制滑动。原J2/整体切线/14预算/严格内部事件/原子commit均保留。
27新PASS（10独立接触、3实际小块、14来源/trace）＋1原RED保留；
旧Q7/Q7R39和Q6的31兼容PASS。

唯一实际4增量45.687s接受16区间/7事件至.1125s，原来是.05s/6区间。
物理字段独立部分读回PASS，不能称完整.2s程序通过。8增量/空间/生产NOT_RUN。
下一.125s目标触发negative normal；.416s纯重放确认一次trial、
iteration0尚无Newton修正。ID19 N由accepted+212.145N到初猜-6.972N，
gap全0、机械残差.1853；不是已平衡开离。REJECTED诊断与accepted严格分存。

下一步只查可行初猜/有界局部延续：参考ConeBlock.branch/state的初始位移设置，
复用当前整体切线和约束，不裁负反力、不增加Newton或越过原commit。
先新增test_process_plastic_friction_feasible_predictor.py，独立SPD线弹性
消元参考须同时覆盖“初猜负N但平衡正N”和“固定接触平衡仍负N”；数值见REPORT。
之后才用冻结.1125状态做禁止提交/有限局部对照；没有正法向根证据前，不直接加
脱离再接触模型。若确需正常开离再单独设计法向主动集，不能强制拉住材料。
材料能量告警-4.3466mJ与旧ep/闭循环/空间FAIL保留；温度/热/FKS/真实AFSD不加。
完整路径成本需重新估计，未获共同.2s前不做加密比较。真实旧AFSD仍2us/191步。

## 当前Q7R：有界14预算已接入，事件定位仍BLOCKED（2026-10-02）

先读PROCESS_PLASTIC_FRICTION_Q7_RECOVERY_REPORT/CONTRACT。以下Q7为历史。
新增RecoveryBlock/独立reader，复用旧求解、事件、commit与输出；旧源码未改。
新18测试PASS（6核心含实际.025预压事件、12元数据）；旧Q7的21及Q6的31兼容PASS。
唯一完整尝试10.198s，越过原Newton12阻塞后报cone event not bracketed；
accepted仍.05s/6步，与旧Q7全部30字段逐字节一致。无.05后接受事件。
只读部分读回PASS，≠整段PASS；原材料能量告警与全部FAIL保留。8增量NOT_RUN。

禁止提交的事件trace已实测：37root成功，34strict中点皆失败，width已<1e-12s
而lower=None。ID21探测越界6.095e-5N在9e-5N内未释放，但strict不允许负余量。
独立audit/完整trace已保存，accepted指纹未变。下一步不再加Newton或二分预算，
应区分BOUNDARY_ONSET与INTERIOR_EVENT；独立基准必须覆盖
圆锥边界起始、多接触方向改变、微小带符号残差、拒绝后无状态污染。
之后才做最小事件处理修复，保持原9e-5N/1e-11m等物理验收，不裁剪接触力/耗散，
不把REJECTED转accepted。保留本轮14预算失败，不能仅更改strict或无限细分试到通过。
补丁只限Q7 companion的初始分支分类/事件入口；先新建
tests/review/test_process_plastic_friction_onset.py的解析弹簧/耦合多接触回归；
建议独立参考与源码接口见RECOVERY_REPORT。新修复本轮NOT_RUN。
完整共同.2s和读回通过后才做8增量；热/旋转工具/送料/生产仍不进入本轮。

## 当前Q7：J2＋双切向摩擦已接入，预压段有效、完整路径BLOCKED（2026-10-01）

先读PROCESS_PLASTIC_FRICTION_Q7_REPORT/CONTRACT/REFERENCES及BRANCH_DIAGNOSTIC。
新PlasticFrictionBlock复用Q2实际材料更新/切线和Q6圆锥/整体耦合/事件；独立原子commit
统一累加ep/P/R及摩擦账。没有改生产src或旧失败，没有添加硬化/热或缩放。
新21测试PASS（7独立参考+3实际塑性事务/标签+11I/O）；旧Q6兼容31PASS。
实际4增量在6.026s内接受6步至.05s预压完，2接触事件，8接触因塑性横向膨胀滑移。
P=.239646J,Q=.00528460J；能量诊断Ws-Es-P=-.00428535J告警保留，不补热。
下一.0625s候选Newton12未达1e-8，BLOCKED。独立部分读回PASS不代表整段PASS。
8增量、空间/热/FKS/旋转工具/动态AFSD本轮NOT_RUN；未比不同失败终时。

新只读diagnose_process_plastic_friction_q7.py从accepted诊断数据重建：原12次精确
复现失败；固定14上限第13次分支root1.638e-10通过，但ID21/feature1强制粘着
超圆5.691N仍必须事件分段。两个诊断candidate均REJECTED，不能当成新accepted。
所有状态指纹不变；不是生产restart。不能只提高迭代后直接跳过接触验收。

下一步只做有界求解恢复/贯通既定塑性-摩擦路径，不加新物理：
1. 在Q7入口显式区分旧12次基线和新的有界恢复预算（最多14次的证据已得），
   或先做有界区间细分；物性/残差门槛不变，保留旧失败源版本及数据。
2. 复用Q6 advance_cone的粘着可行性与事件分段，特别捕获ID21越界；
   增加“分支收敛≠可提交”、失败后材料/耗散不污染的回归，不无限加Newton次数。
3. 只有完整共同.2s成功才做8增量敏感性。每次外timeout<=60s，先估成本；
   保留材料能量告警并区分求积/本构问题，不给热塑性自由能认证。
4. 之后再选实际几何/固定法向之外的接触扩展，暂不同时接热、供料或旋转工具。
旧Q4ep严格FAIL、闭循环/空间风险保留。真实AFSD仍2us/191步无持续沉积，
三工况/长跑/GPU未放行。以下Q6及更早章节均为历史，不继承其整段PASS为Q7PASS。

## 当前Q6：双切向圆形摩擦锥完成受限基准（2026-10-01）

先读PROCESS_FRICTION_CONE_Q6_REPORT/CONTRACT/REFERENCES.md。以下Q5及更早为历史。
新增scripts/probe_process_friction_cone.py，不改src/旧模型/旧失败证据。
27点/9接触/54自由DOF；Y和Z联合圆形Coulomb条件，整体法切耦合自然映射求根。
逐接触raw向量功、独立接触定律验收、粘滑事件分段；无逐轴裁力/补热。
3组均到制造.2s，各8事件；接受24/40/24区间，进程15.90/16.69/15.72s。
新31测试PASS（15核心+5guard+9I/O+2比较）；3组实际独立读回PASS；旧Q5兼容23PASS。
初轮共享DOF方向误差度量ERROR及两个mask漏检/一次测试时点错误保留。
固定原force/slip阈值，方向检查改用量纲一致传播，非调整物理；细节见合同附录。
仍是小应变弹性/固定法向，不能把它称为工具转动、热塑性流动或AFSD沉积。

下一Q7只接“已有J2塑性 + 已验证双切向摩擦”，在同小块/固定法向范围内。
预期新companion scripts/probe_process_plastic_friction.py、独立配置/测试，不直接改生产。
复用probe_process_plastic.py的trial/tangent/material状态以及Q6cone_response/cone_rows/
advance_cone/vector_contact_work；实现单次原子commit统一记录ep、plastic_work、
friction_work与accepted/rejected状态，不能在试探期间累积塑性或把两种耗散混作热。
先独立弹塑性材料点+带摩擦的串联系统/整体耦合参考、卸载/反向和功分账，再预算<=60s小块。
Q7先明确常屈服J2制造假设，不冒充FKS；不同时引入温度、供料、旋转法向或工具几何。
旧Q4严格ep尾差FAIL和闭循环/空间风险仍保留；Q7若触发须报告而非松阈/裁值。
正常脱离再接触、空间n5、有限应变、热/FKS、持续供料/旋转方孔、可重启生产仍待开发。
真实AFSD仍停留2us/191步；不得把Q6的.2s准静态制造时间当作AFSD物理进展。
本轮无三工况/长跑/GPU。先完成机械耦合最小例，再决定真实几何接入，不重比仓库。

## 当前Q5：短切向摩擦/反向受限基准已完成（2026-10-01）

先读PROCESS_FRICTION_Q5_REPORT/CONTRACT/REFERENCES.md；以下Q4及更早均为历史。
新增companion scripts/probe_process_friction.py，复用Q3/Q2实际材料trial/tangent/配对力。
27点、9X+接触，明确仅X=0三分量夹持，Y定向Coulomb、Z无摩擦，法向预压保持。
整体法切耦合解，不按P/m或壁速指定摩擦方向；粘滑边界事件分段、逐接触rawQ记账。
4/8每相增量全到制造.2s，各6事件；完整进程8.30/9.17s，外timeout60。
新核心13+guard3+读回7=23PASS。旧Q4兼容18PASS/1原ep严格FAIL，未改阈值/日志。
微小负Q约-2e-19J原样保存，不裁热。最终材料ep0，摩擦维持的受载残余位移
不是塑性残余变形。两档敏感性小不等于空间/本构/论文验证。

下一Q6：同companion小基准扩展到二维切平面Coulomb锥，固定法向/弹性先行。
预计新scripts/probe_process_friction_cone.py及配置/tests；以Q5 branch/state/commit
共享接口扩展tangent_force/relative_increment为二维向量，保留整体法切/接触耦合。
不可逐轴独立限幅；方向由相对运动求出。先独立斜滑、坐标旋转不变、零滑移、
反向和功参考，再<=60s三维小块。不同时加热/供料/旋转几何或改生产默认。
J2+摩擦、脱离再接触、空间n5、真实工具/热/有限应变还需后续验证，Q6当前NOT_RUN。
原AFSD仍2us/191步无持续沉积；原闭循环/空间FAIL保留；三工况/长跑/GPU未授权。

## 当前Q4：双壁共享接触与部分脱离已完成受限基准（2026-10-01）

最新PROCESS_CORNER_Q4_REPORT/CONTRACT/REFERENCES.md；以下Q3及更早章节均为历史。
只新增companion，src/旧配置测试/失败证据未改。27点、18接触行、15个接触粒子，
公共边3粒子各两个feature；共享整体trial/tangent/commit，逐壁算力/间隙/功。
弹性交错、J2对称、J2交错8/16四组全部到同0.2s；不是AFSD过程。最长进程6.136s，
每组outer timeout60s。J2交错先X约.132143s释放、后Y约.166071s释放，材料不跟壁回零。
核心12PASS/1FAIL：卸载ep最大增4.918e-12超过原1e-12严格断言。FAIL/阈值保留，
不裁塑性/热、不反复调参。增量加密净功差缩小，但ep尾差变大；不能称全面收敛。
四组独立读回通过，旧Q3兼容19PASS；新reader事件清单完整性有RED→最小检查修复。
读回最终6PASS，初轮缺失事件1FAIL日志保留；新测试总18PASS/1严格FAIL。
完整证据以Q4报告/日志为准。生产默认未替换。

下一Q5仅设计/实施短切向接触：在同一companion共享整体柔度下增加法向/切向
状态接口和分面相对功，先独立粘着/滑移、移动壁零滑移、反向基准，再预算内小块。
预计新scripts/probe_process_friction.py及配置/tests；复用Q4 CornerBlock或Q3
ContactBlock/Q2 trial/tangent/commit，不先改src生产contact，不移植未验证旧限幅。
必须先定义独立参考、单位、适用域和单次<=60s预算，不能仅要求Q>=0而漏掉反向。
Q5现在NOT_RUN；本轮不同时接工具旋转、热、供料或生产三工况。
任意法向/新碰撞/空间n5/热/FKS接触/有限应变/可重启生产均未验证；旧空间及本构
闭循环风险保留。真实AFSD仍2us/191步、无持续沉积；长跑/GPU/生产未授权。

## 当前Q3：单边法向接触及脱离已跑通（2026-10-01）

先读PROCESS_CONTACT_Q3_REPORT/CONTRACT/REFERENCES.md。只新增companion，src和旧证据未改。
27点/54永久自由DOF；弹性及J2两档增量全COMPLETE至制造问题0.2s；不是AFSD。
J2在约.132142857s自行卸载至零接触力并脱离，终材料面残余缩短27.142857um，
没有跟壁回零。用整体约束切线解接触耦合，不用P/m；事件分段后独立核对壁/材料功。
最长完整进程2.157s；核心13PASS、读回6PASS，复用Q2的23项本轮也PASS（不重复计新测试）。
初轮事件metadata提交回滚FAIL已最小修复、日志保留。未改物性/容限/积分让旧失败变绿。

下一Q4只做两个相交壁面的法向耦合/脱离：扩展companion的接触候选/法向映射，
先独立多接触约束参考，再三维角部。共享整体材料trial/tangent，逐壁存力/间隙/功，
提交/回滚保持统一；先定义预算<=60s，不同时加入摩擦、热、供料、旋转工具。
通过双面后再做短切向运动/摩擦和实际工具几何。Q4现在NOT_RUN，生产默认不替换。
本轮空间敏感性NOT_RUN；Q2空间误差、原本构闭循环风险继续保留但不无限纠缠。
固定参考、小应变、制造J2/nu0限制明确；新变体不能冒充有限变形AA6061热加工。
真实AFSD仍2us无持续沉积；三工况/长跑/GPU未授权。以下Q2及更早章节均为历史。

## 当前Q2：独立受载塑性已接通，下一单元只做单边接触（2026-10-01）

先读PROCESS_PLASTIC_Q2_REPORT.md/CONTRACT.md。本轮新增companion driver，src未改。
J2三组与原FKS一组均完成0.2s压缩—回程，27/125点、16/32增量，最长进程2.717s。
复用PySPH/原material.update/成对corrected算子，真实整体平衡求自由位移；
先trial/tangent后一次commit，分别记反力功、塑性功和返回诊断，不补热。
核心15PASS、读回/篡改回归8PASS；最小修复reader漏查账本帧数，RED日志和旧FAIL均保留。
空间反力偏差约22.11%不单独卡下一机制，但不称空间收敛。原FKS率下限432/432触发，
低率材料规律、真实AA6061参数和热学仍未验证。0.2s是制造基准，不是AFSD过程时间。

Q3下一步：新scripts/probe_process_contact.py/新配置/tests，复用Q2 trial/tangent/commit；
将双向端面位移约束替换成刚性压头的单边法向接触，先单面压入—脱离/反力—间隙—功，
再两个面的整体耦合/约束/回滚；接触柔度必须来自整体切线，不用单点P/m冒充。
先预定义独立参考、单位和<=60s预算，再运行。Q3当前NOT_RUN。
不一起接旋转摩擦/入口/热，不先替换生产solver或本构，不为了精度反复调参。
完整AFSD仍2us/191步、无沉积；全116回归本轮NOT_RUN；生产/长跑/GPU未授权。
下面Q1及更早章节是历史，由本节决定当前下一步。

## 当前Q1：独立准静态受载基准已跑通（2026-10-01）

用户已确认按粗精度工程近似继续；不要求先将所有误差修绿。
最新PROCESS_TIME_Q1_REPORT.md/CONTRACT.md。仅新增独立driver/config/tests，src未改。
27/125/343点三维块到同0.1s加载终时，4/8增量，最大完整进程0.637s；
13项回归/读回PASS，反力与功成立。不是AFSD已推进0.1s，也没有替换生产积分器。
应力偏差26.55/40.98/34.13%非单调，作为风险记录，不因粗精度阻断下一机制开发。

下一单元Q2只做受载塑性/挤压：独立J2参考及真实physical_dt下的材料状态提交/回滚，
再在相同独立driver显式接现有FKS；保留其温域/率正则化与b1来源局限。
预计改scripts/probe_process_time.py的小型材料接口或伴随driver、新Q2配置/tests，
src/material.py、mechanics.py和AFSD默认solver不先替换；不并行做全局接触/入口/热。
先实际测试塑性耗散、卸载与不可逆应变，再做整体多面接触柔度和反力功。
旧应力/闭循环FAIL保留；基本有限值/不可穿墙/守恒检查仍保留。
Q1没有求惯性/接触/热/FKS/大变形；新AFSD/三工况/GPU均NOT_RUN，不授予生产。
下面Q1以前的“先取得方向确认”已由本轮继续请求完成；历史计划不覆盖本节。

## 当前E3-4：原2us窗口完成，停止微秒续跑（2026-10-01）

优先读CONTINUATION_E3_REPORT.md、PROCESS_TIME_DECISION.md。
未改src/配置；从合法48步/.5us续到191步/2us，新增143步用413.07s，单次exit0。
链核验192帧；首次3个自由棒粒子塑性，累计塑性功6.0569e-9J；
它们初始在肩面下，供料穿越/孔外候选体积仍0，不是沉积或70%完成。
账面残差升至3.3109%，旧材料/边界精度FAIL保留，不能因为没报警视为高精度。
新增链/观测/热账28项PASS；完整116项本轮NOT_RUN（引用111PASS/5旧FAIL）。

下一任务不再针对476/439，也不继续无休止续跑：先确定工艺时间方法可行性路线。
当前dt约10.49ns、2.803wall秒/步，转1度约41.23小时、送料50um约145.86天。
建议取得独立隐式/准静态基准方向确认后，在scripts/tests内做固定三维受载块的
力平衡/反力功/试探状态回滚验证；mechanics.py/material.py/solver.py默认暂不改。
新方法必须复用PySPH设施、保留真实时间的率依赖，不能把P/m局部接触直接当整体柔度。
接口/源码范围/基准和结构性缺项见PROCESS_TIME_DECISION.md；新变体目前NOT_RUN。
CFL/空间同终时敏感性未跑且预算超限，不以异终时称收敛；生产/GPU仍未授权。
以下E3-3及其“下一步”均为历史；HANDOFF和本节优先。

## 当前E3-3：0.5us真实集成通过（2026-10-01）

先读SHORT_WINDOW_E3_REPORT.md。未改src/物理，单次48步到.5us，122.24s，
所有接受诊断/逐帧/源配置/身份/7项读回检查通过；ep和穿越质量0，尚无沉积。
停止针对476补丁。下一单元核对本次合法checkpoint和预算后续到原2us：
余约143步、预测364s，加20%约437s；建议内540/外570s，当前NOT_RUN。
0.5us处截短步的时间网格须记录，短窗口/续跑不是dt收敛。
不自动启动长跑/生产/GPU；旧材料/边界/E2失败不变。后文是历史。

## 当前E3-2：局部接触尺度修复已通过（2026-10-01）

先读CONTACT476_SCALE_REPORT.md。只改contact_cone.solve的逐面正标量预条件，
QP/物理/容限不变；冻结476参考力及完整kick、新64组/旧432组全部通过。
完整回归111PASS/5旧FAIL。初版15组失败证据保留，不能隐藏。
真实2us NOT_RUN：预检首步3.576s、191步预测683.10s>600s；
不能据局部PASS将旧E3真实FAILED改判完成，也不能恢复旧REJECTED。
下一单元是预先界定预算内的集成窗口或分离初始化/稳定步计时，
不是继续针对476调参；不改变原物理或把短窗口当作同终时收敛。
下面E3-1及更早章节均为历史。生产仍未授权，宏观沉积尚未展示。

## 当前E3-1受限集成（用户澄清后，2026-10-01）

实测结果：仅一次2us尝试，21步/.220296us后candidate22/476双面cone证书失败。
不是E2应力精度拦截；应力精度工作降为独立风险轨，不阻止受控实现探索。
下一轮只修近停滑/多面尺度的数值稳健性，先冻结476独立参考和原接触回归；
不得直接调高证书容限、增加迭代或改变物理来过关。详情ENGINEERING_INTEGRATION_E3.md。
ep=0、进给面穿越0，尚无塑性流动/沉积；源和参数未改，拒绝末态不用于重启。
下面是本轮运行前安排，现已执行；不是授权立即再跑或延长终时。

用户优先工程过程重构；E2应力FAIL保留但不单独阻断短程APPROXIMATE集成。
本轮按ENGINEERING_INTEGRATION_E3.md：重跑前置回归并估成本后，仅原固定问题2us一次。
不改物理或容限、不切后端、不再自动执行下述离线E2分解；保留作为后续精度工作。
遇接触/能量/温域/穿透/NaN失败即停；不自动换CFL、材料或加稳定化重试。
生产仍禁止；实现探针不等于沉积/混合已出现，下一步由实测状态决定。

## 当前E2结论（2026-10-01）

ENGINEERING_BLOCK_E2.md：5档力驱动三维小块到共同1.963961us，数值约束通过但应力验收FAIL。
时间差异比.2572，空间应力28.66/30.71/27.53%，不单调且最细>20%；不可进入E3/生产。
下一轮仅做离线误差分解：解析位移应力恢复、解析应力/加速度节点残差、实际应力与E*sym(Du)之差。
固定原n5/7/9数据及X<=L/4物理带，不加时间算例、不调参数、不改生产core或重写通用SBP。
92.74%应力平方误差在支承前三层是SUPPORTED边界嫌疑，尚非唯一根因证明。
7新回归、E1 6项、空间4项及5组文件读回通过；旧边界/弹性闭循环失败证据保持。
新输出/日志/配置完整；实际调用和小变形变体局限见报告，不能冒充生产分裂积分器验证。

## 当前路线：工程近似 E1 → 受载基准 E2（2026-10-01）

用户已接受工程近似，不要求逐式完整复现；ENGINEERING_APPROX_CONTRACT.md为新目标。
E1本轮只用现有material.update完成独立三维等体积加载/卸载/反向塑性材料点基准，
不替换FKS、不修改src/、不声称AA6061已标定。证据与命令见ENGINEERING_MATERIAL_E1.md。
旧闭循环/边界FAIL保留；上游实际loop限制新增RED，禁止原样移植impact3d塑性/能量。

下一工作单元E2：先冻结非AFSD受载可变形小块的物理域、自由面、壁/外载、物性、
dt/空间网格、应力/位移/反力和功误差定义；独立连续体参考不得来自被测SPH算子。
复用mechanics.py中的corrected梯度/伴随内力对及material.update，
粒子/NNPS/ID/输出沿用当前PySPH；优先在scripts/独立驱动，生产solver.py不切换。
如需测试专用材料接口，明确注入benchmark law，不修改FKS默认参数/生产JSON。
tests/review/新加基准：不可将所有自由点规定运动，不得用负内力造外力或残差热。
先量化旧边界误差对小变形受载响应的影响及同几何细化趋势，再决定有限范围可用性。
E1不解除一般变体积Jaumann能量风险；热耦合/大变形/AFSD须另过关。
本轮真实AFSD、完整回归、GPU、生产三工况与新空间敏感性NOT_RUN。
以下保留为历史，不能据其“下一步”越过当前范围。

## 2026-09-30第三轮：方法类门禁，不再逐粒子推进

最新报告METHOD_CLASS_REVIEW.md；合同CONTACT_CLASS_CONTRACT.md。
局部接触432固定锥问题、439独立锥/完整冻结kick和旧事件基准通过。
完整回归116项111PASS/5FAIL，历史反例与当前材料/边界FAIL分开解释。
下一工作单元：离散边界弱形式与牵引接口相容性，先独立常应力/非仿射基准，
不把外力设为负内力强行过patch。固定域3/5/7探针未证明节点一致性。
本构保持原Jaumann APPROXIMATE；有限应变能量验收BLOCKED，自由能候选仍仅设计。
两项方法门禁解决后，再定义非AFSD三维变形/摩擦/导热耦合基准；当前NOT_RUN。
本轮无新真实AFSD2us/敏感性/生产；不以局部接触PASS授权长跑。

## 2026-09-30第二轮：454已修复，下一步仅439双面锥

最新证据见CONTACT454_REPAIR_REPORT.md、CONTACT454_REPAIR_CONTRACT.md。
454的近零投影/时间可表示性修复通过独立基准；107项101PASS/6FAIL。
唯一新真实试跑接受26步（0.272746us），候选27/439双面锥慢收敛而拒绝。
先运行test_contact439_open_risk.py；该冻结子问题有独立闭式KKT解。
后续最小范围仍仅contact_kick.py::cone_solve：保留T_aᵀ(P/m)T_b，
研究有原容限残差验收的局部active-set/KKT处理，不能只加迭代或放宽门槛。
再过454、反转、静/动墙、多面执行器、功恒等式和完整回归后才重新评估真实运行。
材料/空间/物理参数不改；旧RED和冻结不覆盖；生产与CFL/细网格仍BLOCKED/NOT_RUN。

## 2026-09-30第一轮历史：下一步只修局部RED

先读COUPLING_FAILURE_REVIEW.md。原负热首发32/485已定位，局部修复有独立参考；
event_aware_v1真实固定问题在3/454近零滑移锥迭代失败，禁止生产/继续CFL序列。
下一最小范围：contact_kick.py的cone_solve/kick；用稳定切向基和一致的
近零事件/KKT尺度处理，保留P/m多面耦合、反力功、每子区间原始Q。
先跑test_new_candidate_red.py，再跑全部run_coupling_review.py，保留原失败。
未解决前不增加迭代上限/放宽门槛/截热，也不换材料或空间算子。
闭循环非零功已解释；FINITE_STRAIN_VARIANT_DESIGN.md仅设计和NOT_RUN测试规格。
常应力边界节点patch FAIL须保留，corrected伴随对不改；其修复另行设计。
修复接触后重测成本，先同一dx=.5mm/CFL=.1到2us；只有达标才按既定顺序细化。

## 下列为上一轮历史计划（非本轮扩范围授权）

先读SCIENTIFIC_REVIEW.md与VALIDATION_CONTRACT.md，旧计划不构成生产授权。
1. contact.py/solver.py：多力/多面一致粘滑与离散功，补外圆/顶面/棱角；不得裁负热。
2. material.py/solver.energy/output.py：预应力压缩自由能一致性，独立循环/压缩/旋转基准。
3. mechanics.py/stable_dt：完整性修正谱稳定界/自由表面；同终时dt与同物理BC空间敏感性。
4. 解材料/启动来源与守恒入口；通过后再编译加速、性能测量及三工况。
5. 精确ROI/截面面积算子及独立数据验证；只读可视化先验hash，不改变物理。
现有本轮命令、反例与成本见审查报告；默认600s/次，无长算例授权。

# AFSD implementation plan

## Current handoff: implemented DEBUG path, next gates (2026-09-24)

The preparation-only plan below is historical. The current implementation is in
src/afsd_backend/pysph/{kernel,material,mechanics,thermal,contact,geometry,particles,
integrator,solver,config_adapter,application,output}.py. Read PHYSICS_STATUS.md and
PARAMETER_SOURCES.md first. The equation-to-test map was created before coding in
EQUATION_CODE_TEST_MAP.md. No upstream source was changed; PySPH remains the only backend.

Actual acceptance: scripts/run_physics_checks.ps1, RunId 20260924_physics_verified,
exit0; 19 old +21 new tests. Three 657-particle real 3D coupled probes, dt-half,
restart, deliberate alarm, wall limit and independent time-series readback ran.
This is HOST_FKS_DEBUG, NOT a cleared paper configuration or convergence result.

| Next source/interface | Required work and tests | Blocking gate |
|---|---|---|
| material.py::update, solver.py::evaluate, integrator.py | Audit stress-work conjugacy, first-order thermo-plastic splitting and contact switching; dt/2, dt/4 with identical physical end time; material work and coupled elastic/plastic energy oracles | dt-half changes peak plastic strain ~69%; no time convergence claim |
| config_adapter.py, configs/debug, material.py::FKS | Obtain P40/P41, resolve b1/rate/T normalization; temperature-dependent k/Cp tables and bounded interpolation tests; parameter-range sensitivity | Paper values unresolved; floor=1/s and dimensionless b1 are DEBUG only |
| contact.py, geometry.py | P44 triangle-face neighbor adapter, robust closest-point edges/corners, persistent stick/slip state, shear cap; oblique impact, wall-work balance, rotating-bore corner and thermal contact area tests | Current finite-plane memoryless sliding is not paper mesh/stick law |
| particles.py, solver.py | Finite vs infinite inlet interface: persistent new IDs, injected mass/momentum/enthalpy/actuator work ledger; deformable-deformable separation/joining model tests | append is only identity-tested; initial shared-neighbor interface is pre-bonded |
| geometry.py, config_adapter.py | Real dimensions, clamp/contact extent, startup/preheat/dwell/lift/termination with sourced program; nonzero dwell/lift currently fail closed | Current bar is scaled and preheated; no macro deposition |
| mechanics.py, kernel.py | Surface completeness, tensile instability/stabilization, heterogeneous h and resolution convergence; objective 3D patch and conservation after every change | Only short uncorrected-kernel bulk/patch checks passed |
| equations/ (new, same PySPH backend) | Translate tested host operators into compiled PySPH Equation stages; keep exact scalar oracle, stepwise field/digest comparisons, CPU/GPU capability separate; redo performance estimate | Host loop extrapolates ~3.2e4 wall hours/physical second at nominal paper N |
| output.py, schemas/, tests/backend | Versioned wider-ID/checkpoint migration, full source/parameter provenance, failure injection, field-time staggering tests | Windows upstream long is int32; no silent ID wrap or incompatible restart |
| analysis/ (not started) | Predeclare Fig8 area operator, ROI/section/AS-RS mapping, initial tracer selection, uncertainty and experimental source | Qualitative NOT_RUN; quantitative BLOCKED |

Do not extend simulated duration or launch paper scale before the numerical and
performance gates above. Every trial remains bounded to 600 seconds by default;
no long or paid computation is authorized by a terminal objective alone.

## Archived preparation-stage plan

## Stage result and boundary

The measured backend gate is closed: use only PySPH at
`third_party/PySPH`, pinned to
`1cd3dea9cfaac7d9f1bda2a2491ce4b42a31dbe5`. WeldForm CPU and
WeldFormGPU remain evidence/reference checkouts; do not develop either in
parallel.

This preparation stage did not implement a new constitutive, thermal or
contact model and did not run an AFSD case. It established:

- an isolated PySPH environment and successful upstream generic, 2D-solid and
  3D-solid bounded runs;
- three provenance-bearing case configurations with a fail-closed unresolved
  parameter gate;
- manifest/particle schemas, a canonical field contract and test entrypoint;
- an independent 118-point `GEOMETRY_FIXTURE` with exact VTK/ID/label
  round-trip;
- 19 passing scaffold tests, logged in
  `logs/tests/20260924_final_scaffold.log`.

The upstream `impact3d.py` one-step pass proves an executable 3D neighbor,
velocity-gradient and full stress-rate call path. It does not prove a correct
3D plastic law: that example calls `VonMisesPlasticity2D`, whose yield
invariant omits z-related terms.

## Selected-backend invariants

1. Keep `third_party/PySPH` at the pinned commit. Prefer project adapters and
   custom Equation classes; if an upstream patch becomes unavoidable, preserve
   it as a small reviewed patch series and record a dirty-tree digest.
2. Every solver launch resolves a case with
   `afsd_scaffold.config.load_resolved_case()` and calls
   `assert_physics_ready()`. No missing physical value becomes zero.
3. Calculation state stays SI. Display conversions are downstream-only.
4. `tracer_type` is passive metadata and must never affect neighbors, force,
   stress, material, thermal or contact calculations.
5. Output is keyed by an immutable `particle_id`, not a ParticleArray index or
   transient ordering. All fields follow the same sort/permutation.
6. A `DEBUG` run may use a reviewed, explicit required-parameter subset.
   `PAPER_COMPARISON` always requires the complete resolved tree.
7. PySPH's default JIT cache currently lands in
   `<LOCAL_USER_PROFILE>\.pysph`; the runner must record the cache/compiler state
   and either adopt a documented project-local mechanism or declare this
   reproducibility exception. Never silently reuse an unknown cached module.

## Next-stage source and interface map

All `src/afsd_backend/` paths below are planned new project code. No complex
physics is implemented in this stage.

| Planned path | Interface / responsibility | First required tests |
|---|---|---|
| `src/afsd_backend/pysph/config_adapter.py` | Convert one fully resolved SI case into typed backend settings; return a complete list of consumed config paths and provenance | Config-to-backend golden files; unknown/unit-mismatched/unresolved values fail; consumed-path coverage |
| `src/afsd_backend/pysph/application.py` | Single PySPH `Application` entrypoint, explicit `dim=3`, deterministic seed/thread settings, bounded run class and failure status | Startup gate; max-step termination; non-zero exit on NaN/missing field; command/commit/config digest in manifest |
| `src/afsd_backend/pysph/particles.py` | Create substrate/feedstock arrays and immutable `particle_id`, initial position, material, tracer and component properties | Unique IDs; sorting/update/restart/inflow mapping; tracer-independence test |
| `src/afsd_backend/pysph/geometry.py` | Build reviewed deformable feedstock/substrate and rigid tool/base geometry; maintain laboratory/tool frames | Counts/bounds/normals; square bore remains open; tool pose convention; fixture is never accepted as solver geometry |
| `src/afsd_backend/pysph/output.py` | Adapt PySPH state to `run-manifest.schema.json` and canonical particle records; atomic status transition; optional VTK series | Schema validation; NPZ-to-canonical round-trip; tensor ordering; checksum; failed/incomplete run; restart continuity |
| `src/afsd_backend/pysph/equations/kinematics.py` | Wrap/review `VelocityGradient3D`, kernel-gradient correction and genuine 3D neighborhood requirements | Non-coplanar neighbor gold set; affine patch; zero-distance case; 1/many-thread repeatability |
| `src/afsd_backend/pysph/equations/solid.py` | Auditable full 3D stress/objective-rate update; initially reproduce a verified elastic baseline | Uniaxial, three shear planes, volumetric/deviatoric split, rigid rotation of prestress, energy consistency |
| `src/afsd_backend/pysph/equations/plasticity.py` | Full 3D return mapping and material state; later a separate FKS flow-stress Equation once inputs are resolved | 3D J2 invariant; elastic unloading; return residual; plastic strain/work; independent FKS scalar curve and limit tests |
| `src/afsd_backend/pysph/equations/thermal.py` | Temperature state, k(T), Cp(T), particle conduction, plastic-work source and explicit source ledger | Uniform-temperature invariance; two-body direction/conservation; analytic slab convergence; plastic work heating; no silent clipping |
| `src/afsd_backend/pysph/contact.py` | Deformable-particle to rigid-triangle contact, friction, moving tool/base and thermal/contact heat ledger | Plane normal/penetration; separation; stick/slip; edge/corner cases; force/work sign; two-sided heat balance |
| `src/afsd_backend/pysph/integrator.py` | Explicit paper-compatible forward-difference/leapfrog state ordering and stable-step reporting | Constant acceleration; convergence order; momentum/energy drift; thermal/contact coupling order; accepted-dt trace |
| `scripts/run_pysph_baseline.ps1` | Re-run the pinned upstream bounded baseline without changing physics; isolate outputs and save exit codes | Generic and `impact3d` smoke; output readback; reject wrong commit or dirty tracked source |
| `tests/backend/` | Backend unit/benchmark tests, kept separate from paper cases | Every capability above has a numeric oracle and non-zero failure |

Relevant pinned upstream review surfaces are:

- `pysph/examples/solid_mech/impact3d.py` — verified 3D example wiring;
- `pysph/sph/solid_mech/basic.py` — momentum, full stress-rate and energy
  equations;
- `pysph/sph/solid_mech/hvi.py` — existing 2D-named/2D-invariant plastic
  projection that must not be reused as complete 3D plasticity;
- `pysph/sph/integrator_step.py::SolidMechStep` and
  `pysph/sph/integrator.py` — verified PEC/EPEC state sequencing;
- `pysph/base/*nnps*.pyx` and `pysph/base/tests/test_nnps.py` — 3D
  neighbor implementations and brute-force reference tests;
- `pysph/sph/rigid_body.py` — rigid/DEM proof-of-concept reference, not an
  AFSD deformable-solid/triangle contact implementation.

## Implementation sequence and stop gates

### 1. Reproducible baseline wrapper

Re-run the exact selected baseline from
`logs/upstream/pysph/00_exact_commands_and_status.log`. Record environment
freeze, compiler/JIT cache, upstream commit, command, exit code, output count
and checksum in a run manifest. Add the missing test-only dependency only in a
separate reviewed environment update, then run the selected genuine-3D NNPS
test against its brute-force oracle.

Stop if the original `impact3d` one-step run or NPZ readback no longer passes.

### 2. Identity, configuration and output adapter

Implement `config_adapter.py`, `particles.py` and `output.py` before
physics changes. Reuse the current config and output contract tests. Upstream
NPZ files in `outputs/upstream/pysph/` are evidence artifacts, not
schema-conforming project results.

Stop if sorting/restart changes the
`(particle_id, initial_position_m, material_id, tracer_type, component_id)`
mapping.

### 3. Existing 3D elastic foundation

Write small numerical tests around the actual PySPH Equation sequence used by
`impact3d.py`. Establish a trustworthy 3D neighbor, gradient, stress and
integrator baseline before adding plasticity. A successful 3D visualization is
not a test.

Stop on failed affine patch, objective rotation, energy sign or reproducibility
checks.

### 4. Full 3D plasticity, then FKS

Do not extend `VonMisesPlasticity2D`. Add an independently tested 3D material
state/update. Only after P41/P40 or an explicitly approved approximation
resolves FKS dimensions, reference strain rate and temperature normalization
may the FKS law be enabled.

Stop if the yield residual, plastic work or limiting behavior lacks an
independent numerical oracle.

### 5. Thermal and rigid contact modules

Implement conduction and heat sources with separate ledger fields, then rigid
triangle contact/friction, then their thermal coupling. PySPH gas artificial
thermal conduction and rigid-body sphere DEM are not acceptable substitutes.

Stop if either pairwise energy balance or contact work sign fails.

### 6. Paper-compatible integration and coarse AFSD DEBUG

Implement and convergence-test the required leapfrog/forward-difference
sequence. Only after modules 1–5 pass may `geometry.py` create a coarse
`DEBUG` AFSD case. Measure particle count, accepted dt, steps/s, memory,
simulated duration and projected wall time before any expensive run.

Do not start `APPROXIMATE` until every approximation is named and bracketed
by sensitivity runs. Do not start `PAPER_COMPARISON` while any required base
configuration value is `UNRESOLVED`.

## Validation ladder

1. Pinned upstream generic/2D/3D smokes and NPZ readback.
2. 3D NNPS and identity-permutation tests.
3. Elastic 3D patch/objectivity/integration tests.
4. Full 3D plastic material-point tests, then FKS tests.
5. Thermal analytic/conservation tests.
6. Rigid contact/friction/heat tests.
7. Small coupled 3D benchmark with mass, momentum and energy ledgers.
8. Coarse AFSD `DEBUG` performance probe.
9. Named sensitivity runs for approved missing inputs.
10. Only then, predeclared Fig.8/11/12 measurement operators and a
    `PAPER_COMPARISON` attempt.

## Blocking evidence/decisions

| IDs / area | Blocker | Evidence or decision required |
|---|---|---|
| A02–A04 | FKS dimensions/reference state; k(T), Cp(T) | Obtain/review P41/P40 or approve sourced ranges; independent formula/property tests |
| A05–A06 | Tool/substrate/feed geometry; startup/path | P40/supplement or named geometry/program approximations with sensitivity |
| A07 | Friction, thermal contact and heat partition | P40/P44 or sourced model/ranges; conservation and work tests |
| A08 | EOS, h/spacing, stabilization, CFL and scaling | Trace and test every adopted value; unscaled benchmark first |
| A09–A12 | Sampling, Fig.8 raw data, coordinates/AS-RS and measurement operator | Predeclare ROI/projection/selection; preserve digitization uncertainty and undefined zero denominator |
| PySPH plasticity | Upstream 3D example invokes a 2D yield invariant | New full 3D update and material-point oracle |
| PySPH thermal/contact | Required solid thermal and deformable-solid/rigid thermal contact paths are not available | Implement/test modules before AFSD coupling |
| Integrator | Verified upstream PEC/EPEC is not the paper's asserted leapfrog path | Explicit state sequence, convergence and stable-step evidence |

WeldForm's missing GSL and WeldFormGPU's CUB neighbor failure are recorded
candidate blockers, not work items for the selected-backend phase.
