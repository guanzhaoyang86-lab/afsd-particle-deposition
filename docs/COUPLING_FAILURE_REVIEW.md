# AFSD 首次耦合失败诊断与最小接触修复审查（2026-09-30）

后续更新：见CONTACT454_REPAIR_REPORT.md。本文保留第一轮历史事实；
454已在第二轮独立修复，新的真实阻塞是候选27/439双面锥慢收敛。
第二轮107项101PASS/6FAIL，仍APPROXIMATE/BLOCKED，未达到2us。

## 结论

已定位原始首发负耗散；建立独立参考并实现显式选择的局部接触修复变体。
**不是稳定耦合修复完成，更不是论文复现。event_aware_v1 仍 BLOCKED。**
最终回归101项：95 PASS、6 FAIL、0 ERROR；RED没有xfail或删掉。
只试跑一次新的固定dx=.5mm/CFL=.1问题：接受2步，第3候选步因局部接触锥迭代不收敛而拒绝。
未达到2us，CFL=.05/.025、dx=.25mm均NOT_RUN；不比较不同失败终时。

唯一后端仍为PySPH commit
1cd3dea9cfaac7d9f1bda2a2491ce4b42a31dbe5，上游clean，许可证/既有设施保留。
没有GPU、并行、生产、界面或材料/空间重写。本轮不重新比较后端，不改参数匹配论文结果。
所有文内相对路径相对于项目根；先读本报告和COUPLING_DIAGNOSTIC_CONTRACT.md。

## 1. 源码实际调用顺序

不是简单的“mechanics → contact → constraints → material → thermal”。
由上游LeapFrogIntegrator.one_timestep和实际host callback确认：

    stable_dt
      → PySPH LeapFrogStep.stage1：位置前半drift，速度尚未kick
      → solver.evaluate(t+dt/2, dt)
          → NNPS → corrected_operators：L、drho、B
          → material.update：S_new/ep_new试算（此时未提交粒子状态）
          → rho_mid、sigma_mid → internal_force
          → legacy：接触力按v_old计算
          → legacy：顶cap x/z约束改候选加速度，并记反力功
          → 只读逐接触诊断
          → legacy：以同一个v_old算摩擦力，
                    但以包含所有作用的(v_old+v_candidate)/2算热
                    ← 原始首发异常在这里
          → 提交S/ep/率字段
          → 导热 + 塑性热 + 接触热 → 温度
      → PySPH LeapFrogStep.stage2：
          v/rho/e scratch整kick → 位置后半drift
      → time/count更新 → 全局报警检查 → 接受/输出

关键位置：solver.py evaluate:52、material试算:62、内力:65、
legacy接触:75、热复算:92、材料提交:98、热提交:106；
integrator.py HostLeapfrog:18–32；
third_party/PySPH/pysph/sph/integrator.py:468；
integrator_step.py:716、722。引用位置以本轮冻结源码为准。

新event_aware_v1只替换接触/执行器候选计算块（solver.py:67–73）：
domain_kick → 共享自由DOF迁移率的多面Coulomb → 零滑移事件/转向子区间
→ 每个实际接受子区间冲量与中点速度做功 → 返回等效整kick加速度及热量。
粒子、NNPS、上游积分设施、材料和径向导热未替换。
子区间只解冻结几何/内力的局部速度问题，不暗示材料和位置也做了子步求解。
这项分裂近似及其全局能量误差仍必须在后续真实试跑中检测。

## 2. PHASE 1：原始首个负耗散事件（完整数据已落盘）

主证据：outputs/coupling_20260930/phase1_verified/
- report.json：observer off/on前31个接受态所有粒子属性及账本逐位相同。
- contacts.jsonl：每一个接受kick和被拒候选的全部活跃接触。
- last_two_kicks.json：完整第31个接受kick + 第32个REJECTED候选。
- first_negative_event.json：首个触发原守卫的接触事件。
- provenance.json：源码逐文件hash、配置文件hash/解析配置hash、上游pin。
- last_accepted.npz：第31步接受态。
- REJECTED_partial_integrator_state.npz：前半drift后真实部分更新数组。
- REJECTED_candidate_fields.npz：候选末位置/速度/应力/rho/ep等预测，
  temperature_not_advanced明确不是求解后的末温度。
- REJECTED_manifest.json：禁止重启和禁止作为物理末态。

诊断只额外读取/重算力，不向求解器回传任何数值。
原contact.py、material.py、mechanics.py、geometry.py、integrator.py均未改。
最终event_solver_smoke进一步将当前legacy第31步与该checkpoint逐位比较，PASS，
随后仍在第32步出现原异常。原物理和失败没有被诊断代码“修好”。

### 身份、时间、几何

- step=32；particle_id=485；contact_id=485:BORE_Y_MINUS。
- 阶段：半drift后、约束后候选速度已知、热复算守卫前；整速度kick尚未执行。
- 已接受时刻3.251964305249681e-7 s；dt=1.0490071451121287e-8 s；
  接触几何评价时刻3.3044146625052875e-7 s。
- 粒子中心[-4.9997827704e-4,-9.9997843670e-4,3.7500035447e-3] m。
- 法向[-1.0381124828e-5,0.9999999999461161,0]；
  gap=2.4802675369940404e-4 m，penetration=1.97324630059597e-6 m。
- 前一接受态、当前中点、候选末态预测均活跃；features始终只有BORE_Y_MINUS，
  无换面、无第二接触、actuator=false。

### 速度时间层与冲量

诊断中的v_after_internal等为**虚拟部分冲量和**，不是源码执行了三次顺序kick；
真实上游只有一次合成kick。这个区别不省略。

| 量 | [x,y,z]，m/s |
|---|---|
| v_old | [0.0342645607277,-0.155753608831,0.00201708381804] |
| 虚拟只加内力 | [0.0252718050856,-0.181170558943,-0.000814296870747] |
| 加全部接触后 | [0.0245435286750,-0.140844884265,-0.00283138068878] |
| 约束后候选 | 与上一行相同，只有舍入差 |
| 实际使用的壁速度（粒子中心处） | [0.0335352491048,-0.0157073028289,0] |
| 初始切向相对速度u0 | [0.000727857784673,7.555982514e-9,0.00201708381804] |
| 候选端切向相对速度u1 | [-0.00899301949767,-9.335765800e-8,-0.00283138068878] |
| 热计算使用的相对中点切向速度 | [-0.00413258085650,-4.290083774e-8,-0.000407148435373] |

全部冲量单位N s：

- internal=[-3.0350550292e-9,-8.5782206628e-9,-9.5559098246e-10]。
- normal=[-1.4128625511e-13,1.3609917754e-8,0]。
- friction=[-2.4565200233e-10,-2.5501440948e-15,-6.8076578859e-10]。
- other_contacts=[]。
- constraint=[4.65e-25,-3.31e-24,3.10e-25]，只是代数差分舍入。

沿初始滑移方向投影的速度增量m/s：
internal=-0.005715740732，friction=-0.002144388977，normal≈1.32e-23，
another contact=0，constraint≈1.33e-18；初始滑移速率0.002144388977。
原独立m|u0|/dt限幅恰好让摩擦单独停到零，却没有纳入内力造成的速度变化。
内力单独已经使旧方向投影反号；摩擦保持旧方向进一步推动该方向反号。

**不能把这个向量事件误称为真实一维过零：**
候选直线速度轨迹中，沿u0投影在kick的0.2728185229处变号
（候选时间3.2805831632e-7 s），到中点已负。
整个切向向量在该直线路径上最小模仍为0.001480157216 m/s，
发生于0.1428367624处；并没有完整向量降到零。
这些是被拒离散候选的时间层诊断，不是已接受的物理轨迹。

### 做功、尺度与接触点

- Q_raw=-Jt·u_mid=-1.292349487965555e-12 J（负耗散，未取abs或截断）。
- 法向接触实验室系功=-2.018344702876596e-9 J。
- 摩擦实验室系功=-6.945611547285534e-12 J。
- 内力功=1.1832898346454383e-9 J；
  约束功≈5.04e-25 J；ΔK=-8.420004797784418e-10 J。
- 局部kick功恒等式残差=7.37e-25 J。
  该local_energy_residual是**动能kick恒等式**，不是完整热力学残差；
  本候选还没执行热更新，且材料不存在一致有限应变自由能。
- 舍入估计3.8727483347e-24 J；负值/舍入≈3.337e11。
  原热守卫等效J阈值约1.049e-20 J，负值也远超此阈值。
- 真实几何投影点为[-4.9997570224e-4,-1.2480051904e-3,3.7500035447e-3]m，
  该点壁速度[0.04132723938,-0.01570722194,0]m/s；
  用L外推材料接触点速度也已记录，明确DIAGNOSTIC_RECONSTRUCTION_ONLY。
- 实际模型的力和速度均在粒子中心配对。只把壁速度改到几何投影点会引入伪共转摩擦。
  此处没有改施力位置；有限半径真实接触运动学仍属APPROXIMATE。

对六个问题的直接回答：
1. 第32候选kick，485，BORE_Y_MINUS。
2. 旧方向投影在合成候选区间内翻转；不是实际执行了先内力再接触，也不是整个切向向量过零。
3. 内力是首发驱动；旧摩擦方向/孤立限幅使之恶化。
4. 未换feature。
5. 源码一般确实在摩擦力计算后施加cap约束；但此粒子无执行器、无有效约束冲量。
6. 明显离散负耗散，不是浮点噪声。

## 3. 根因分级

| 级别 | 结论 | 证据 |
|---|---|---|
| CONFIRMED | 旧力按kick前滑移定向/孤立限幅，热却按全部冲量后的中点；首发内力改变方向，导致负Q | 原始逐接触JSON、contact.py force/evaluate、31步逐位一致 |
| CONFIRMED | 本构闭循环非零功不是只改Es字段就能消除 | 实际update四档时间步+独立积分0.146111J |
| CONFIRMED | 新变体第3候选锥迭代2000次仍两点振荡 | new_candidate_RED.json完整末四次迭代 |
| SUPPORTED | 新失败为近零滑移/极小子区间下投影消减与力收敛条件不匹配 | 系数4.0313e-21s、法向相对速度0.1016m/s，切向量约1e-14m/s；见后文 |
| POSSIBLE | 其他粒子/更长时间的多面约束、几何激活切换、分裂积分能量误差仍可能触发失败 | 当前测试覆盖有限，未把局部PASS推广为全局保证 |
| RULED_OUT（只限原首发事件） | 其他接触、换面、执行器、纯舍入负值是该事件起因 | 单一面、actuator=false、冲量分解与尺度 |
| RULED_OUT（代码/基准范围） | 标签进入力热、全基板刚化、封孔、预设自由轨迹、负热裁剪 | 未改几何/粒子标签通路；新旧被动标签基准PASS |
| BLOCKED | 稳定2us、同终时敏感性、完整材料自由能/边界一致性及论文生产验证 | 新候选RED、闭循环与常应力边界FAIL |

独立单粒子Coulomb参考进一步说明：在原第32步冻结载荷下，
切向外力小于摩擦承载能力；正确向量滑移先转向，约2.49830511885ns后停滑，
随后可粘着，而不是原候选继续错误滑移。
独立热Q=6.62994166314e-13J。
新局部解默认精度的热相对误差5.0141e-5、停滑时刻相对误差1.6182e-4；
更紧数值精度分别3.1816e-7、1.1162e-6。
参考由极坐标闭式积分推导，未调用被测contact函数，
文件first_particle_independent_reference.json和契约附录给出公式。
该精度检查不是调材料/摩擦参数，更不是拟合实验比值。

## 4. PHASE 2 独立回归

### A：真正一维反向滑移

m=1kg、u0=1m/s、Fother=-2N、muFn=1N、总时长1s。
测试中直接写入解析值：
tstop=1/3s；uf=-2/3m/s；Q=7/18J；Wother=1/9J；ΔK=-5/18J。

原算法测得uf=-2、Q=-0.5、Wother=1、ΔK=1.5，RED保留。
新局部变体准确分为停滑前/后两区间，粘着不可行时恢复相反方向滑移，
上述参考以绝对1e-10 SI容限PASS。没有把真实反向自动判失败。

新增多面+执行器独立例：x/z约束、y自由，两面各muFn=1，Fother_y=-4，
u0_y=1，dt=1；tstop=1/6，uf_y=-5/3，Q=14/9，Wother=22/9，ΔK=8/9。
接触迁移率非对角面块共享y响应=1，约束方向响应=0。PASS。

### B：完整弹性闭循环

实际调用material.update，只有测试law返回无限屈服来关闭塑性；
没有改production FKS/Jaumann。thermal/contact关闭。
F=I,S=0 → F=exp(eA) → .99^(1/3)exp(eA) → .99^(1/3)I → I，
A=diag(1,-1/2,-1/2)，e=.001，m=1kg，rho0=2700，G=26.3GPa。
四段各1s，分别16/32/64/128子步；V按J变化，记录全部子步F/S和sigma:D Vdt。

四档总功分别约0.14611111111119、0.14611111111135、
0.14611111111152、0.14611111111141 J。
独立精确偏功=(1-.99)*G*(A:A)*e²/rho0=0.146111111111111J。
压缩/恢复体积功相消到约1e-13J。
最终S=0，F误差≤2.33e-14；Es回到0，峰值14.611111111111J。
material没有求解内能变量；pa.e为未使用scratch，记录null/NOT_SOLVED，不冒充物理内能。

**解释：非零极限功来自当前Cauchy/Jaumann常G偏应力率与变体积功积分，
不是缩小dt即可消失的时间误差，也不是单独Es记账错误。**
反向走该闭路径会改变净功符号，不能把它当成可任意分配的耗散热。
回归PASS的含义是成功识别并独立验证这个缺陷；本构热力学验收仍FAIL。

## 5. 最小patch范围、接口和限制

| 文件/函数 | 改动与接口 | 状态 |
|---|---|---|
| kick_diagnostics.py / KickRecorder.record/accept/reject | 只读旧kick观察者；记录所有活跃面、前一kick、哈希、接受/被拒状态 | 已验证非干扰 |
| solver.py / __init__, evaluate, state | 可选diagnostic_recorder；contact_variant显式legacy/event_aware_v1；保留旧分支；状态记录变体 | 默认仍legacy；新变体BLOCKED |
| contact_kick.py / cone_solve | 共享约束后迁移率P/m、多面Coulomb圆盘；失败记录末迭代，不pinv/不放宽门槛 | 新近零反例FAIL |
| contact_kick.py / kick | 恒方向停滑分段/粘着可行性/恢复滑移；向量转向自适应中点；返回每子区间冲量、Q、反力功 | 指定独立基准PASS |
| contact_kick.py / domain_kick | 复用surfaces、原kn/法向罚力、原热分配与导热面积；返回加速度/热/账本/逐粒子局部轨迹 | 完整链短步PASS，2us试跑FAIL |
| output.py / restore/checkpoint | REJECTED拒绝恢复；先写pending再原子替换，硬中止不覆盖上个有效态 | 回归PASS |
| application.py / main | 显式--contact-variant；600s硬watchdog；FAILED保存被拒态、最后有效checkpoint不覆盖 | 实际输出完整性PASS |
| scripts/diagnose_first_coupling_failure.py | 限40kick/600s、observer逐位对比 | 原异常定位PASS |
| scripts/run_coupling_review.py | 统一101项回归；每个子命令≤60s；总exit1保留RED | 实测95/101 PASS |

新变体每个接受子区间严格使用Q=-Jt·R[(v_before+v_after)/2-vwall]；
圆盘投影约束摩擦承载力，不是截负热；Q原值均保存。
共享C_ab=R_a P R_b/m，多面不是各自孤立m|vt|/dt。
执行器反力及功显式保存；力/壁速度仍同位置。
最多4096局部区间、2000锥迭代，不收敛明确拒绝。
没有新增切向弹性储能、材料自由能、空间稳定化、人工反弹或能量残差热。

material.py、mechanics.py、geometry.py、contact.py、integrator.py、
thermal.py、particles.py、全部配置和原测试未改。
旧物理源码快照见logs/coupling_20260930/source_before/。
这些数值基准并未解决论文材料、边界、完整工具面等原缺项。

## 6. 测试矩阵与原记录保留

最终可重跑证据：logs/coupling_20260930/final_suite_v2/summary.json，各命令独立log。
汇总101项95 PASS/6 FAIL；脚本整体exit1。此前final_suite的100项94/6记录保留；
仅新增“checkpoint中断写入不损伤上个接受态”回归，无新增真实试跑。

| 项目 | 实际状态 |
|---|---|
| 新事件接触12项：反向、壁、共转、耦合、功/热、首发局部独立参考、REJECTED重启/原子checkpoint | 12 PASS |
| 实际材料闭循环四档，符合独立非零功参考 | 回归1 PASS；本构一致性FAIL |
| 仿射/自由表面畸变邻域、功/总内力/总力矩、条件数拒绝、单云零模/弱边界功 | 4 PASS |
| 非零初始应力刚体旋转（材料点+SPH链） | 原科学套件内PASS |
| 27粒子畸变自由表面云，81 DOF，零模SVD | 6个刚体零模、额外0；只限此云，不证明全局无零模 |
| 常应力强/节点边界patch | FAIL，独立牵引的相对残差0.326157187666 |
| 新变体标签/ID/排序/重启、legacy原快照/异常、变体切换拒绝 | 3 PASS |
| 原一维摩擦反向 | 1 FAIL，保留RED |
| 新真实候选454独立停滑参考 | 1 FAIL，保留RED |
| 原科学审查 | 35 PASS/3 FAIL（原raw两项与Es反例） |
| 原后端、原配置输出骨架 | 21 PASS、19 PASS |
| 失败试跑的3帧VTK/ID/hash/最后有效checkpoint | PASS，只证明输出完整性；manifest仍FAILED |
| DOP853独立检查 | NOT_RUN，隔离环境无SciPy；保留ImportError日志，改用闭式独立参考，无安装 |
| 新变体CFL .05/.025，dx .25mm | NOT_RUN，上档没达2us |
| 有限应变自由能新material-point实现、长时/生产三工况 | NOT_RUN |

首轮新空间条件数fixture将z压至1e-5，没有实际超过cond=1e6，
所以“应拒绝”断言失败；不是守卫失效。该核原点径向斜率不为零，
近重合竖向点导致最小特征值按dz而非dz²主导。
改测试fixture为1e-9后触发明确拒绝，容限/生产几何未改，原失败log保留。
常应力边界FAIL没有借弱仿射功PASS掩盖；本轮不修空间离散。

早期诊断入口的ValueError捕获和stdout NumPy序列化错误、测试fixture缺feed、
缺SciPy均保留对应attempt日志；物理结果只引用成功生成的phase1_verified
和final_suite，不把入口错误宣称为物理失败或测试通过。

## 7. 唯一新真实试跑及新BLOCKED事件

试跑前门控/预算：logs/coupling_20260930/preflight_gate.md。
固定物理配置configs/debug/starved.json没有改值：
dx=.0005m，CFL=.1，层厚3个粒距，N657；
正常初始dt≈1.04907e-8s，目标约191步。
试跑前首步实测0.8724357s，线性外推166.64s，硬上限600s。
后续单位回归重复测量另有preflight_timing.json，不替代首次预注册记录。

真实结果位于outputs/coupling_20260930/event_dx05_cfl01/：
status=FAILED，接受2步，物理时间2.0981304868235094e-8s，
wall=1.5001573s，已完成步median=0.51603895s，峰值工作集89,346,048bytes。
末态VTK/checkpoint均为第2步，不含被拒第3步。
REJECTED_candidate.npz及REJECTED_contact_debug.json明确隔离。
审计发现该run的辅助trace布尔标记写反；原JSON不覆盖，另存
TRACE_METADATA_ERRATUM.json：last_completed_contact_callback_trace是接受步2，
failed_local_context及其嵌套已完成局部trace才属候选3。时间字段可独立区分。
写入器已改为保守标记“不能保证当前候选”，不据此重跑/更改任何物理数据。
未启动任何第二次2us试跑，也没启动半/四分dt或细网格。

新候选：particle454，BORE_X_MINUS，dt=1.04905432172e-8s，
自由DOF全开，单面，内部切向力0.00744458184N、摩擦容量0.39425806322N。
独立闭式参考tstop=5.98547901766ns，Q=8.16487022434e-12J，停滑和粘着可行。
实际局部中点锥迭代在近零滑移阶段不收敛：
- 中点系数4.03127786413e-21s（对应子区间约8.06256e-21s）；
- 力更新两点振荡1.95263748579e-11N；
- 原接受阈值7.8851612644e-13N；
- 1997–2000次迭代逐项保存于new_candidate_RED.json；
- 基底包含0.10162575m/s法向量，而切向已约1e-14m/s。

这属于数值可解性/近零事件处理BLOCKED，不是允许忽略的物理负耗散，
更不能通过增加迭代数、放宽容限、硬设零速度或截热宣布通过。
下一最小工作仅限contact_kick.cone_solve/kick：
将面内未知量放到稳定正交切向基中、核查近零事件和尺度化KKT残差；
保持完整P/m与多面耦合，不偷换回独立m限幅。
先使这个独立RED和现有解析/功/耦合回归共同通过，再评估是否重跑同一固定问题。

## 8. 可重跑命令与阶段顺序

在项目根，使用现有.envs/pysph/python.exe；不安装系统依赖。
输出目录必须是新目录，不覆盖先前run。

    .envs/pysph/python.exe scripts/run_coupling_review.py --output logs/coupling_20260930/recheck_01

应返回exit1（6个真实RED），不是“全绿”。此命令不启动2us生产/试跑，
只含最多两步新变体及31+被拒一步legacy的回归。

    .envs/pysph/python.exe scripts/diagnose_first_coupling_failure.py --output outputs/coupling_20260930/legacy_recheck_01

最多40kick/600s；诊断捕获成功exit0，但report的物理状态仍PHYSICS_FAILED_DIAGNOSTIC_CAPTURED。
单独关键基准：
tests/review/test_contact_events.py；
tests/review/test_coupling_benchmarks.py --legacy-red 或 --cycle；
tests/review/test_new_candidate_red.py。

以下只记录本轮已执行命令，当前BLOCKED，不建议再次启动：

    .envs/pysph/python.exe scripts/run_afsd_debug.py --case configs/debug/starved.json --contact-variant event_aware_v1 --output outputs/coupling_20260930/event_dx05_cfl01 --end-time 0.000002 --steps 10000 --output-every 1 --wall-seconds 600

接触门槛修复后，仍先dx=.5mm/CFL=.1到2us，再.05/.025；
三档均成功达相同2us才比dt；之后才能dx=.25mm且物理边界/刚度尺度不变。
现在没有同终时收敛证据，不作生产成本/沉积形貌外推。
早先粗糙论文尺度性能外推不作为授权或本轮建议。

PHASE4仅给设计及独立材料点验收规格，见FINITE_STRAIN_VARIANT_DESIGN.md；
未实现/未接生产。当前材料、边界与来源缺项继续APPROXIMATE/BLOCKED。
旧SCIENTIFIC_FREEZE.json不重写；新授权修订的哈希清单单独保存，
文件完整性PASS绝不等于科学验收PASS。

