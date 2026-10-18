# 更新构形整体候选：独立接口审查（2026-10-03）

审查对象为新增 `src/afsd_backend/pysph/updated_trial.py` 及其真实调用链。
未修改旧空间、本构、接触、热学或冻结配置；此审查没有运行新的真实AFSD时间步。
用户允许工程近似不等于授予75%：该文件只记录代码/独立基准/剩余风险，
全局实际运行状态以本轮运行报告为准，不自动继承此处单测PASS。

## 实际执行顺序与时间层

`accepted PA独立克隆 → 候选端点域检查 → xmid/rhomid/Tmid → 新3D NNPS/M/B/L(vbar)
→ 同一S0/ep0的原FKS/Jaumann/J2 → 中点压力/应力及伴随内力
→ 从真实v0作旧event-aware多面/轴向约束kick → before/mid/after feature守卫
→ 原热学分账 → 事件平均位移/速度与连续方程map → map域检查 → NOT_ACCEPTED输出`。

- 候选中点位置为 `(x0+x1)/2`，密度和FKS温度分别为端点算术中点。
- L采用独立未知量 `mean_velocity_m_s`，不是端点速度平均；接触结果以各接受子区间
  `Σdt_a(v_before_a+v_after_a)/2dt` 返回真正事件平均速度，驱动 `x1_map=x0+dt*vbar_map`。
- 每次材料更新从克隆的同一 `S0/ep0` 开始；塑性功乘同一中点体积 `m/rhomid`。
- 内力使用同一中点几何、密度、B和g，与梯度保持原伴随对。
- 热增量从T0起，内部成对导热用T0，FKS及接触导热反馈Tmid；原真实导热行和时间步界仍检查。
  return耗散、法向势能残差、机械残差均未被加入热源。
- 原始z0只选择轴向执行器；tracer/AS/RS及缓存actuator标签不参与受力或热学。
- 该接口不提交ParticleArray/累计账/checkpoint。`accepted=False`、`restart_allowed=False`。

## 独立联合基准结果

唯一新增测试文件：`tests/review/test_updated_trial_integration.py`。
最终6项PASS，外进程0.7674334s、测试内部0.426s；外60s保护。
证据：`logs/updated_trial_20261003/integration_tests_final.log`。

| 项目 | 独立参考与证据 | 状态 |
|---|---|---|
| 当前三维邻域/核矩阵 | 27真实PySPH点；测试独立全有序点对径向导数求M/B，不调用被测kernel/gradient；改变候选几何和密度后B与原构形有可测差异 | PASS |
| 三维仿射与功/净力矩 | 独立矩阵A，L=A误差<=1e-10/s；同构形Σv·F+ΣVσ:A及净内力/力矩按原尺度容限 | PASS |
| 连续方程 | 各向同性aI，独立 `rho1/rho0=(1-1.5a dt)/(1+1.5a dt)`；非任意有限步指数解 | PASS |
| 非零初应力刚转 | 独立sin/cos的QS0QT；坐标旋转角θ的中点梯度对应材料角2tan(θ/2)，两者有限步误差不混淆 | PASS |
| 标签/重排/纯候选 | 实际342点、192板+150棒；换tracer/wholebar/ASRS/actuator缓存并重排，候选物理结果一致；重复/拒绝不改accepted属性 | PASS |
| 反向摩擦事件 | 实际旧boundary/有限底座，m1/u0=1/Fother=-2/muFn1/dt1；vbar=-1/18、uf=-2/3、Q=7/18、Wother=1/9、ΔK=-5/18；端点假热=-1/18不能代替 | PASS |
| 真实checkpoint来源 | `short_07us/plastic_accepted.npz`，step62、650.3936329ns、342点、epmax1.8228662e-5；独立核artifact和33源SHA及配置，导入保留非零v/S/ep/T | PASS：仅导入，未推进 |

最后两行属于独立文件的第5、6个联合测试，并非额外计数。
独立参考不调用被测contact/material/kernel来生成预期值。
制造仿射/旋转候选是接口探针，不是给真实自由材料规定运动。

### 初次FAIL完整保留

`integration_tests_initial.log`：5PASS+1FAIL，外1.3000336s。
刚转fixture先使用 `(x1-x0)/dt` 构造平均速度；近等坐标相减生成约1e-11/s量级的虚假
对称梯度，经26.3GPa真实模量和dt放大，最大应力差3.42823341e-6Pa，超过5e-7Pa容限。
仅修测试构造为独立Cayley skew矩阵作用于xmid；另验 `dt*vbar=x1-x0` 至3e-19m。
未改solver、本构、物性、dt或5e-7Pa容限。最终同一6项PASS。

## 接受层必须补查的接口风险

1. **candidate与map不可混作一个接受状态。** 当前feature守卫针对被evaluate的候选端点；
   map另作域/穿透检查，但没有独立完整feature切换守卫。
   若全局层最终保存map，须对实际保存端点重查feature，不能因为位置残差小就忽略
   非连续feature掩膜。接受的S/ep/热账也必须来自同一次最终同史evaluate。
2. **冻结位置仍是近似。** 事件平均速度修正位移/接触功，但 `(x0+x1)/2` 一般不等于
   非线性事件轨迹的时间平均位置。法向力/壁速仍冻结该构形，不是精确大步接触积分。
   法向势能时间残差必须原样记录，不得补热。
3. **角动量账不能无条件继承旧漂移法。** 对事件平均位移，严格端点恒等式为
   `ΔL=m*xmid×Δv + m*Δx×(v0+v1)/2`。第二项在普通端点中点漂移时为零，
   在反向/多面事件后通常不为零；旧boundary仅记录冻结xmid×外冲量。
   全局层应报告实际端点角动量残差及上述运动学项，或按事件位置独立积分力矩。
   不得虚构外矩来补账，不可仅凭空间内部净力矩PASS宣称整步角动量精确守恒。
4. **求根精度不是时间精度。** 有界Picard可检查此变体的非线性一致性，
   不证明任意大dt、沉积过程或dt/空间收敛。全局根尚须实际运行与独立读回。
5. **热学仍受显式部分的真实时间步界。** Tmid材料/接触反馈不是全隐式导热；
   不能因为机械全局迭代收敛就关闭旧导热稳定检查。
6. **既有科学局限保留。** Jaumann/Cauchy闭循环自由能缺陷、有限棱边、热启动/初始干涉预载、
   板棒预结合及有限料棒/持续入口均未被本候选接口解决。

## 未授予的验证

本独立审查阶段，全局真实50ns试步、多步过程推进、同终时dt/空间敏感性、
持续供料/肩面承载/沉积输运、三工况及论文定量对照均未由本测试文件执行：NOT_RUN。
不能把6个联合测试等同于完整复现，不能把冻结哈希等同于物理正确性。

## 全局runner接受接口只读审查

随后已只读审查 `scripts/probe_updated_stage.py` 与 `stage_iteration.py`，没有运行真实probe/root。

- `solve_stage` 对五个未知场求同一个map固定点；.75只混合未提交候选。
  最大32次evaluate、明确时间限、无自动dt改写/再启动；根收敛先返回NOT_ACCEPTED。
- runner保持最后被evaluate的candidate作为保存状态，不把它snap成未经evaluate的map；
  因此上一节candidate/map feature风险在当前保存接口中已经规避。
- 保存S/ep/D及应变率来自同次最后e；接受检查用保存c.T重算mcΔT而非只看thermal内部map温度。
  FKS必须是原参数默认实例，不能将制造law输出转为真实接受态。
- 检查实际step62来源：旧冻结manifest、全产物、33源、配置与checkpoint都核验；提交前再次核哈希。
  状态保留原始非零应力、塑性和速度，没有清零惯性后转静力。
- checkpoint名为新变体snapshot，显式 `restart_allowed=False`，不冒称已验证常规续算；
  `valid_properties` 排除继承的旧scratch/未更新累计数组，真实逐粒子P/Q等在e记录中另存。
- 成本probe也保存完整未提交候选；仅同源码、配置及产物hash的实测耗时可放行一次有界root。

**运行前发现的确定接口错误：** 接受函数曾读取normal snapshot的`potential_J`，
而真实返回键是 `normal_potential_J`。该错误会在根已收敛后抛KeyError导致拒绝。
已通知主代理修正精确键；不是物理参数或容限调整，尚未启动真实物理来追错。
最终修复和测试状态以主报告/后续日志为准。

**需显式输出的差异：** `boundary.ledger.delta_kinetic_J`对应边界map速度端点，
保存候选c.v仅在1e-9m/s残差容限内与之相同。候选dK与该账的差应原样记录，
不要将边界内部功身份式PASS表述为保存状态精确同式成立。
同理，候选mean与接触事件mean之间的内力功差、以及前述角动量运动学项需单独记录。

## 唯一真实global attempt的独立只读结果

实际运行完成后，新增 `scripts/assess_updated_stage.py`，只读取已落盘证据，
不调用UpdatedTrial、acceptance、contact、material、mechanics或thermal求解函数。
首次读回完成，无修复重跑：外2.5160054s，内部1.5214866s，外60s保护。
输出 `outputs/updated_trial_20261003/independent_readback.json`；
日志 `logs/updated_trial_20261003/independent_readback_initial.log`。

**状态是PASS_BLOCKED_EVIDENCE_ONLY_NOT_ACCEPTED，不是物理试步PASS。**
原global attempt在固定32次evaluate后因迭代限REJECTED；没有新接受checkpoint。
所有38个登记源码和原始输入/产物hash前后不变，旧接受史仍是step62来源；
本次失败不覆盖上一变体已存在的700ns轨迹。

独立读回重建最后一组完整未提交candidate/evaluation的8131点对、M/B/L、空间功/净力矩、
7566接触子区间、原始摩擦Q、时间平均速度/位移、冲量、塑性增量工作、径向导热及实际温增账。
它不把未收敛候选状态作为物理结果，不执行第二次求根。

| 检查量 | 最后完整未提交候选的证据 |
|---|---|
| 子区间动量残差最大值 | 9.51258e-24 N·s |
| 空间功共轭残差 | -2.22045e-16 W |
| corrected矩阵条件数最大值 | 2.2150576 |
| raw摩擦热总量 | 5.24192596e-9 J |
| 最小raw摩擦热 | -1.08311933e-25 J，处于独立复算的浮点roundoff bound内，未裁剪 |
| 候选温增−真实热源总差 | 3.39072344e-10 J；因候选未收敛，不授予热平衡接受 |
| 候选动能−contact map动能差 | -3.65381713e-12 J，原样记录不补热 |
| 候选mean−event mean内力功差 | 7.03049521e-12 J，原样记录不补热 |

### 主导残差不是当前接触face异常

| 变量 | 最大残差 | 粒子/方向 | 当前contact features |
|---|---|---|---|
| 端点速度 | 5.0427445671e-4 m/s | PID268，棒，z | 无 |
| 区间平均速度 | 2.5213722836e-4 m/s | PID268，棒，z | 无 |
| 位移 | 1.26068613e-11 m | PID268，棒，z | 无 |
| 密度 | 1.43589341e-5 kg/m3，相对5.31814539e-9 | PID270，棒 | 无 |
| 温度 | 1.116958856e-6 K | PID269，棒 | 无 |

只读32行iterations显示：最大归一化残差由首轮9.04754e7降到第8轮2.37362e5，
其后回升/震荡，第16–32轮停在约5.04e5，不是只差最后一次接触kick。

- CONFIRMED：现有.75混合Picard在该固定50ns真实系统上未达到冻结容限；接触事件本身完成。
- SUPPORTED：下一主要问题是整体自由棒材料/空间/热耦合映射的非线性求解，
  不应继续将主导失败归给一个接触face。
  主导粒子没有直接接触并不排除邻域传力及全局接触耦合影响。
- POSSIBLE：塑性切换/全局切线或映射收缩性导致振荡；本轮没有Jacobian/谱证据，不能定论。
- NOT_RUN：接受分支snapshot读回、额外dt、第二次global attempt、多步/同终时敏感性、生产。

下一工作应研究同一物理系统的有界整体非线性求解/切线，保留本次REJECTED与独立证据；
不能通过放宽容限、加热补差、清零自由材料惯性或隐藏重复试跑把本轮改判成功。
