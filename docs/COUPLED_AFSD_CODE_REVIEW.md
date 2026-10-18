# 集成 AFSD 适配器：独立源码与证据审查

2026-10-02；范围为 `coupled_solver.py`、`coupled_boundary.py`、`coupled_thermal.py`、
其实际旧调用路径与已落盘测试日志。没有修改冻结运行源码/合同，没有启动新物理算例。
本审查不重新认定论文参数，不将代码作者或上一报告的 PASS 当作宏观验证。
审查依据是源码调用、字面能量关系及当前日志；本次独立重跑物理/材料闭循环均 **NOT_RUN**。
短窗运行时开始本审查；完成后只读核对了最终 manifest、history 与外墙钟日志，见第6节。
没有代替独立全轨迹 reader 宣布物理账验证通过。

## 1. 实际调用链，而非按函数名猜测

```text
run_coupled_short_window.main
  └─ 原 run_coupled_afsd.run（只绑定显式短窗 settings/source/预算）
      └─ IntegratedSolver(initial_particles(p))，从 t=0 新建状态
          └─ advance：stable_dt + release 边界检查 + 保存接受态
              └─ HostLeapfrog / 固定 PySPH LeapFrogIntegrator
                  ├─ stage1：半漂移
                  ├─ evaluate(t+dt/2, dt)
                  │   ├─ 3D NNPS → corrected L、drho、B、原始径向 g
                  │   ├─ material.update(Sold, epold, L, Told)：Jaumann/J2/FKS
                  │   ├─ sigma_mid = (Sold+Snew)/2 - p(rho_mid) I
                  │   ├─ internal_force(sigma_mid, B)：伴随内力
                  │   ├─ evaluate_boundary：实际有限面 → 多面 kick + z-only 约束
                  │   │   └─ 同时记录接触/约束冲量、相对中点功、壁功
                  │   ├─ 用实际 kick 区间中点速度计算 Lmid、空间功身份式
                  │   └─ advance_thermal：P、原始 Qfr、接触导热、内部成对导热
                  └─ stage2：速度/密度 kick + 后半漂移
              └─ 检查 before→mid→after feature、有限值/域/功/热 → atomic commit
                 任一步失败：保留 REJECTED 候选，恢复全部 PA/账本/计数
```

源码定位：`coupled_solver.py:141,255`，`integrator.py`，
`third_party/PySPH/pysph/sph/integrator.py::LeapFrogIntegrator.one_timestep`。
这里 material 在 contact 之前；constraints 在同一个 contact kick 内，而不是之后覆盖自由速度。
继承旧 `SolidSolver` 的初始化和 `stable_dt`，但**不调用旧 evaluate/domain_kick 的 x,z 双轴夹持路径**。
旧 `particles.make_array` 给过顶层 x 初速度，新入口 `initial_particles` 显式清零全部 u/v/w，
再仅给原始顶层25点 `vz=-0.00212 m/s`；初始动能和动量计入账。

## 2. 四关目前可以授予的层级

| 关卡 | 代码实现 | 独立/组件基准 | 本审查可直接确认的真实轨迹 | 不授予的宏观结论 |
|---|---|---|---|---|
| 1 旋转开放孔 + 轴向送料边界 | PASS，有限面、中心同点近似；顶层仅 z 受控 | PASS：边界14项，含字面反力功、双面反向、共享 mobility、开孔、同点共转 | 67步/700ns完成；壁与轴向执行器外功非零；进料面穿越仍0 | 连续入口、实际棒材完整输送、孔口棱边/外缘完整接触、持续沉积 NOT_RUN/BLOCKED |
| 2 可变形基板与棒 | PASS，192基板+150棒，真实3D邻域和9分量应力；无底层材料夹具 | PASS：独立3D/仿射/非零应力旋转/功与净力矩；集成12项包含自由响应与原子恢复 | 真实自由材料响应，末态基板最大位移约15.77nm，不是二维复制 | 大变形可靠性、自由表面定量精度、同几何空间收敛 NOT_RUN/未解决 |
| 3 塑性流动 | PASS：原 Jaumann/J2/FKS 进入同一调用链，未换成制造 Y | PASS：独立J2/FKS公式；制造 Y=90MPa 的热塑分账 | 第62步首次非零ep；末态2个自由棒点塑化、基板0点；允许称“局部塑性响应” | 不等于持续塑性流动、板面塑化或材料沉积 |
| 4 热耦合 | PASS：导热、塑性功、摩擦热、接触导热同一实例 | PASS：热12项；独立系数、符号、单次分配、dt界、拒绝裁剪 | 同一短窗摩擦/接触/塑性热均非零，真实微小温变；全轨迹账另交独立reader核验 | 预热过程、稳态温场、实验温度吻合、塑性热主导机制 NOT_RUN |

“代码 PASS”限定为本表声明的 APPROXIMATE 范围，不是完整论文关卡通过。
旧657点2us轨迹中少量塑化粒子的证据不能冒领给新342点初态。

## 3. 关键科学/边界审查

### PASS（有明确范围）

- **三维与自由度**：`mechanics.pairs` 为 `LinkedListNNPS(dim=3)`；修正矩阵和应力均3×3；
  梯度 `Lraw @ B` 与 `internal_force(..., correction=B)` 成对；病态矩阵显式失败，无 pinv。
  `coupled_boundary` 对所有材料点调用同一力学 kick；`axial_actuator.control(z0,...)`
  仅原始顶层 z 不自由，全部 xy 自由。tracer/ASRS/可修改 actuator 标签不选择物理自由度。
- **开放方孔**：`geometry.surfaces` 的肩面排除方孔内部，四个孔壁随工具旋转；
  `contact_rows` 调 `assert_outside_tool`，不存在实心圆柱填孔或穿工具后投影修补。
- **功共轭接触**：接触力施加点与用于壁功的刚体速度均取粒子中心；投影面点仅诊断。
  多面共享 `P/m`，而非各面独立限幅；约束反力和功显式保存。自由轴反力有验收。
- **原始热量不裁剪**：`Qfr=-Jt·u_mid` 按接受子区间累积；舍入界只影响拒绝判断，
  不执行 abs/max 改热量。温度和非正屈服直接拒绝，不改写为上下限。
- **不隐性缩放**：粒子 `m=rho0*dx^3`；沿用 `rho0=2700 kg/m³`、`G=26.3 GPa`；
  法向刚度来源 `penalty_factor*G*dx²/contact_reference_spacing` 明示；短窗仅改观察终时。
  声学/法向/热学稳定条件仍在，未放大dt或物理时间。
- **标签与重启**：可变标签不进入材料/力/热；集成测试交换标签后物理逐位相同，
  重排按ID回读，非零应力/温度 checkpoint 续步相同；REJECTED、源或参数不符拒绝重启。

### 保留风险：不能通过本轮 PASS 消除

1. **材料自由能一致性：FAIL 保留。** `coupled_solver.energy` 仍用
   `Es=sum(V*S:S/(4G))` 加 EOS 体积能。独立关系已足以说明限制：纯各向同性压缩时
   `dev(D)=0`，无旋转/塑化则 S 不变且 `S:D=0`，但 V 变化使 Es 改变。
   所以此 Es 不是该 Jaumann/Cauchy 更新的有限变形功势。
   旧闭弹性路径的非零功0.146111111111J及细化不消失证据保持，见
   `FINITE_STRAIN_VARIANT_DESIGN.md`、`COUPLING_FAILURE_REVIEW.md`；本轮未重跑。
   更小全局残差或 `DIAGNOSTIC_WITHIN_SCREEN` 不能反证这一已知结构问题。
2. **时间层/分裂：APPROXIMATE。** 材料与连续方程用 kick 前 L；接触改变速度后
   Lmid 只用于真实区间功诊断，并未回代重求应力。因此空间伴随功身份式通过，
   不等于整个有限步应力—能量积分精确。日志分开记录 material_old_work、midpoint_stress_work，
   热塑采用 Told。当前没有同终时dt收敛证据，不能称二阶全耦合。
3. **正常储能时间层：APPROXIMATE。** 中点几何的 Fn/法向/壁速在整个 kick 冻结；
   `Rn=Uafter-Ubefore+Wnormal_relative` 单列，不加热。
   有储能的有限面切换由 before→mid、mid→after 守卫拒绝，
   但这不是完整孔口/肩外缘最近点接触，也不认证任意大位移路径。
4. **预结合与初始预压：ASSUMPTION。** 棒板共享邻域，相当于预结合，不预测新接合/脱粘；
   初始每孔面2µm干涉产生约5.26e-5J法向储能，初期自由速度可以来自它的释放。
   “棒运动起来了”不等于轴向送料产生了稳态搅拌；必须同时报告初始储能、壁功/反力与进料面穿越。
5. **材料与热参数：APPROXIMATE/BLOCKED。** FKS 保留无量纲 b1 解释及 `1/s` 低率正则化。
   `max(rate,rate_floor)` 是已显式登记的率模型，不是没有使用任何 floor；其效果须保留计数。
   热软化项确为 `1-(1+exp(-c1*T*))^(-1/c2)`。率floor、b1、常数k/cp、热分配、热启动
   不因代数测试通过就成为经论文或实验唯一确认的参数。熔化及材料温域外均不支持。
6. **空间精度：未解决。** corrected 功/矩和仿射测试不能替代受载边界、自由表面、
   畸变点云/零能模态及定量网格收敛。原 raw Table-1 两项 FAIL 仍实测保留，未改回 raw。
7. **步长执行是 fail-closed，不是自动热步自适应。** 继承 stable_dt 给声学/法向及热预估；
   `advance_thermal` 另外以真实内部+壁面行和检查 `dt<=.25*m*cp/row`，超界拒绝并恢复，
   并未自动缩子步重跑。本工况声学dt远小于热限制，但不能因此声称任意热配置都自动选稳步。

## 4. 热/总能量账核对

令 P 为实际 J2 返回的 `V*qnew*dep`；return 项 `V*1.5*G*dep²` 只诊断。

```text
材料热增量 = chi*P + partition*Qfr_raw + dt*Qcontact_W + dt*Qinternal_W
未热化塑性功 = (1-chi)*P
流向刚体热库 = (1-partition)*Qfr_raw - dt*Qcontact_W
内部成对导热总和 = 0（舍入尺度）
总诊断 R = (K+Es+U+Eth+未热化塑性功) - E0 - Wwall - Wact + Qwall_out
```

`boundary.contact_power_W` **仅**为带符号的接触导热；摩擦只在 thermal 分配一次。
法向相对功、法向时间层残差、numerical_return_dissipation、总残差均不进入热源。
初始动能/接触储能/热参考明确扣除。PySPH 继承的 `e` 不被伪称已求解的认证总内能，
输出真实温度、应力、ep 等；VTK 没有虚构完整自由能。
热源账通过只认证给定源→温增的一致性，不能认证未知材料/摩擦源本身与实验准确一致。

## 5. 证据路径与本审查的使用边界

- `logs/coupled_afsd_20261002/boundary_tests_initial.log`：14 PASS。
- `thermal_tests_initial.log`：12 PASS；独立制造参考 Y90MPa、P=.09J、
  chi=.9 时热=.081J、非热化=.009J、T=400+1/30K；不是 AFSD 参数。
- `solver_tests_source_final.log`：12 PASS；真实小步、标签、ID、热分账、原子失败、checkpoint。
- `old_mechanics_review.log`：14 PASS；`old_gate_b.log`：5 PASS；
  `old_contact_events.log`：12 PASS；`old_contact_review.log`：10 PASS。
- `old_open_risks.log`：2 FAIL（raw 边界旋转伪应变、非等密度功不共轭），原记录保留。
- `short_window_tests_initial.log`：10 PASS；仅预算、参数/证据不变和原入口绑定，不运行物理。
- `outputs/coupled_afsd_20261002/cost_probe/`：5接受步至52.4527456ns；
  此范围温度/摩擦/接触热非零，塑性ep/P为0；不能与后来的正式短窗混报。

即使短窗 manifest 最终 COMPLETE，也只表示登记终时内的数值接受及可重放证据；
独立 reader 的 PASS 进一步证明保存字段/身份式一致，**不是**独立求解器交叉验证、
论文定量准确度、完整几何接触、持续沉积、2/4µs达标、时间/空间收敛或三工况生产授权。
建议后续根据一次短窗的真实机制证据选择最小过程缺口，不继续为局部比例追绿；
任何更长时窗、换材料/边界或连续供料方案须另行显式登记，不能由可视化偷偷改物理。

## 6. 已完成短窗的只读终态补核

证据：`outputs/coupled_afsd_20261002/short_07us/{manifest.json,history.jsonl}`，
`logs/coupled_afsd_20261002/short_07us.log`。本节核对保存值及含义，非独立重新求解。

| 量 | 实际保存值 |
|---|---|
| 状态 / 步数 / 终时 | COMPLETE_REGISTERED_INTERVAL_NOT_SUSTAINED_DEPOSITION；67；7e-7s |
| 内 / 外墙钟 | 131.6388125s / 132.2734607s，低于事前431.764629s成本估计 |
| 首次塑性 | 第62步，650.3936329ns，1点；后续末态2点 |
| 末态塑性粒子 | 自由棒ID266/269，ep分别1.822866198757164e-5 / 1.45797157667538e-5；基板0点 |
| 塑性功 / 塑性热 | 各8.829178011658406e-9J（本配置chi=1）；不能另把return加热 |
| 入材料摩擦热 / 接触导热 | 2.2785158386994034e-7J / -1.4041484426862866e-7J |
| 温度范围 | 723.1499889980487–723.1500631928834K；初态723.15K |
| 基板最大位移 / 进料面穿越质量 | 1.5773176142441807e-8m / 0kg |
| 全局能量诊断残差 | 4.1161720927001887e-7J |
| 初始总能量 / 初始法向势能 | 5.260001896075826e-5J / 5.260000000000826e-5J |
| 壁功 / 执行器功 | 4.4277447180574865e-8J / 1.2426020399430737e-8J |
| 法向时间层累计残差 / numerical return | -1.96602797750244e-10J / 2.68683500392314e-9J；均不补热 |

必须保留解释：全局残差为初始能量的0.78254%，保存状态因此在既有screen内；
但同一残差约为本窗壁功加执行器功的7.26倍，**不能**写成“能量误差不到1%，因此准确”。
初始预压储能支配该分母；材料Es非自由能和显式时间层不一致的风险仍在。
该比值是已保存账本的尺度比较，不是论文误差百分比，也不是新的事后验收阈值。
热源账残差约-2.84e-16J只说明热增量分账一致，不能抵消这个机械/本构风险。

当前允许的结论是：四机制已在一个可复跑真实3D短轨迹共同起作用，且出现少量局部塑性。
不授予持续塑性输运、宏观沉积或三工况生产资格；完整2/4µs和敏感性检查仍NOT_RUN。
