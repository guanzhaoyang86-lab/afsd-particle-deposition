# 1–4关整体集成：事前冻结合同

2026-10-02，用户明确要求整体推进旋转孔/送料、可变形基板、塑性与热学，不再逐一孤立交付。
分类ENGINEERING_APPROXIMATE_COUPLED_AFSD_LOCAL_DEBUG，仍PySPH固定pin；不重搭项目或优化GPU/并行。
本阶段采用唯一集成适配入口，复用现有SolidSolver顺序及其粒子、3D NNPS、积分、本构、导热和输出设施。
旧源码、参数、FAIL/警告不覆盖；新增命名adapter，不偷偷将新约束变成旧生产默认。

## 已有证据与本轮要完成的边界

旧continuation_e3/to_2us已实际跑191步、657粒子、塑性/热同时非零：3个自由棒粒子塑化，
ep最大1.69159e-5、塑性功6.0569e-9J；没有基板塑性或送料面穿越，不能称持续沉积。
Q7的0.2s是27点冻结参考构形/制造J2准静态组件、无热无旋转，不能改名作为现成长AFSD求解器。
这次从新初态启动，不能绕过旧checkpoint来源校验。

目标是四机制在**同一个**可重跑3D实例中接通，真实记录有无塑性、热和材料传力。
“代码接通”“运行接受”“发生塑性/传热”“持续宏观沉积”分别判定，最后一项不预先授予。
若短窗仍没有可持续塑性或未到终时，明确PARTIAL/BLOCKED；不通过人为应力、屈服调小或自由初速度制造结果。

## 固定物理与几何

configs/benchmarks/coupled_afsd.json继承configs/debug/starved.json及base的SI/来源。
变更只有显式DEBUG几何和边界：基板4×4×1.5mm、方棒2.5×2.5×3mm；dx=.5mm，h=.65mm；
总342点=192基板+150棒，均可变形，无材料底层夹具。原rho/G/nu及接触刚度尺度保持。
工具半径2.7mm、开放旋转方孔宽2.496mm（原每面2um干涉）、肩面z=1.5mm、工具顶z=2.5mm。
底座z=-1.5mm刚性。层高1.5mm有3dx；缩小横向几何/棒长只属DEBUG，不与旧657点冒称网格收敛。
棒/板共享SPH邻域，仍预结合近似，不预测接合生成/脱粘。工具外缘/孔口棱边最近点接触未完整实现。

omega=300rpm=31.4159265359rad/s，走速2.12mm/s，轴向进给2.12mm/s。
只有原始顶0.5mm棒层的z受控，xy自由；初态其vz=-feed，其余速度0。不同于旧锁x/z路径，已显式命名。
控制由原始z0选材料层，不用tracer/ASRS/可改actuator标签。不预置自由材料转动、径向流或横移。
release_time=1s只表示注册微秒窗内保持驱动，不声称复现启动、停留、抬升/终止程序。
初始T=723.15K、工具733.15K、底座703.15K是既有热启动假设，不是假装已模拟的预热。
FKS b1无量纲解释、1/s低率floor、常数k167W/(mK)、cp900J/(kgK)、mu.3、接触h1000W/(m²K)、
摩擦入材料比例.5、塑性热比例chi1均原近似；单位/出处承继PARAMETER_SOURCES，未新增论文事实。
标记带2–3mm为局部诊断，非论文3–4mm带；wholebar标签独立。实验室+z轴CCW/+x走速，AS为y<0。

## 真实调用与新适配接口

半漂移 → 3D NNPS/corrected L → 旧Jaumann/J2/FKS update(Told)
→ 中点应力/adjoint内力 → 实际有限face contact_rows → 每粒子event-aware kick+z-only约束
→ raw径向核成对导热+实际接触导热+塑性/摩擦热 → 旧temperature_step
→ 速度/密度kick → 半漂移 → 验收/原子commit。
不是整体隐式热力求解；是保留的显式一阶材料/热分裂，仍受实际稳定dt限制。
机械校正梯度与伴随力成对，导热仍用原始径向g以保持正导热系数及成对守恒。

coupled_boundary.evaluate_boundary复用contact_rows/kick；所有力和壁速取同一个粒子中心作用点。
face投影点仅诊断，未重构面速度或假装有伴随分配；这是已声明中心同点近似。
Fn=kn*overlap，U逐实际feature求和；before/mid/after及feature on/off均记录。
空间/时刻切换丢失非零势能时不能伪装热；未支持孔口/外缘储能feature切换需停止并报告。
kick共享约束mobility，反力不能落在自由xy；无contact粒子也有内力/约束完整记录。

coupled_thermal.advance_thermal只接收物理源：
Qplastic=chi*P，Qfr_material=partition*Qfr_raw，Qcontact=dt*hA*(Twall-T)，成对导热内部净0。
boundary.contact_power_W严格只含接触导热，不预加摩擦；摩擦只在thermal加一次。
未热化塑性功(1-chi)P单列；return_dissipation、能量残差不进入热。
Qwall_out=(1-partition)*Qfr_raw-Qcontact（可有符号），
Ediagnostic=K+Es+U+Eth+(1-chi)P；R=E-E0-Wtool-Wact+Qwall_out。
原Es仍非认证自由能；其残差仅WARNING，不能借此让NaN/温域/穿透失效继续。
初始K/P/L/U/热参考明确，不能按“初态静止且无储能”计账。
接受状态、REJECTED候选、source/config hash分存；失败恢复全部PA字段及ledger/time/step，拒绝态不可restart。

## 事前数值容限与独立检验

- 字面反力/功单点atol1e-12SI；FKS/J2解析参考不调用被测函数生成。
- 热制造参考：m=.0027kg、cp900、T400K、G30GPa、Y90MPa，dev应变diag(.002,-.001,-.001)。
  qtrial180MPa、dep.001、Snew=diag(60,-30,-30)MPa、P.09J、return.045J。
  chi=.9后Qplastic=.081J、非热化.009J、Tnew400.0333333333333K；制造Y不进入AFSD配置。
- 纯导热手算系数1152/35W/K、两点温差100K，功率±115200/35W；检查从热向冷、总和0。
- 混合热源账误差1e-12J+1e-9*sum(abs(source energies))；超限FAIL，不用舍入热修补。
- 标签置换/清零不改变物理；排序/重启x/v/T/rho atol1e-12，S atol1e-6Pa；ID唯一。
- 质量相对1e-12；密度相对变化<=.05；穿透<=.1dx；speed<=20m/s；全部finite；T在[room,melt)。
- 最大位移50um为有限feature短窗域限制；cond(M)<=1e6，不pinv。
- dt同时受声学CFL=.1、法向.1sqrt(m/kn)、实际热行和.25mc/row限制。
- 空间功1e-12W+1e-10*sum|F||v|；净力1e-12N+1e-10sum|F|；净矩1e-14Nm+1e-10sum|r×F|。
- kick功账1e-15J+1e-9*sumabsworks；自由轴反力1e-17Ns。原始负摩擦热按旧事件算法处理，不裁零。
- 总能量诊断沿原5%交换/初能尺度+1e-7J作WARNING；不声明热力学完整，不给能量误差补热。
- 源/config/输出hash严格一致；失败原记录保留。没有为符合预期曲线的调参循环。

## 成本与唯一实际集成尝试

短测试/独立读回外部60s；5步成本探针外180s；完整试跑内部540s、外570s，均<=600s。
冻结可选终时2或4us。成本探针后按1.5*最慢实测步时*预估步数+30s选择不超过450s的最长窗，
选择只基于成本，不看塑化比例或预期图；若均超预算则完整运行NOT_RUN，不拆段绕预算。
只有一次冻结几何/参数的正式尝试，失败记录完整；不反复修改CFL/接触/材料追绿。
完成后从接受序列独立统计棒/基板自由响应、塑化粒子/连续发生塑性步数、真实温变、热账、
接触features、质量/密度/穿透/能量及ID，判断1–4关究竟达到什么层级。
单个完整实例达到终时不证明收敛；本轮三工况、时间/空间敏感性、长期持续沉积仍NOT_RUN。
没有GPU/并行优化、质量/时间/模量缩放、自由粒子动画或论文误差百分比。
