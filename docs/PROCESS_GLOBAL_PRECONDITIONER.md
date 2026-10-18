# 五场全局候选结构预条件器：事前合同

2026-10-03。状态：组件实施与独立测试待运行（NOT_RUN）；不授予物理步接受或生产。
本文件的推导与以下容限在首次测试前写入。旧 UpdatedTrial、stage_anderson、
材料、空间、接触、热学与物理接受标准均不修改。

## 所解问题与单位

真实五场候选 y=(x1,v1,rho1,T1,vbar)，真实残差 d=F(y)-y 来自原 UpdatedTrial。
vbar 是接触事件实际平均速度，不是恒等式 (v0+v1)/2。
这里只构造新的候选提议 y+P d；外层仍必须重新计算并用原容限验收 F(y)-y。
不得以 P d 的大小接受物理状态。

冻结 accepted 几何、真实 m/rho、三维 corrected/伴随弹性算子给 K [N/m]，
其物理质量 M [kg] 和时间 dt [s] 不变。A=M+dt^2 K/4 [kg]。
仅光滑无事件的冻结线弹性、端点中点消元近似给该 A；它不含 FKS/温度/
rho/空间重建/预应力/接触切换导数，不是完整 Jacobian。
体积刚度仅是近似压缩响应，并未声称消除了真实 rho 方程。

以 bar 指 vbar 字段、f/c 指自由/规定轴：

    delta_bar_c = d_bar_c
    A_ff delta_bar_f = (M d_bar)_f - A_fc delta_bar_c
    delta_x = d_x + dt (delta_bar - d_bar)
    delta_v1 = d_v1 + 2 (delta_bar - d_bar)
    delta_rho = d_rho
    delta_T = d_T

dx[m]、dv/dbar[m/s]、drho[kg/m3]、dT[K] 各自保留直接残差；
因此不会把已有运动学缺陷删除，或把两个独立速度字段强行合并。
后两式直接回代，没有冻结/削弱/裁剪热、密度真实残差。
dt 与系数 2 只描述光滑近似块的交叉修正；有事件时依旧保留 d_x 与 d_v1，
实际接触运动学只由后续 F(y) 再评估决定。

P 为可逆三角块：delta_bar=0 意味规定轴 dbar=0，自由块 M dbar=0，
进而 dbar=0；其余 delta 字段为0当且仅当其直接 d 为0。因此候选提议没有新增根。
这证明方程根集合不变，不证明非线性收敛或事件分支唯一性。

## 规定轴与诊断

由原 control(z0,t+dt/2) 取得原执行器自由轴，validate_interval 禁止跨 release。
只使用当前残差所要求的规定轴修正，不反复加入完整 feed 速度。
返回完整粒子ID/自由轴/质量、原五场d、五场delta、M dbar、约束值、
自由块移项rhs、A delta_bar-rhs 全轴“预条件反力”、自由线性残差。
该反力是代数预条件量，绝不是实际接触/执行器冲量，不进入热或外功账。

## 首次运行前冻结的验收条件

- 线性后向误差 eta_inf = ||A_ff delta_f-rhs_shifted||inf /
  (||A_ff||inf ||delta_f||inf + ||rhs_shifted||inf) <= 1e-12。
  空自由块或全零分子分母定义 eta=0。依据是双精度小型直接求解的后向稳定性；
  此标准不是完整非线性残差或物理误差容限，也不将 cond(A) 与空间 cond(Mi) 混用。
- 独立2DOF字面有理数参考、单位变换与五场交叉块：rtol=2e-12、atol=2e-13，
  数值均为 O(1) 制造量。SI粒子组件以逐字段单位缩放后的相对差 <=2e-11；
  空间刚模修正相对误差 <=2e-11。标签交换应逐位相等。
- 事件反例使用独立给定 u0=1、u1=-2/3、vbar=-1/18 m/s；不得强改成1/6。
  这是解析反向例的位移/时长，不是调用接触函数生成参考值。
- 输入和 accepted 粒子必须不变；0残差精确给0修正；非法shape/NaN/非正dt拒绝。
- 独立两粒子系统保留3轴，弹性仅联结其x轴；字面
  Mx=diag(2,3), Kx=[[4,-4],[-4,4]], dt=1，Ax=[[3,-1],[-1,4]]。
  dbar=(1,2)给delta=(14/11,20/11)，有规定第二轴dbar=0.2、第一轴dbar=1
  给delta=(11/15,1/5)，预条件反力第二轴=-8/15 Ns（自由轴0）。
- 测试单进程外部上限60s；最多两轮有证据实现修复，首FAIL留存。
  不运行 AFSD、真实新物理步或长时生产；真实残差/预算/提交由外层另立合同。

## 接口与解释边界

`StagePreconditioner(accepted_pa,p,dt,time)` 建构冻结 accepted 状态结构近似。
`.apply(candidate,mapping)` 只返回 proposal/correction/diagnostics，分类
`STRUCTURALLY_PRECONDITIONED_CANDIDATE_PROPOSAL_NOT_PHYSICAL_MAP`。
`StagePreconditioner.from_response(response, particle_ids)` 仅为独立组件制造输入。
返回量均独立copy；不修改粒子，不输出应力/热/塑性，不提交状态。

有限孔口、入口注入、隐式热尚未由本接口接入。
宏步精度/真实接触摩擦反向/宏观输运仍 NOT_RUN，不能由组件PASS自动继承。

## 首次 direct-block 运行结果（保留旧实现）

`tests/review/test_stage_preconditioner.py` 首次16/16 PASS，内0.069s、外0.405192s、exit0；
日志 `logs/process_global_20261003/preconditioner_tests_initial.log`。
无修复、无AFSD运行。事件测试只是独立给定真实反向路径平均值不被覆盖，
不是重测整个摩擦实现。rho直接残差虽不改根，但在较大dt的压缩耦合下
缺少连续性凝聚，因此下面注册改进候选预条件器；旧PASS不删除、不改名。

## 连续性凝聚变体：第二项实现前注册

新增 `ContinuityStagePreconditioner`，不修改已通过16项的direct类。
令 E=trace(D) [1/m]，V=accepted m/rho [m3]，rho_ref=p[rho0] [kg/m3]；
K_b=bulk [Pa]，K_bulk=E^T diag(V K_b) E [N/m]。
冻结于accepted构形、线性化rho约等于rho_ref的近似子块为

    (I+dt² M^-1 K_dev/4) delta_bar - C delta_rho = dbar
    dt diag(rho_ref) E delta_bar + delta_rho = drho
    C = dt/4 M^-1 E^T diag(V K_b/rho_ref)

压力闭合 d sigma/d rho = -K_b/rho_ref I 与伴随internal_force
给速度方程中的 +C rho 映射项，所以残差左侧为 -C。
消元rho得到原正定全弹性 A=M+dt²(K_dev+K_bulk)/4，且

    rhs = M dbar + dt/4 E^T (V K_b/rho_ref * drho)
    delta_bar_c=dbar_c
    A_ff delta_bar_f=rhs_f-A_fc delta_bar_c
    delta_rho=drho-dt rho_ref E delta_bar

delta_x、delta_v1依然用上文**原始**dbar的交叉修正，不用已增广rhs假装dbar。
delta_T=dT。此变体仍仅候选预条件，未改变真实中点有理连续性公式、FKS、
非均匀真实密度或接受标准；不把rho_ref替换写入粒子。
所有块单位相合：增广rhs为Ns，回代delta_rho为kg/m3。
输出完整原d/压缩rhs增量/密度回代项/规定轴/后向误差；代数反力不记物理账。

原后向误差<=1e-12及所有容限沿用，不为此变体放宽。
独立压缩参考：m=4kg,V=2m3,K_b=3Pa,rho_ref=5kg/m3,dt=2s,
E=(1,0,0)/m，c=.15、b=10，dbar=.5m/s、drho=2kg/m3。
手算delta_bar=(dbar+c drho)/(1+bc)=.32m/s，delta_rho=-1.2kg/m3。
参考常数必须字面写在测试，禁止调用被测预条件器生成。
另用实际27点corrected/伴随旧力路径构造冻结线性弹性与连续性map，
在无约束情况下比较预置解析仿射根并核原始完整map残差；不是AFSD。
新变体首次测试前仍NOT_RUN；仍每次外部60s、至多两轮证据修复。

### 连续性变体首次实测结果

`tests/review/test_stage_continuity_preconditioner.py` 首次10/10 PASS，内0.044s、
外0.3775748s、exit0，无修复；日志
`logs/process_global_20261003/preconditioner_continuity_tests_initial.log`。
实际27点冻结线性子块的原始五场map与预置仿射根相符，未调用被测E/K生成参考力。
累计26项组件测试，不等于26个真实过程或完整原接触算法测试。
真实全非线性UpdatedTrial、有限工具、注入、宏步误差与真实过孔仍由外层另验；本组件没有运行。

源码接口：`src/afsd_backend/pysph/stage_continuity_preconditioner.py`
的 `ContinuityStagePreconditioner(accepted_pa,p,dt,time).apply(candidate,mapping)`。
与direct类相同返回proposal/correction/diagnostics，另含pressure_density_rhs_Ns、
inertial_rhs_Ns、density_backsubstitution_kg_m3及reference-density近似标记。
完整命令见 `logs/process_global_20261003/preconditioner_COMMANDS.md`。
