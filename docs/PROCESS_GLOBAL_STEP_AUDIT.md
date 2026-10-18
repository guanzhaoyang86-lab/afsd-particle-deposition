# 同史整体步独立审查：结构预条件与一次更大步诊断

日期：2026-10-03。范围：实际源码与旧stage004归档；一次纯读回/标量推导，
**未运行AFSD、未调用UpdatedTrial/material/contact生成新的物理candidate、未改旧文件**。
本审查不是“下一步已通过”。新方向的真实根、接触事件成本、同终时误差全部NOT_RUN。

## 结论与建议

建议预先登记 **dt=500ns** 的一次有界整体求根诊断，不先选5us，不以失败后减步追绿。
选择依据是实际物性/构形的尺度差，而不是猜哪档会PASS：500ns已经跨过显式声学限制，
足以测试结构预条件；但尚小于独立罚接触振子周期的1/8左右。5us则跨过一个罚接触周期，
且旧内力孤立外推末速23.35m/s超过原20m/s域守卫。这个外推不是新map，不授予5us实测FAIL。

即使500ns根通过，**也仅获得原离散方程更大步的根，不证明启动瞬态或过程步准确**。
冻结弹性最高模态的中点相位比约0.274，已经提示显著的高频欠解析。
不应再把几百个微秒根堆成“接近结束”；需以一次结果决定联合残差/接触几何的下一有限修复。

## 实际调用顺序及限制

```text
stage_anderson.solve_stage
  initial_candidate: x1=x0+dt*v0；其余四场继承
  UpdatedTrial.evaluate
    validate_interval(t,dt,release)
    五场shape/finite；candidate _check_state(x1,v1,rho1,T1,vbar)
    xmid/rhomid/Tmid -> NNPS(dim=3) -> corrected B -> L(vbar)
    material.update(S0,ep0,L,Tmid,dt) [每次同一原史]
    sigma=.5(S0+S1)-p(rhomid)I -> paired internal_force
    evaluate_boundary(xmid,v0,...,fint)
      原z0顶层控制 -> 每粒contact_kick.kick [多面同P/m]
      事件区间 -> 实际v1与积分平均vbar、带符号Q、反力/功
    before->mid / mid->candidate feature guards
    advance_thermal(T0显式pair；Tmid接触；本次P/Q)
    连续性有理式 -> 五场mapped
    mapped _check_state [这里抛错则raw residual尚未返回]
    返回candidate-minus-map及全物理journal
  同原五场容限判断；整体acceptance；最后才commit
```

源码：`src/afsd_backend/pysph/updated_trial.py`、`stage_anderson.py`、
`coupled_boundary.py`、`contact_kick.py`、`coupled_thermal.py`、
`axial_actuator.py`、`coupled_solver.py:_feature_transition_guard`。
这里没有以函数名假设先后；以上顺序来自evaluate实际调用行。

| 守卫/预算 | 当前真实实参和实施位置 | 影响 |
|---|---|---|
| 最大速度 | 20m/s；_check_state同时检v1及独立vbar | candidate和mapped都检；并非只约束最终接受态 |
| 原始位移 | 50um；距粒子x0/y0/z0的总范数 | 不是距当前接受态的本步增量 |
| 密度 | rho>0且相对2700kg/m3不超过5% | candidate、mapped及连续性分母/分子另须正 |
| 温域 | 293.15K<=T<855.15K | FKS与热各自仍检域，不裁温/屈服 |
| 穿透 | 0.1dx=50um | normal_snapshot/contact_rows；深埋工具实体另拒绝 |
| 空间 | corrected moment正定、cond<=1e6 | 无pinv；与结构矩阵A条件数不是同一个门槛 |
| 热步 | dt<=.25min(mc/(pair_row+wall_row)) | stage004=292.157883us；不是500ns/5us首要标量限制 |
| feature | 双半段检活动ID消失/出现且延拓仍有正罚储能 | 真正unsupported finite-edge transition仍拒绝，不能用位移小掩盖 |
| 执行器 | 原z0的最上.5mm，仅z规定-feed；release=1s | 不允许跨release；x/y自由，25轴约束 |
| 接触局部 | 每粒max_intervals=4096；rtol=2e-7、atol=1e-12m/s | 4096对attempt/records分别检查；不是全局7746记录的上限 |
| 原Anderson | 每阶段32次、β=.75、depth6、QRcond<=1e10 | 任一trial域异常直接StageFailure；无域内线搜索重试 |

当前有限几何为底座、平肩及四孔壁，未接新finite_tool_geometry。
法/切力和壁速均在粒心同点；face projection仅诊断。
单粒kick内部原内力/法向力/法向和壁速冻结，滑移事件会分段但不更新这些空间量。

## 已运行的独立尺度诊断

脚本：`scripts/assess_global_step_scales.py`。只用上游NPZ读取器、JSON/gzip和NumPy数组算术，
没有import被测trial/contact/material/preconditioner用于参考。旧输入与manifest SHA核验。
命令（单线程环境下，外进程timeout=60s）：

```text
.envs/pysph/python.exe -B scripts/assess_global_step_scales.py outputs/process_global_step_audit_20261003/scales.json
```

真实日志：`logs/process_global_step_audit_20261003/scales.log`。
首跑exit0；脚本0.160309s、外0.484556s；新增物理步0。
报告保存输入hash和脚本hash，不覆写旧归档。

读取stage004：t=850.393633ns、342点、vmax=.349929263m/s；
最大原始位移=.138646040um，穿透1.972596247um，最大密度偏差5.5494e-5；
T范围723.149986530–723.150068597K；旧区间内力/m最大4.669777e6m/s2。
该旧区间共有7746接触子区间，单粒最多550（PID341），不是尚余3546个全局名额。

独立材料尺度：G=26.3GPa、nu=.33、rho=2700给K=68.586274510GPa，
cL=sqrt((K+4G/3)/rho)=6195.964672m/s；dx/cL=80.697684ns，h/cL=104.906989ns，
旧CFL=.1的h/cL尺度=10.490699ns（只作独立尺度，不代替实际stepper控制）。
m=3.375e-7kg、kn=657500N/m给omega_contact=sqrt(kn/m)=1.395761e6rad/s、周期4.501620us。

| 指标 | 500ns | 5us |
|---|---:|---:|
| 工具转角 | .0009度 | .009度 |
| 轴向命令位移 | 1.06nm | 10.6nm |
| x0+dt*v0弹道预测最大原始位移 | .192859um | 1.763995um |
| omega_contact*dt | .697880 | 6.978804 |
| 旧力单独v0+dt*fint/m的max速 | 2.336482m/s | 23.350070m/s |
| dt/旧热界 | .0017114 | .0171140 |
| 冻结弹性omega_max*dt | 10.027787 | 100.277867 |

**上述孤立力外推省略接触/约束、重算材料与中点几何，因此不是actual-map域检查结果。**
它只提示5us的初始非线性残差可能在守卫外，不证明真根在域外。
转角小也不意味着罚接触或材料动态准确。

最高弹性频率来自旧`response_probe/report.json`在50ns的A自由块条件数1.25139126655，
真实质量均匀且保留自由x/y刚体平移，因此lambda_min(A)=m，
omega_max=sqrt(4(cond(A)-1))/dt=2.005557345e7rad/s，周期313.288739ns。
这是未预应力冻结弹性K的频率推断，不是已测完整非线性Jacobian谱。
线性中点每步相位为2atan(omega*dt/2)，与实际omega*dt之比在500ns=.274025、5us=.030931。
因此root PASS也必须与准确性/能量限制分开。

## mapped-domain与求根域：CONFIRMED设计耦合，第一真失败仍NOT_RUN

原`evaluate`在返回残差前要求mapped也满足全部物理域；数学上这是要求c属于D且F(c)属于D。
Newton只需在可定义残差的candidate域D求R(c)=c-F(c)=0；F(c)本身在离根时并非接受状态。
一个域内candidate因F(c)速度大而报错，不等价于D内无根。
当前Anderson更会将此首异常当整个attempt终止，无法借线搜索回到可评估点。

最小后续方向（若一次真实journal证实这正是阻塞）：

1. 保留旧UpdatedTrial不动，单独residual-evaluation policy区分candidate域与mapped诊断域；
   只对**未投入物理求值的mapped位置/速度门**考虑返回原值+违规证据，不clip。
2. candidate密度/温域、constitutive合法性、真接触负耗散、NNPS/feature/thermal稳定性仍必须失败；
   不把这些物理求值异常统一当“Newton可穿越域”。
3. raw残差照原tol；最终candidate **和map** 全部回走原domain检查后才accept。
4. 独立标量manufactured反例需证明：D内有根但F(c0)在D外、线搜索不伪造状态；
   再测真根在D外时必须拒绝。所有reject/line-search candidate留档。

这不是要求现在预先放松门槛；确切首触发器必须以真实attempt为证。只读诊断未定位新失败。

## 结构预条件：rho凝聚、独立vbar与完整raw接受

审查到新增`stage_preconditioner.py`当前使用A=M+dt²K/4，
delta_vbar=A^-1 M defect_vbar，再修正位置/末速提议；rho/T仍用原map直接defect。
K来自接受构形的真实D与伴随力，保留真实M/G/K。
此操作不生成新的物理内力，也不得将接触kick内P/m直接改成A^-1。

### 独立推导和限制

无接触、冻结空间、线弹性、rho约等于参考值时：
S_mid=S0+G dt dev(Dvbar)，连续性rho1/rho0=(1-a)/(1+a)，a=dt div(vbar)/2。
其严格导数d(rho1)/d(div)=-rho0*dt/(1+a)^2。
在a=0附近消元rho，sigma_mid对vbar的变化对应(dt/2)C D delta_vbar，
因此mean-velocity残差的弹性块是M+dt²K/4。
这仅解释该预条件器尺度，不是完整FKS/Jaumann/变密度/几何/接触切线。
真实rho消元也不能忽略D随rhomid体积和B改变；独立vbar更不能强制为(v0+v1)/2。

当前实现未把rho残差带入速度Schur RHS，也未按新delta_vbar回代rho。
一个不用被测函数的标量压缩反例为F_bar=c*r、F_r=-b*bar（b,c>0，s=bc）：

```text
dbar=c*r-bar; dr=-b*bar-r
精确线性Newton：delta_bar=(dbar+c*dr)/(1+s)
                delta_r=dr-b*delta_bar
当前块预条件： delta_bar=dbar/(1+s); delta_r=dr
```

β=1当前迭代矩阵为[[s/(1+s),c/(1+s)],[-b,0]]，谱半径sqrt(s/(1+s))->1。
这说明大dt时仅把bulk K加进A却保留rho直接提议，**不保证**已经清除连续性/压力的刚性反馈。
可以保留它为受原raw验收的可逆近似提议，但不能称完整凝聚Newton。
建议将此2×2参考加入测试，并决定有限一次加入一致Schur RHS/回代，或明确等同缺失切线。
本审查只推导，不修改他人正在实现的接口。

预条件变换在A自由块可逆、约束轴correction保存时是块三角可逆；
所以它不会数学上创造原始defect不为0而预条件defect为0的根。
但接受必须检查**实际重新evaluate的原五场**，不能只用变换后的norm或一个装配残差。
现有`probe_updated_stage.acceptance`原tol、ID/原史、功共轭、反力/动量、实际事件均速、
raw摩擦/热账足以防这种代数假PASS，前提是source/classification和真实journal无替换。

### 原acceptance不是宏步准确性认证

原acceptance只**返回**total_energy_diagnostic_residual、normal_temporal_residual、
angular_residual与额外event-average角动量项；没有用p中的energy_relative_limit/.absolute_limit拒绝它们。
这是原短步分支已命名NOT_ENERGY_CERTIFIED，不是本轮才引入的隐藏改变。
因此“原完整raw接受”足够区分真根/假根，却不足以在500ns授予
宏步守恒/精度/自由能认证。新探针必须至少逐项显式报告这些量和NOT_VALIDATED状态，
不把同史根收敛当时步敏感性结果；不得用残差补热。

## 有限下一步与状态

- PASS（只读核验）：旧stage004来源/实参、5个输入hash、尺度计算/脚本首跑。
- CONFIRMED：双层domain检查可能阻止返回原残差；热不是500ns/5us的标量首限；当前K不是全切线。
- SUPPORTED：500ns是区分全局结构刚性与域guard的合理单次诊断尺度，5us暂不作为首跑。
- NOT_RUN：500ns实际首失败/根、完整Jacobian、rhoSchur独立回归、宏步同终时敏感性。
- BLOCKED：真正料流仍受50um总位移域/未集成有限棱边与入口限制；新根PASS不能消除这些。

一个真实attempt仍须事前独立回归和成本probe，按实测candidate耗时设置总evaluate预算，
内540/外570s或更紧，超过600s不得启动。不以多轮更换dt/参数继续追绿。
本审查没有许可去惯性、隐藏数值耗散、缩放物性、扩大接触事件预算或删除历史FAIL。
