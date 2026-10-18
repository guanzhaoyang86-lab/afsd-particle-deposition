# 整体求根器独立复核与过程时间推进决策

日期：2026-10-03。范围：`stage_structural_root.py`、`stage_preconditioner.py`、
`PROCESS_GLOBAL_ROOT_CONTRACT.md`、`PROCESS_GLOBAL_THERMAL_VARIANT.md`及对应root回归。
未改root/测试，未运行新的AFSD。材料/有限几何/热变体自身正确性不因本审查自动PASS。

## 审查结论

**允许按当前合同进行一次事前登记500ns的有界“求根诊断”，不允许把该步当过程准确性认证。**
合同明确一次、固定dt、原质量/模量/力学容限、完整物理接受、真实费用probe和<=600s上限，
并明确500ns不解析最高弹性模态，范围是合理的。无须为这个诊断先解决全部宏步物理；
但其PASS/FAIL后必须围绕已证障碍转向过程方法，不能延伸为更多同尺度短跑即“完成”。

## 源码逐项审查

| 检查 | 状态 | 证据与边界 |
|---|---|---|
| 实际五场raw收敛 | PASS，源码+制造回归 | 每次trial.evaluate之后直接residual_norms；scaled<=1才返回，未用preconditioned norm代替 |
| 预条件器仅改未接受candidate | PASS，源码+制造回归 | apply结果仅用于迭代f；last物理map不被proposal替换；数值诊断另置numerical_proposal_only |
| 极小/零proposal假PASS | PASS，独立参考测试 | 原速度map=-2u+.002、root=1/1500；proposal=0或1e-20仍按raw超预算FAIL |
| 最佳锚点及增长保护 | PASS，源码+测试 | >2*best raw仅拒绝加速candidate；恢复best的x/f并清历史，保留原dt |
| 独立event-average未知量 | PASS，制造回归 | literal mean=-.004不等于endpoint-average；未由预条件器强制约束成两端均值 |
| 原史/输入不污染 | PASS，root接口范围 | initial_candidate拷贝，root无commit；制造shared初数组逐字段不变 |
| 故障现场对应当前candidate | PASS，源码+测试 | 每循环last=None；evaluate异常无旧eval冒充；proposal异常则保留本次真实eval与失败row |
| dt与evaluate预算 | PASS，测试范围 | dt原样传每次trial；maximum_evaluations精确；墙钟前后检查；外硬预算仍由driver保障 |
| NaN与域异常 | PASS，测试范围 | 不用proposal替换坏map，保存异常上下文并退出；无秘密line search/缩dt |
| 真AFSD500ns根 | NOT_RUN，本审查时 | 必须以随后真实manifest判定，不能用10个制造测试代替 |
| 宏步时间误差与材料输运 | NOT_RUN/BLOCKED | 本root并不提供这些认证 |

`last['numerical_proposal_only']`会增加输出量，应计入same-source cost probe；
root的时间检查不能中断一个正在进行的trial，因此内/外硬deadline仍必要。
root接口本身不核preconditioner的dt/ID/source是否与trial匹配，driver必须在构造/冻结时绑定；
即使错误预条件在raw过门前通常只影响求根效率，也不能依赖该性质免除来源检查。

本次审查初版StagePreconditioner为rho/T直接defect，未做rho Schur RHS/回代，
其局限与独立2x2反例已详见PROCESS_GLOBAL_STEP_AUDIT。
随后若加入新Schur类，需要该类独立回归与接口/源hash重绑；当前root测试不能自动认证新类。

## 独立复跑证据

先读原首跑日志，再独立执行同一最终测试文件；没有改测试/放宽容限。

```text
.envs/pysph/python.exe -B tests/review/test_stage_structural_root.py
```

外层subprocess timeout=60s，OPENBLAS/MKL/OMP线程均1。
`logs/process_global_20261003/root_tests_independent_review.log`：10/10 PASS，
测试0.051s，外0.2329721s，exit0。是同10个测试的复核，不增加为20项不同验证。
原首跑`root_tests_initial.log`10/10 PASS保留。

独立参考包括：

- `u=-2u+.002` => u=1/1500 m/s，不调用求根器生成参考。
- 五场线性根v=1e-3*(221/146)，x=.002+1e-5*(152/73)，rho=2700.27，T=723.02，mean=-.004。
- m=1kg、k=4N/m、h=1s，x0=.002m、v0=.001m/s的中点振子解：
  x1=.0005m、v1=-.004m/s、mean=-.0015m/s；使用实际结构proposal，但力和参考为literal独立式。

这些验证根/代数，不是热态FKS、旋转多面接触或宏步准确性验证。

## 500ns：如何报告准确性，不能怎样报告

前次只读证据`outputs/process_global_step_audit_20261003/scales.json`给出冻结弹性
omega_max*h≈10.0278。标量中点振子的每步数值相位为2atan(theta/2)，与真实theta之比≈.274。
这不是“实际AFSD整体误差72.6%”：当前状态在该最高模态的幅度和非线性模态分配没有测，
不能把单模态相位比冒充真实位移/应力误差。它足以否定“隐式root PASS即动态准确”的推理。

若500ns返回原完整接受，需要单独输出：

1. raw五场根精度、动量、event均速、原始Q、热账、feature与域等接受证据。
2. normal_temporal_residual、angular_residual及额外event-average几何项、
   total_energy_diagnostic_residual；注明旧Es并非自由能，不能补成热。
3. 新T为BE热变体，不与旧显式热根混称同方程；参数不变不等于时间离散没变。
4. 状态为`NUMERICAL_SOLVABILITY_DIAGNOSTIC; TEMPORAL_ACCURACY_NOT_VALIDATED`，
   不凭该接受态赋予生产/持续输运许可。

若driver仅写snapshot，应有明确published acceptance索引才允许后续重启；
“可重载NPZ”与“可作为新物理解合法接续”仍是不同证据。

## 接下来不是再堆50ns/500ns，而是一个过程积分器的有限选择

当前真实初态处于非零预压/应力/动能释放，不能直接删惯性或清零速度。
另一方面，隐式中点对保守高频振子不衰减幅度；把h抬到毫秒只是给高频严重相位错位，
不等于已进入可信的慢过程。真正的过程计算需要同时选择：
**保留哪些瞬态、容许何种显式登记的数值耗散、怎样在同一区间解析接触/材料/热。**

### L-stable惯性变体值得评估，但不是一句“改BE”

可以设计独立命名的真实惯性隐式过程变体，保留真实M/G/FKS和实际工艺时间；
若采用L-stable方法，应将算法耗散与塑性耗散、摩擦热、外功分开存储。
这是时间离散选择，不是质量/时间缩放，也不自动授权用数值耗散抹除有意义的预压瞬态。
需在合同中冻结适用目标，并用慢加载敏感性证明过滤高频没有改变所关心的料流/反力/P/Q。

独立标量说明：对y'=i*omega*y，BE放大率1/(1-i*theta)，其模为1/sqrt(1+theta²)。
theta≈10时每步幅度约.099，而中点幅度为1；这是**算法滤除**，不是物理阻尼被发现了。
对固定线弹性无外功系统，BE离散能量损失可独立推导为

```text
D_alg = .5*delta_v^T M delta_v + .5*delta_x^T K delta_x
E1-E0 = -D_alg
```

该恒等式仅适用于固定线性M/K与一致BE运动学，不能直接把本项目Jaumann/FKS总能量差
都命名为D_alg，更不能将该差加进温度。非线性模型无法闭合的剩余项必须另列未解释残差。

**不能直接用BE端点摩擦功替换用户要求的事件区间Q=-Jt·u_mid。**
接触滑移反向、stick可行性和多面约束必须仍有一致的实际子区间速度/冲量；
简单先做BE内部力再套旧整步冻结kick，会再次失去整体相容性。
需要联合残差与其结构/接触切线，而不是把P/m替换为W后宣称完成。
更高阶L-stable方案也需明确其多stage材料史/接触事件语义，不能机械套RK系数。

### 一个可控工作包，而非无期限局部修补

下一工作包建议只完成一个**单变体联合过程步**及一组关键独立基准：

| 同时必须覆盖 | 独立证据/停止条件 |
|---|---|
| 惯性与时间滤波 | 快慢双频弹性振子；慢受迫响应解析解；报告高频算法损失，不把它当热；同终时步细化 |
| 接触区间 | 原1D反向reference u=-2/3、Q=7/18；加弹性耦合/多面执行器后全J/W，不仅Q>=0 |
| 几何事件 | 真有限孔口/肩缘跨越时的法向势能及工具功；不能只证明50um域内旧face从未切换 |
| 本构/热 | 独立非共轴非零应力更新和BE热源闭合；实际率/T耦合，同史不累积 |
| 原始ID与供料 | 只在接受事务入口插点，完整外部质量/动量/焓账，passive标签不入方程 |

原root在一次500ns结果中若失败，先由记录判断是rho压力反馈、mapped域、接触预算或其它；
只做该证据支持的有限修复，不反复换dt找到绿色样本。
若root成功，也不应立刻再做大量近邻短窗；以本次成本和上述联合基准确定
能否进入一次预算内、同时含可观察输运量的过程探针。

过程h最终应由允许的几何/接触/材料/热截断误差与成本共同决定，而非只由声学CFL或root稳定性。
工具1ms转1.8度只是几何量级，不是自动许可1ms步；当前50um安全域仍无法过孔，
必须在有限接触/事件覆盖后重订几何可达域，不能悄悄删除guard。

如果目标必须保留全部启动高频瞬态且不可负担解析其时间尺度，则应如实报告成本阻塞；
如果要省略启动/采用准静态，需独立新初始化与授权，不能从旧动态checkpoint伪接续。
本审查没有改选这类物理模型，也没有开启新长期/付费/生产任务。

## 达到的层级

新root：**代码审查+10项独立参考回归复核PASS**。
500ns：**合同允许一次求解诊断，实际结果另见manifest；本审查未跑**。
过程可用L-stable/联合事件变体：**仅设计建议，NOT_IMPLEMENTED/NOT_RUN**。
持续塑性输运、过孔/沉积、过程步敏感性、论文定量验证：**未被以上证据证明**。
