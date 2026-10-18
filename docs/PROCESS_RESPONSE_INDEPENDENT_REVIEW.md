# 实际342点结构响应探针的独立只读审查

2026-10-03。审查对象：scripts/probe_process_response.py、PROCESS_RESPONSE_PROBE_CONTRACT.md、
outputs/process_completion_20261003/response_probe/report.json及同名日志。
本审查未重跑矩阵/物理、未修改已哈希源码/合同/结果；不把读取既有PASS当成新增测试。
结论：**现有报告对“冻结弹性组件”的限定正确，可进入全局候选预条件接口开发，不能据此放大真实dt或生产。**

## 源码逐项支持的结论

| 项目 | 审查判断 | 证据与边界 |
|---|---|---|
| 真实接受态来源 | SUPPORTED_BY_CHECKS | 脚本核原manifest来源/产物hash、accepted标记、完整particle_state_hash；最后再核输入/源码不漂移。本文未另重算全部hash |
| 真3D算子 | CONFIRMED_SOURCE_PATH | build_frozen_elastic调用既有PySPH三维NNPS/corrected B，342点/8131pair，D为3078×1026；不是输出坐标复制 |
| 自由轴与质量 | CONFIRMED_SOURCE_PATH | current control(z0,t)只驱动25个z轴；1001自由/25规定轴。m原样重复三轴，没有因病态额外固定自由DOF |
| 接触映射同点 | CONFIRMED_SOURCE_PATH | 当前活跃法向填入对应粒心的3轴；104行，没有换用面投影点速度或施力位置 |
| 异粒子响应 | SUPPORTED_MEASURED_COMPONENT | W跨不同particle_id元素最大非零；这只是冻结线弹性响应，不是已测真实材料接触传递力 |
| 实际弹性力相合 | PASS_RECORDED_SINGLE_DIRECTION | 独立旧corrected_operators+字面弹性stress+internal_force路径，力差4.1210178e-16N，原判据通过 |
| 全耦合宏步/摩擦/热 | NOT_RUN | 脚本不调用material.update、contact kick、thermal或积分advance，new_coupled_steps=0 |

27点组件14回归提供独立手算、能量、刚模、有限差分方向导数等先验支持；
但本342点探针本身**只额外核一次冻结线弹性方向力差**，不是又在实际构形运行了一整套能量/客观率/接触测试。
合同中的“力/功相合”可由该伴随代数路径理解，报告没有独立记录实际时段的功增量，不能把它写成真实过程能量PASS。

## 数值结果如何解读

| 量 | 50ns | 1ms |
|---|---:|---:|
| A自由块condition | 1.2513913 | 100556507.9301 |
| 构造/谱检查耗时 | 0.597622s | 0.593669s |
| W求解耗时 | 0.025629s | 0.024945s |
| max跨粒子W / kg^-1 | 27463.5604 | 13194.0851 |
| max逐元素W-W^T / kg^-1 | 2.91038e-11 | 1.36514e-6 |

共同D/K装配0.1594011s；全部内1.4774996s、外1.8175799s；峰值RSS164294656B。
总成本包含旧证据hash及两个冻结矩阵响应，**不是1ms耦合求解的单步成本**。
它只消除“342点冻结D/K装配必然不可负担”的猜测，不证明全局非线性迭代次数、事件预算或工艺运行成本。

### 风险1：两种条件数不能混用

- 空间`M_i=sum_j V_j r_ij tensor grad W`是每粒子3×3一致性矩阵，原`cond(M_i)>1e6`继续硬失败。
- 此处A_ff是1001×1001、以kg计的质量-弹性线性响应块，`A=M+dt²K/4`。
  其1e8不是空间邻域退化证据，也不是违反原M门槛；不能因此pinv、增加夹具或改真实质量。
- A的较大condition表明真实质量控制的软模与较硬弹性模相距较远；会放大某些方向的求解误差。
  `cond(A)*machine_epsilon`约2.2e-8仅是舍入敏感性尺度，不是此次实际误差或新的接受阈值。
  下一接口应记录自由块线性后向误差
  `||A_ff dv_f-rhs_shifted||/(||A_ff|| ||dv_f||+||rhs_shifted||)`，以及真实非线性残差。
  不能仅凭Cholesky成功或某一个condition数判全局求解通过。

### 风险2：零模、W对称性与缺失证据

K在事前1e-10相对cutoff下有6个近零值，与刚体模式相容；该计数不提供每个模态的物理识别，
也不能排除其他构形/自由表面存在零能模。A的质量项正定不会把这些问题变成“材料已稳定”。
原K最小值-1.2319093e-8N/m未裁剪，组件相对尺度检查通过；本文不另改判负刚度模型。

1ms W不对称绝对差比50ns大，不能只看1.36e-6的绝对数断言失败。
报告未保存W全矩阵、完整对称归一尺度、特征谱或线性后向误差，
因此本审查不能独立重算W精度或证明其满足后续摩擦锥所需精度。
现有报告没有承诺这些，所以不为追绿改原结果；下一真正使用W的联合测试需同时保存并检查这些量。
不能把两个dt下某个最大交叉元素的下降解释为实际材料耦合变弱，两者既非同一个完整本构切线也非物理力。

### 风险3：104行只包含法向

这里Gmap没有每面的两个切向基、摩擦锥、法向罚势导数、接触激活/开离或旋转法向导数。
也没有材料FKS切线、预应力和几何刚度、rho/T耦合导数。
因此不能把现成104×104 W直接交给旧粒子摩擦kick并宣布“effective mass已补齐”。
若未来构造完整接触Schur块，法/切行必须共用同一结构响应，并保留非零执行器规定运动的仿射项。
W只给增量响应；它不是实际接触速度、接触力、normal spring储能或外功本身。

## 最小合理global接口：先改求根路径，不先改物理残差

建议下一新接口明确命名`STRUCTURALLY_PRECONDITIONED_CANDIDATE_ROOT`，
复用当前同史UpdatedTrial和接受标准；不把1ms当其首次真实运行dt。
这是验证结构响应可用于全局候选修正的最小单元，不声称仅该单元就得到完整过程求解器。

```text
accepted state + candidate five fields
  -> unchanged actual residual + branch/event journal
  -> frozen elastic structural preconditioner (candidate correction only)
  -> bounded globalized Newton/Krylov or nonlinear-correction proposal
  -> re-evaluate unchanged actual residual/domain/physics checks
  -> ROOT_CONVERGED_NOT_ACCEPTED, then existing full acceptance/atomic commit
```

接口应隔离为三件事：

1. `evaluate_candidate(accepted, y, dt)`：仍给原五场缺陷、材料/热/接触原始账、feature和branch签名；
   每次同史，无Newton累计。预条件模型不得生成应力/温度/热或作为最终接受态。
2. `build_preconditioner(pa, dt, free_axes)`：只生成此次冻结弹性近似块，可从accepted或指定candidate重建；
   重建位置必须在日志声明。约束自由轴来自真实执行器，不从K特征向量删除自由度。
3. `apply_preconditioner(residual, variable_convention)`：单位和端点/平均速度约定明确。
   保存rhs、修正量、规定轴修正和线性残差；外层必须重新评估真实五场残差，不能以预条件残差代替。

### 不能弄错vbar未知量

当前UpdatedTrial的`mean_velocity_m_s`是实际接触事件平均速度，
并不普遍等于`(v0+v1)/2`。候选v1对力反馈也不等同于vbar，不能按字段名硬套普通中点推导。

在**冻结光滑线弹性、无事件**的独立基准中，
`r_bar=vbar-mapped_vbar`的近似机械块可以写成
`A delta_vbar = -M r_bar`，由`A=M+dt²K/4`给预条件修正；
相应`delta_x=dt*delta_vbar`、`delta_v1=2*delta_vbar`仅是该光滑近似块的关系。
体积Kbulk项还隐含线性化连续方程消元；若保留rho块，需要明确块消元/回代约定。
**这些关系只能服务预条件器，不能覆盖跨停滑/反向事件的真实运动学或真实五场残差。**
对旧UpdatedTrial的完整导数，材料率、温度、rho、NNPS/B及接触事件变化都还需包含。
不能把弹性A直接标成完整Jacobian。

对同一个接触冲量增量J，组件`contact_response(..., at='endpoint')`给端点dv响应，
普通中点才是W/2；真实事件平均响应必须由实际事件路径线性化/分段确定，不能无条件使用1/2。
规定轴输入也应是当前未知量的**修正量**，不是每次重新灌入完整送料速度。

## 下一次验证的最小范围与停点

先用少量联合独立参考验证新global接口，而非马上跑过程时长：

- 2DOF真实惯性弹簧＋一受约束接触，独立消元解析根；比较预条件前后同一残差方程根一致。
- 非零初速度的滑移反向例保留uf=-2/3、Q=7/18；证明预条件器未把vbar强改成端点平均。
- 实际小3D current-geometry trial：多个方向的完整残差差分，只把弹性A标近似，
  测根迭代/evaluate成本、旧接受态不变、标签/ID交换不影响力学。
- 非线性根PASS、线性后向误差PASS、物理接受PASS、时间离散精度PASS分别报告。
  热隐式组件、有限工具、入口即使各自通过，也不能自动继承为本global接口的联合PASS。

以上新global接口/联合测试均NOT_RUN。先冻结变量/单位/限次及总预算，再实施一个有界组件；
不能以探针1.8s预测全耦合1ms步，更不能先按这个数启动0.1s或三工况。
此审查未发现要求回改本次已哈希探针的错误；发现的是其允许解释和下游使用的明确边界。
项目真实物理时间仍850.393633ns，无新增沉积证据。
