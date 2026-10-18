# 冻结构形弹性结构响应组件（事前合同，2026-10-03）

范围：保留PySPH三维NNPS、corrected B及伴随内力，只新建当前构形弹性预条件/响应组件。
不是完整FKS、预应力/几何/接触切线，不接入生产，不改变任何旧源码/物理/失败。
真实质量、剪切模量、体积模量和物理dt由调用者给定，不增刚、不缩放、不删自由DOF、不pinv。

## 独立公式与接口

从实际当前邻域构建线性梯度D，`L_i=D_i u`；
`sigma=2G symdev(L)+Kbulk tr(L) I`，`V=m/rho`，
`K=D^T diag(V C_elastic) D`，内力`f=-K u`。
组装`A=M+dt² K/4`；Z只是明确规定自由轴的选择矩阵，不能为了可解临时增加约束。
对冲量rhs与规定速度增量dc，求`Aff dvf=rhsf-Afc dc`；反力冲量为`A dv-rhs`。
接触行Gmap的端点冲量响应`W=Gmap Z Aff^-1 Z^T Gmap^T`；中点响应明确为W/2。
该W只含当前冻结线弹性，不是完整非线性接触effective mass认证。

计划接口：`build_frozen_elastic(pa,shear,bulk)`、`ConstrainedResponse(mass_per_dof,K,dt,free)`，
`response(impulse,prescribed_delta_velocity=None)`、`contact_response(Gmap,at='endpoint'|'midpoint')`。
所有数据是脱离输入的副本；不更新ParticleArray材料史或位置，不生成物理checkpoint。

## 测试前冻结容限与独立reference

- 两DOF手算：m=(2,3)kg、K=[[4,-4],[-4,4]]N/m、dt=1s；
  A=[[3,-1],[-1,4]]kg，A^-1=(1/11)[[4,1],[1,3]]kg^-1。
  Gmap=I时跨接触响应1/11不是零；中点是上述矩阵的一半。
- 相同系统约束第二DOF dv2=.2m/s，rhs=(1,0)N·s，则dv=(.4,.2)m/s，reaction=(0,.4)N·s。
- 质量-only：K=0，W=Gmap diag(free/m) Gmap^T；不能用人为小刚度代替零刚度。
- 任意约束与全约束：所有轴都保留在结果中；约束冲量/自由方程独立回算。
- 实际3×3×3 PySPH小块：dx=.5mm、rho=2700kg/m3、G=26.3GPa、
  Kbulk=2G(1+.33)/(3(1-2*.33))Pa。该制造弹性组件不是AA6061塑性AFSD。
  任意仿射位移使用独立3×3矩阵A，解析应变/应力/能量；小扰动当前三维自由表面邻域也测试。
- 位移场L参考、均匀应力参考相对2e-11；张量绝对容限按浮点尺度设置：梯度2e-12，Pa为1e-5。
- 代数手算/反力：atol2e-13、rtol2e-13；对称性/正性判据为矩阵范数尺度1e-12。
- 真SPH刚平移/刚转：以非零刚体字段测试K u的范数，界`2e-11*||K||2*||u||2`；
  不能要求浮点稠密矩阵乘积字面零，也不能直接清零结果。
- 独立伴随力/能量：归一残差<=2e-11，绝对下限1e-12N或1e-15J；
  有限差分方向导数用冻结的实际corrected_operators + internal_force独立路径，
  不调用新组件生成参考；中心差分epsilon=1e-6，归一误差<=2e-8。
- K零模统计：按`1e-10*max(abs(eigenvalue))`识别数值近零，本质是诊断，不截断K。
  刚体6模以外的零能模必须报告；A因真实M正定可解不等于K无零能模。
- 输入非有限/非法质量/不对称或明确负刚度/无效自由轴拒绝；退化三维邻域继承cond>1e6失败。

每个测试进程外60s，最多两轮有证据修复，保留所有日志。新物理AFSD运行NOT_RUN。
目前以上测试全部NOT_RUN；完成后只在下文追加实际结果，不追改容限和reference。

## 实际实现与首跑结果（随后追加，2026-10-03）

已新增`src/afsd_backend/pysph/structural_mobility.py`和独立
`tests/review/test_structural_mobility.py`，未修改旧物理/配置/driver。
上文“目前NOT_RUN”是事前状态，以下是真实首跑结果，不删除原合同。

- **14项PASS，0修复轮**；内0.132s、外0.4530272s、exit0，外部timeout60s。
- 日志：`logs/process_completion_20261003/mobility_tests_initial.log`。
- 2DOF非对角接触响应、质量-only、任意/全约束反力、输入保护均按独立手算通过。
- 实际27点PySPH邻域：三维仿射梯度/应力/能量、刚平移/刚转、旧伴随力独立路径、
  中心差分方向导数、扰动自由面/不均匀rho、标签独立与ID重排、退化邻域失败均通过。
- 81×81 K在事前1e-10相对尺度下有**6个数值近零特征值**，与3平移/3转动刚模相容；
  该样本未检出额外零模，**不证明所有SPH自由表面或大变形构形没有零能模**。
- 原始最小K特征值`-1.3708585838962082e-8 N/m`、
  非对称范数`2.043491705502407e-8 N/m`是矩阵乘法相对尺度内舍入，原数未裁剪。
  K保持原矩阵；未加epsilon刚度、未对称化回写、未伪逆或删刚模。
- 真实质量构造的dt=1e-7s A自由块cond=`2.2551298894894423`，全部81轴仍自由。
  A的正定来自质量项，不能当作K强制性或生产大步稳定认证。

实际命令（设置OPENBLAS/OMP/MKL_NUM_THREADS=1；外层用同一Python subprocess timeout=60）：

```powershell
& .\.envs\pysph\python.exe -B tests/review/test_structural_mobility.py
```

新增接口的`reaction_impulse_Ns`是**线性化结构方程**的规定轴反力冲量，
不是已经求得的真实AFSD执行器反力。自由轴尾差保留为`free_equation_residual_Ns`。
W是当前弹性预条件器响应，法/切接触律、预应力、FKS切线、旋转feature、热耦合、
事件路径及密度/几何导数均未在此求解；这些联合检查仍NOT_RUN。
中点W/2只对当前增量未知量约定成立，不能将其无条件粘贴到旧kick替换P/m。
刚转测试验证的是弹性算子刚模，不冒充带初应力的客观材料更新测试。

未运行真实AFSD或任何物理时间推进；不生成accepted状态、VTK或checkpoint。
可将本组件交给下一过程步的预条件/切线组装工作，不能据14PASS结束整个AFSD项目。
