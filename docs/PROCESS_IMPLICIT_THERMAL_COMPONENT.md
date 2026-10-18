# 保守隐式导热组件：事前合同（2026-10-03）

分类ENGINEERING_APPROXIMATE_IMPLICIT_THERMAL_COMPONENT_NOT_INTEGRATED。
当前updated_trial实际使用coupled_thermal显式T0导热；本组件不改它、不跑真实AFSD。
几何/质量/密度/导热系数/接触热导在当前候选区间内冻结，是明确时间近似，不是全耦合宏步认证。

## 方程与单位

沿用thermal.py原对称径向系数gij[W/K]>=0：
gij=mi*mj/(rhoi*rhoj)*4ki*kj/(ki+kj)*[-rij·grad_iW]/|rij|²。
机械修正后的非径向梯度不允许混入。每无向边只出现一次，异密度仍质量配对守恒。
Ci=mi*cpi[J/K]，Hi=sum hwall[W/K]，Bi=sum hwall*Twall[W]。
壁温Bi/Hi必须来自真实接触热库；Hi=0时Bi必须为0。

令L为对称图Laplacian，求温增而不是大绝对温度：

    [diag(C)+dt*(L+diag(H))] dT
      = Qplastic*chi + Qfriction*partition + dt*(B-H*T0-L*T0)
    T1=T0+dT

Q均为实际区间能量[J/粒子]；机械return或残差不能作为热源。
Qi_pair=dt*sum_j gij*(T1_j-T1_i)，Qi_contact=dt*(Bi-Hi*T1_i)。
每粒子Ci*(T1-T0)与上述来源闭合；内部热对和为0。
Hi非负且Ci正使矩阵SPD；BE无显式导热CFL限制，但一阶时间误差仍需同终时检查。
参考：[PETSc TSBEULER官方说明](https://petsc.org/release/manualpages/TS/TSBEULER/)
及[FEniCSx热方程教程](https://jsdokken.com/dolfinx-tutorial/chapter2/heat_equation.html)（2026-10-03访问）。
仅使用方法依据，不安装/复制这些软件，不更换PySPH后端。

## 事前独立测试与容限

- 两体C1=C2=1J/K、g=1W/K、T0=[300,400]K、dt=1s：BE解[1000/3,1100/3]K；
  连续解[350-50exp(-2t),350+50exp(-2t)]；同终时N=10/20/40，误差下降比应在1.7–2.2（渐近一阶）。
- 单体C=2J/K、H=3W/K、Tw=500K、T0=300K、Q=4J、dt=.5s：T1=1354/3.5K。
- 单体无导热：塑性10J、chi=.9，摩擦4J、partition=.25，C=2J/K，ΔT=5K；
  未热化塑性1J、摩擦壁热3J分列，不能漏算/二次分配。
- 对称图总能量、零源最大值原理、异容量、重排、重复边/负热导/错梯度/越固相域明确失败。
- 独立解析温度绝对2e-12K，矩阵/账尺度差1e-12J+1e-10*实际绝对热量；
  重排温度2e-12K；原径向系数手算相对1e-13。这些是代数回归而非论文误差容限。
- raw摩擦能量带符号保留，只接受调用方提供的非负舍入证书范围内负值；不abs/max裁热。
- 只读小规模PySPH邻域接口可核对旧显式导热同系数，不得将该核对作为独立参考的替代。

实现：src/afsd_backend/pysph/implicit_thermal.py；测试tests/review/test_implicit_thermal.py。
测试实际结果另列，不从本事前合同推断PASS。旧材料自由能并未由热求解器修复。

## 首次短测试与一次离线接线检查登记

thermal_tests_initial.log：11项全部PASS，外0.221689s，无修复。
同终时1s的连续解析温度误差为1.308515/.665417/.335520K（10/20/40步），支持一阶而非大步精确。

接下来只对旧4个接受stage分别重算热分量：保持每stage原始dt、原始T0、真实中点NNPS/密度及原热源，
壁热导和h*Twall从逐接触face独立重建。舍入证书从实际接触子区间累加，不随结果放宽。
用各旧stage初态分别读回，**不把新温度继续用于下一旧stage，不推进全耦合物理时间**。
核对原导热行和界、独立手算pair系数、逐点mcΔT/壁/源账及输入hash；容限沿上节。
与旧显式温度差只报告，不把相同答案当要求，也不将热离线PASS当新AFSD接受态。
入口scripts/replay_process_thermal.py，外60s；最多一次实际离线重放，失败留证后最多2轮修复。
