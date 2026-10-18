# Q1运行前合同：小型准静态受载可行性，不重写AFSD

2026-10-01，用户已要求按工程近似继续，不追求逐项高精度。
本轮只增加scripts下独立driver/config/tests。src、原配置、原失败不改；PySPH唯一后端。
保留真实质量/模量/时间，无GPU/并行/生产/AFSD续跑。每命令外部timeout60s。

## 物理与独立连续体参考

沿用E2立方体[0,.01m]^3、rho2700kg/m3、E70GPa、nu0、alpha1e-5，
参考NNPS与corrected梯度/成对伴随内力、梯形V、m=rhoV、h=1.3dx。
制造弹性材料沿用E2的高屈服接口，实际调用原material.update；非FKS/AA6061标定。
实际生产核、本构、接触、积分器不修改。该driver明确省略惯性，属于新的准静态基准。
不是把现有显式步长直接放大，更不能将本基准0.1s称为完成AFSD0.1s。

新增双向无摩擦roller：X=0约束ux，Y=0约束uy，Z=0约束uz；各面其余方向自由。
剩余三个外表面自然零牵引；不是单边接触，不是全域规定运动。
这与旧E2仅X面roller不同，目的是消除准静态不唯一的刚体模态，需保持所有Q1网格一致。

独立载荷f_x=V*lambda(t)*E*alpha*k*sin(kX)，f_y=f_z=0，k=pi/(2L)。
lambda=3s²-2s³，s=t/.1s，物理加载0→.1s；自由材料轨迹不预设。
连续体静力平衡解：ux=lambda*alpha/k*sin(kX)，sigma_xx=lambda*E*alpha*cos(kX)，
其余分量0；支承总反力Rx=-lambda*70N，Ry=Rz=0。
终态ux(L)=2e-7/pi m；总弹性能U=1.75e-6J；平衡加载外功同U。
离散体力积分本身有梯形误差：Fx=lambda*70*q*cot(q)N，q=pi/[4(n-1)]，
不人为归一到70N。反力应与实际离散体力相反；另报告相对连续体偏差。

## 试探更新与约束求解

参考梯度D沿用probe_boundary_sbp.extract_original（只提取，不调用SBP投影）。
nu0小应变弹性刚度块K_ab=E/2[delta_ab sum_k Dk^T V Dk + Db^T V Da]。
它只是原material.update/原配对力残差的固定弹性迭代矩阵，非新的生产材料。
试探梯度=(D u_trial-D u_accepted)/physical_dt；每次material.update都从
同一accepted S/ep开始，迭代不能累积塑性/热/功。总球应力E/3 tr(Du)。
试探残差r=f_internal+lambda*f_body；仅自由DOF解Kff delta_u=r_free。
约束反力r_reaction=-(f_internal+f_external)仅作用于约束DOF。
最多8次修正；无pinv、正则化、人工阻尼、模量/质量缩放或按参考答案设置外力。
仅收敛候选一次commit。失效/过期候选拒绝，accepted状态逐数组不变。

## 最小检查与精度政策

- NNPS矩阵cond(M)<=1e6原保护保留；K对称相对1e-11，Kff最小/最大特征值>1e-12。
  这是识别奇异/未约束模式，不是要求连续体解高精度；失败不静默加约束。
- 相对自由残差<=1e-9为线性/准Newton停止准则，不是论文精度阈值。
  ID精确、有限量、无隐式塑性、位移/L<1e-4、无状态污染；只接受本小变形域。
- 独立配对力与矩阵力/虚功检查相对1e-10；反力加外力、力矩独立核验。
- 梯形外功sum .5(f_old+f_new) dot delta_u与
  应力功sum V*.5(sigma_old+sigma_new):Ddelta_u各自记录；
  两者相差不得超过1e-7相对（与求解残差比宽100倍），不能残差加热。
- 与解析位移/应力/能量的偏差只报告，不作为本轮继续方法开发的精度硬门槛。
  已有闭循环/E2/自由面FAIL不变，不为全绿删改；本轮不修材料/边界算子。
- 同物理几何n3/5/7，4加载增量；n5另8增量，均到0.1s。
  先n3单案例预算探针：预计27点/54自由度、4增量、最多32迭代<10s；
  再实测n5建立n7预算，预测>60s则NOT_RUN。禁止用不同终时称收敛。

## 输出与范围

逐接受状态保存持久ID、参考/当前坐标、实际应力/ep、外力/反力、
physical_time、载荷系数、非线性迭代与功账；位移差/dt只标kinematic_rate，
不是求解了动量方程的速度。温度为输入不写成求解场；密度仅reference_density。
最后沿用上游NPZ诊断snapshot，restart_allowed=false；未实现生产重启不得借名通过。
保存source/config hash、报告、真实命令/log。源可继续由HANDOFF恢复。

PASS仅表示此受载基准可在大物理加载间隔上求得离散平衡，
不证明热塑性接触大步稳定，也不授予生产。下一关才是受载材料/FKS或多面接触之一。
