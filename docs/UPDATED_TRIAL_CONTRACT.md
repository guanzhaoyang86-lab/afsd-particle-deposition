# 整体候选与有界全局试步：事前合同（2026-10-03）

目标：把更新构形/材料/真实接触/热学接成同一可求解候选，争取一个真实全局试步。
不以单测数量或文档数授予75%。75%工程里程碑至少要求可负担的过程推进及实际材料输运证据，
而非仅微秒启动、单个求根或动画；未达到就明确缺口，不为了数字改定义。
既有342点/700ns与所有旧FAIL保留；唯一后端仍固定PySPH。无GPU/并行/依赖安装/缩放。

## 命名数值变体与真实方程

新增UPDATED_CONFIGURATION_EVENT_AVERAGE_STAGE_IMPLICIT_APPROXIMATE。
不是静态平衡，不删除质量/惯性；不是原HostLeapfrog的同公式续算，也不冒称完整有限应变自由能模型。
只复用现有粒子/NNPS/本构/接触事件/热与输出接口，旧入口不改。
accepted输入只读，trial与求根迭代不得累计ep/热、改标签或写accepted checkpoint。

candidate包含x1/v1/rho1/T1及独立mean_velocity_m_s，即候选区间时间平均速度vb。
使用xmid=(x0+x1)/2、rhomid=(rho0+rho1)/2，Tmid=(T0+T1)/2；
midpoint PA的速度为vb，重新建立当前3D NNPS/M/B/L。
材料从同一accepted S0/ep0用L、dt、Tmid调用原Jaumann/J2/FKS；压力闭合用rhomid；
同一B/g配中点应力伴随内力。不得继承Q7冻结参考D/B/常rho/制造Y。

原event-aware边界kick以真实v0、质量、候选Fint和冻结中点几何/工具运动求解；
每面共享约束后mobility，只有原始顶层z规定运动。返回v1_map与真正各子区间积分的vb_map。
映射x1_map=x0+dt*vb_map，不以(v0+v1)/2代替反向滑移段的平均速度。
连续方程用Rrho=(rho1-rho0)/dt+rhomid*tr(L)，对应rho1_map=rho0*(1-.5dt trL)/(1+.5dt trL)。
非正密度/分母、NaN/温域/穿透/condition>1e6均失败，不裁剪/伪逆。

热来源仍旧coupled_thermal：塑性P、实际raw摩擦Q、接触导热、原始径向核成对导热。
明确split：FKS和接触导热反馈Tmid，内部导热及热增量起点仍T0，受真实热行和dt界。
没有把return、法向势能时间差或机械残差补成热。此为stage-implicit源反馈，不叫全隐式导热。
检查有限feature before/mid/after非零势能切换，未支持边缘仍拒绝，不能丢失旧面储能。

## 独立联合基准

- 更新三维仿射梯度：独立矩阵A，当前邻域L=A，容限1e-10/s；不同候选几何须重建B。
- 连续方程：独立各向同性L=aI时，rho1/rho0=(1-1.5a dt)/(1+1.5a dt)，不谎称任意有限步等于指数体积解。
- 非零初应力刚转：独立sin/cos旋转参考；材料给定常spin的指数旋转与坐标中点Cayley有限步误差分开。
- 空间功/净内力矩沿旧1e-12+1e-10尺度；纯候选重复/失败/重排/换标签不能改accepted。
- 原热制造参考仅测试law可注入，实际默认FKS；真实输出不得用制造Y/轨迹驱动。
- 反向接触平均速度：原m1/u0=1/Fother=-2/muFn1/dt1基准真实位移=-1/18，Q=7/18；
  不用端点中点1/6位移或-1/18假热替代。保留多面/约束与实际子区间功。

## 全局有界迭代与运行决定

环境无SciPy，采用可审查的有界阻尼Picard解当前候选映射，不安装/实现通用大型Newton框架。
迭代松弛.75只用于非接受试探，不是物理阻尼；最大32次evaluate，不使用pinv、调参追绿或自动时间缩放。
验收原始map残差：x<=1e-12m，v1/vb<=1e-9m/s，rho相对<=1e-10，T<=5e-12K；
反力/空间/热/finite检查也必须通过；这些是新非线性解收敛精度，不是论文误差目标。
求根成功仍须检查温增与真实热源账、事件平均位移、质量/密度/穿透/feature；失败candidate明确REJECTED。

唯一真实数据来源：旧short_07us的plastic_accepted.npz（step62，保留非零应力/ep/T/v，不清零）。
先核原manifest/source/artifact/hash和参数，显式作为新变体的初态导入，不伪造旧hash一致的常规restart。
候选dt固定5e-8s(50ns)，几何/物性/速度/接触刚度/温度完全不改；不声称这是可生产宏步。
先一次纯trial成本（外90s），记录全部接触子区间、dt、时间/内存、source/config/初态hash。
保守成本=1.5*完整trial墙钟*32+30，<=450s才启动一次global attempt，内540/外570。
超预算或局部接触/材料阻塞就记录BLOCKED/NOT_RUN，不在同一点无限修复，不追加别的dt追成功。
全球迭代成功后保存真实新变体接受态与检查账；若不满足任何接受条件仅保存REJECTED且旧接受态原样保留。
单次真实尝试严格<=600s。测试/独立读回各外60s。后续多步、同终时敏感性、三工况、持续入口、沉积NOT_RUN。

## 交付与停点

UpdatedTrial、可重跑全局试步入口、独立联合tests/实际日志、候选/接受/拒绝来源与结果schema。
最多两轮本轮核心实现修复，所有失败保留；不根据结果调整材料或容限。
若未获得过程尺度输运，明确尚未达到有效75%，同时交付本轮实际接通的完整链和下一主要瓶颈。
