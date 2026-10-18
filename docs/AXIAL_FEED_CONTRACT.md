# 局部轴向送料接口：事前十项合同

2026-10-02。ENGINEERING_APPROXIMATE_LOCAL_AXIAL_FEED_NOT_AFSD；不修改旧源/配置/FAIL。
继续PySPH同一pin。此阶段无旋转、热、生产、GPU或并行优化，不授完整复现标签。

## 十项有效输出

1 真实调用/坐标/参数合同；2 原始材料层轴向执行器与释放接口；
3 独立解析反力/功基准；4 有底座法向力的执行器功基准；
5 释放后自由响应及标签/ID不干预；6 真3D局部小棒—刚底座接入；
7 初始K/P、执行器/底座冲量与功、法向势能/离散残差分账；
8 非零态/释放点重启、REJECTED回滚及输出验证；
9 固定物理两档有预算试算；10 独立逐步读回、传力/接触事件、共同终时对照。
代码/实际测试/精度screen/未运行严格分开；十项输出不等于十项科学认证。

## 固定物理问题与调用路径

沿用原64粒子4x4x4、dx=.0005m、h=.00065m、rho2700kg/m3、G26.3GPa、nu.33，
Jaumann/J2/FKS参数不改，corrected梯度与伴随力成对。等温293.15K仅输入，不解温度。
中心x=y=0，世界z=-.01+(k+.5)*dx；实体高度2mm。
沿用geometry的刚底座z=-substrate_depth=-.01m，不把底座误当z=0。
工具静止、肩面z=-.009m，方孔宽8mm保持开放，tool_top=.02m；2mm小棒远离孔壁。
无基板刚化：本探针根本没有基板；全部64点为同材可变形局部方棒，不能冒充完整AFSD。
实际contact_rows只允许BASE(feature0)；底部点不固定，不用role1夹具。
法向n=+z, gap=z+.01m, radius=dx/2, Fn=657500*positive_part(radius-gap)N，
U=.5kn*positive_part(radius-gap)^2；初始恰触，U0≈0，所有位移真实受力求出。
此处底座为平面，其势能与实际feature一致，不复用旧肩面top-cohort假定。
friction=0明确用于隔离反力/法向功，不把已有摩擦当不存在或已全面验证。

固定材料子集z0>base_z+height-depth，depth=.5mm（顶16点）。选择只看原始z0，
不看tracer/ASRS/当前位置，也不采用旧domain_kick的锁x规则。
0<=t<400ns：仅此16点z规定-0.00212m/s，所有x/y及其他48点xyz自由。
t>=400ns：全部xyz自由，保留实际速度，不设零、不再规定材料位置。
初始顶部vz就是-0.00212m/s，必须显式计入K0=1.213488e-11J、P0z=-1.1448e-8Ns；
其余点初始静止。初始速度制备的功不伪造为第0步执行器功。
400ns指令行程仅0.848nm；这是小幅弹性送料边界探针，不是可见压实/沉积。

真实顺序保持：半漂移→3D NNPS/corrected L_before→旧material.update→中点应力/伴随内力
→contact_rows/normal_state→逐点kick（内力+底座+约束一起）→密度/速度kick→半漂移→验收。
kick不把初速偷偷覆盖成命令，它检查初速一致后用受约束mobility。
转折必须落完整步边界；跨转折时间步明确拒绝。法向力仍半漂移几何冻结。
dt5/2.5ns，160/320步至共同800ns。无空间/质量/模量/时间缩放。

## 独立字面参考（不调用被测函数生成）

m2kg,v0=(1,-2,-.5)m/s,Fother=(4,-8,6)N,dt=.25s，free=(1,1,0)：
v1=(1.5,-3,-.5)，Jact=(0,0,-1.5)Ns，Wact=.75J，Wother=5.5J，deltaK=6.25J。
叠加底座Fn=2N,n=+z,mu0：v1不变，Jact_z=-2Ns，Wact=1J，Wnormal=-.25J。
释放后v1不重置，在无contact且相同F、dt下：v2=(2,-4,.25)，Jact=0，
Wact=0，Wother=deltaK=8.5625J。错误锁x/y、反力漏账、释放归零均需测试捕获。
已有摩擦反向和旧法向解析回归保留；本探针不以mu0替换生产摩擦。

## 分账、输出与独立核验

Jn=dtFn*n，Ja=m*deltaV-Jinternal-Jn-Jt，Wact=Ja·v_sub_mid；反力符号为作用于材料，
执行器受到的反力为其负值。输出实际force_on_material=Ja/dt，非另算静力预期。
全局P-P0=累计Ja+Jn+Jt，角动量亦算初值；底座静止Wwall=0，不等于Wnormal=0。
deltaK+Q-Wact-Wwall-Winternal-Wnormal_rel=0；Rn=deltaU+Wnormal_rel不得补热。
全局K-K0+Es+U-U0+Pplastic+Q-Wact-Wwall保留原符号；Es仍诊断而非自由能。
没有无穷供料：系统质量固定，boundary mass flux=0，驱动行程不叫沉积长度/质量流。
每步源/config hash、全64点x/v/rho/S/ep/标签、i/j/g/B、全部kick子区间、原始反力/功、
控制前后phase/selected/active/free_axes、底座gap/Fn/U与真实工具位姿落盘。
接受checkpoint/释放边界checkpoint与REJECTED分存，后者禁止restart。T/e未解不伪造VTK。

## 事前验收、成本与边界

解析字面v/J/W atol1e-12 SI；标签交换/清零物理逐位相同，排序/restart x/v/rho1e-12，S1e-6Pa。
制动/反力不得作用到自由轴：abs<=1e-17Ns；受控z速度1e-12m/s、指令行程2e-17m。
所有点质量相对1e-12，rho变化<=.01，speed<=1m/s，位移<=10um，penetration<=.1dx；
cond(M)<=1e6，不许pinv；非有限立即拒绝。stable_dt同时受原声学界和.1sqrt(m/kn)限制。
空间功误差1e-12W+1e-10*sum|F||v|，净力1e-12N+1e-10sum|F|，
净矩1e-14Nm+1e-10sum|r cross F|；kick能量账1e-18J+1e-9sum绝对功。
全局能量仅screen：abs>1e-15J+.01*(K0+U0+|Wact|+|Wwall|+|Q|)记ENERGY_WARNING，不修补热。
同终时共同帧dt筛查：x1e-8m、v5e-4m/s、S1e5Pa、rho.1kg/m3、功/能1e-12J；
超出只保留SENSITIVE_WARNING，不能当收敛，也不阻止独立接口开发。
全轨迹需有底部实际接触、非驱动层非零速度、xz/yz耦合与rank3邻域；不按期望曲线驱动。
接触首次采样的区间不充当独立声波传播精确时刻；粗粒子算子存在离散传播尾部。

所有短回归/reader外60s；3D单次内170s外180s(<600s)。先单测与两步成本，
唯一粗档完整+独立读回通过后唯一细档。若失败保存原记录，不无限修同一参数。
第三dt、空间加密、长时、旋转孔壁动力传递、热、塑性流动、完整AFSD/三工况：NOT_RUN。
以前材料/空间FAIL和NORMAL_CONTACT敏感性警告继续有效。
