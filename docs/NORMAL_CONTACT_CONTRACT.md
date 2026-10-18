# 法向承载、卸载与开离：十项事前合同

2026-10-02。ENGINEERING_APPROXIMATE，不宣称复现论文。旧PySPH pin、源码、配置、
FAIL与ROTATING_DEFORMABLE_FREEZE全部保留。新companion复用原设施，不做GPU/并行优化。

## 十项各自的可检验输出

1 调用映射/合同；2 normal_penalty单边gap/Fn/U接口及上游一致性；
3 静壁解析碰撞/回弹/开离；4 匀速壁解析冲量/外功；5 抬升工具的几何/壁速一致接口；
6 解除顶层法向夹具的3D SPH适配；7 独立法向功/势能/离散误差账；
8 标签/ID/非零态重启/拒绝回滚；9 固定物理粗细短程；10 独立读回/开离事件/同终时比较。
每项要有代码或数据证据，不能用文档个数代替进展。关卡失败留FAIL/BLOCKED，后项NOT_RUN。

## 已查真实源码顺序与范围

integrator.HostLeapfrog=半漂移→回调→速度/密度kick→半漂移。
回调仍corrected_operators/3D NNPS→material.update(L_before)→中点应力/配对内力
→每点contact_kick.kick含约束→材料提交；等温，没有热更新。
contact_kick.domain_kick现有normal_forces=wall.law.kn*positive_overlap，但无底层xyz
夹具角色，不能直接沿用其入口约束。复用既有逐粒子kick以保留明示底层边界。
geometry.surfaces原壁速z恒0；仅改layer(t)会错失壁功。新接口同步移动layer/tool_top
并加Vz至工具所有接触的中心/面点壁速与translation_velocity；base不跟工具动。

## 法向定义及独立解析参考

signed_delta=radius-gap；delta=signed_delta若>0，否则0。这是单边几何接触律，
不是裁剪耗散或温度。Fn=kn*delta，U=.5kn*delta²，无黏附/吸引/黏性阻尼。
保留signed gap、active、Fn、U与其单位。实际normal law与固定上游
RigidBodyWallCollision(en=1)法向项一致；参考值由字面公式独立生成。

单点m=1kg,k=4N/m,radius=1m，初始z=1m，壁z=Vt、初速vz=V-1m/s。
相对固有频率2/s，接触时0<=t<=pi/2：delta=.5sin(2t)，u=-cos(2t)，
Fn=2sin(2t)，U=.5sin²(2t)J；t=pi/4时delta=.5m,U=.5J；
t=pi/2脱离、u=+1m/s、累计Jn=2Ns。此后Fn=0，自由开离，不拉回墙上。
终时2s：z=1+2-pi/2+2V，vz=1+V。
V=0时初末K=.5J、Wwall=0；V=.2时初K=.32J、末K=.72J、Wwall=.4J。
Qfriction恒0不是将负数裁零，因为本基准mu显式0且仅有法向运动。
dt=.005/.0025s，各400/800步到共同2s；参考必须独立，不能调用被测normal/kick。

实际DKD冻结法向力不严格守物理K+U。记录
Rn=ΔU+sum Jn·(v_sub_mid-vwall)，不把Rn补成热。
对持续接触/匀速平壁的孤立弹簧，有Rn=(kn*dt²/(4m))*ΔK_relative；
修正离散不变量(1-kn*dt²/(4m))*K_relative+U仅作为数值分析，不能改称物理能量。
开离边界仅步采样定位，报告包围区间/时间分辨率，不冒称精确连续事件根。

## 三维短程物理问题（与旧固定z夹具是新命名变体）

继承64粒子/4层/dx=.5mm实体2mm立方体、G26.3GPa/rho2700/nu.33、
原校正SPH/Jaumann/J2/FKS与293.15K规定输入。参数extends旧SI配置，新覆盖显式来源。
底16点xyz夹具；其余48点全部xyz自由，包括顶16点，top标签不再是法向约束。
保留原kn=.05Gdx²/.0005=657500N/m；面牵引刚度kn/dx²=2.63e12Pa/m。
初始layer=.002-1/kn m，每顶点delta=1/kn≈1.520913um、Fn0=1N，
U0=16/(2kn)J。不是旧50um重叠、不是静力预平衡。初始弹簧能明确供给；材料初应力零。
没有规定材料压实轨迹。肩面先保持到200ns，然后仅刚工具以+20m/s抬升，
到600ns总抬升8um；omega20rad/s，traverse0，孔宽8mm始终开放。
这一人为程序只测卸载/开离，不是论文启动/预热/送料程序。该块不在孔内，
不能据此认证旋转孔壁面点功传递；也未设置无限供料。

dt=5/2.5ns，各120/240步到共同600ns，物理几何/刚度/温度/边界不变。
稳定界同时要求旧声学dt及.1sqrt(min(m)/kn)，且抬升转折准确落于步边界。
正常穿透<=.1dx，mass relative1e-12，rho偏差<=.01，speed<=1m/s，
位移<=1e-5m，cond(M)<=1e6，非有限立即失败；无pinv/温度/屈服/耗散裁剪。
Es仍是旧本构诊断；全局K+Es+Pplastic+Q+U-U0-Wwall-Wconstraint保留原符号。
Jn·vwall与Jn·v_material分开，接触纯动能恒等式必须包含正常相对功，
不得错误沿用旧固定z的零正常功假设。

## 事前验收与预算

接口解析代数：力1e-12N、势能/功1e-12J，几何/壁速1e-12SI；无吸引，
Fn/U在开离状态必须严格为0；负signed_delta仍保留。共转/平移同点功检查继承。
单点完整轨迹：末位置abs5e-5m、速度2e-5m/s、最大物理能量残差1e-4J；
脱离采样时间abs<=2dt；累计Jn2Ns abs2e-5、移动壁功.4J abs1e-5。
DT细化比值报告，不反复改参数追求特定比值。
三维代数：空间功1e-12W+1e-10*sum|F||v|；净力1e-12N+1e-10*sum|F|；
净矩1e-14Nm+1e-10*sum|x cross F|；接触动能账1e-15J+1e-9*绝对功和。
自由轴reaction<=1e-17Ns；底夹具速度1e-12m/s/位置2e-17m；不拉住顶层。
全局Es/正常时间层偏差告警为abs>1e-12J+.01*(U0+|Wwall|+|Wconstraint|+Q)，
不据告警范围内授予有限应变能量认证，不用残差补热。
标签交换/清零物理逐位相同；非零态按ID重排/重启位置速度密度atol1e-12，S1e-6Pa。
两完整三维档同终时逐共同帧筛查：position1e-7m，velocity1e-2m/s，stress1e5Pa，
density1kg/m3，Q/W1e-7J；active转变逐ID记录区间，差<=2个粗dt为事件筛查。
对两个刚好同终时结果才比较，禁止不同失败终时叫收敛。

单点/回归外60s内50s；三维每次外180s内170s，均<=600s。
先短回归和成本探针，再唯一粗档，粗档完整+独立读回通过后唯一细档。
源代码/配置hash保存；独立reader不得导入求解器生成参考。
last accepted与REJECTED分存；预抬升checkpoint与重启位置/ID/材料/ledger均需检查。
旧FAIL、历史三维AFSD2us/191步记录不改变；本阶段未求温度不输出假热场。
塑性流动验证、局部送料、热、空间加密、第三dt、生产/长时NOT_RUN。
