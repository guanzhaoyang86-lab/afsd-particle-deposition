# 旋转接触：两项有界推进（2026-10-02）

本轮实现两项：①旋转刚壁的显式接触点/功接口；②复用旧积分/接触的随时间
旋转肩面动力组件。没有改旧物理源码、旧Q7/Onset、原失败、上游pin和生产默认。
这是ENGINEERING_APPROXIMATE方法组件，不是新一轮AFSD沉积，不增加“真实沉积时间”。

## 一、旋转刚壁接口及其独立验证

新src/afsd_backend/pysph/rigid_contact_kinematics.py复用原geometry.py。
contact_rows输出persistent particle_id:feature_id，面编号0基座、1肩面、
2/3方孔X-/X+、4/5方孔Y-/Y+。方孔仍开放，角部两个面保留，不平均法向。
刚体坐标面ID随工具转动不改；基座不跟工具动；标签不在接口输入中。

明确区分：力作用粒子中心x；有限面投影c=x-gap*n；两处壁速V+omega×r。
这不是完整闭合刚体最近点算法，顶部/外侧/真实棱角接触仍不支持。
现粒子只有平移DOF：继续用中心的材料/壁速计算力功，标CENTRE_COLLOCATED_APPROXIMATE。
给定L时可诊断v(c)=v_i+L(c-x)，但没有其转置力分配，**不能**当生产力学已接通，
也不虚构spin。偏心力矩(c-x)×F及功差omega·deltaTorque只是诊断，不补成热。
contact_power保留signed摩擦功率，法向相对功独立，不把法向储能叫摩擦热。

17项独立回归全部通过，涵盖一般三维omega、90°旋转、开放孔/双面角部、
有限差分壁速、ID排列、标签接口隔离、非有限输入/埋入工具拒绝。
明确反例：x=(.45,.2,1.5)m，c=(.5,.2,1.5)m，omega_z=2/s、Vx=.3m/s，
wall(x)=(-.1,.9,0)，wall(c)=(-.1,1,0)m/s。只换壁速不换材料作用点会造出.1m/s伪滑。
F=(-10,3,4)N时偏心力矩=(0,-.2,.15)Nm、功差=.3W；同点与重构后的共转都零滑。
刻意非物理力样本raw摩擦功率=-9.5W也原样返回，不abs/裁零。

一个实际旧kick解析样本通过：m=Fn=mu=1，dt=.2s，wall=(.3,1,0)，
v0=(1.3,1,0)，正常外力平衡；v1=(1.1,1,0)，Jt=(-.2,0,0)Ns，
Q=.18J，Wwall=-.06J，deltaK=-.24J。它是冻结几何样本，不冒充旋转轨迹。

## 二、旋转肩面驱动自由质点的真实积分

scripts/probe_rotating_shoulder_component.py读取新SI制造配置。
复用PySPH ParticleArray、HostLeapfrog、原contact_kick.kick/contact_cone和checkpoint/VTK。
真实顺序由固定上游源码核对：

`半漂移 → 中点位置/时刻 → contact_rows → 原kick(约束mobility) → 速度kick → 半漂移`

没有事后修轨迹或再次覆盖自由速度。omega=20rad/s只属于刚壁；没有给自由点
预设圆周运动或偷偷加向心力。原法向Jn和法向执行器reaction都记录，xy完全自由。
4个制造点质量各.001kg、mu.2、规定Fn1N；z=.0013m固定，肩面z=.0015m，
初始xy=(.01,0)/(0,.012)/(-.01,0)/(0,-.012)m；初速0，终时.0008s，工具转角.016rad。
参数均为独立contact benchmark，不是AA6061质量/材料缩放，也不是论文300rpm工况。
这里的dx仅接触检测尺度，不是SPH空间分辨率。法向载荷规定，不在本轮标定penalty。

**这是平面xy动力学通过三维接触接口演化，不是三维可变形SPH固体。**
此组件无需NNPS或材料内力；并未替换AFSD基板/方棒为刚体。
未求stress/plasticity/thermal/density，不导出伪字段；PySPH checkpoint必需rho/e
明确为scratch，h仅元数据。当前肩面偏移c-x平行omega，故两处壁速相等；
不据这一特例宣称旋转孔壁真实面点功共轭已解决。

## 参考、实际命令和数值结果

独立RK4只写运动方程r'=v，m v'=mu Fn*(omega×r-v)/|omega×r-v|，
同步积分rawQ及壁功；不调用contact/kinematics/integrator生成参考。
reference采用literal物理参数，并校验其与配置一致、保存来源哈希。
1e-7和5e-8s参考步长各到共同.0008s，最小slip=.0400058m/s，始终高于
事前.03m/s域限；它只适用于持续滑动，不能当止滑/反向通用真解。
参考两档位置差1.735e-17m、速度差1.638e-15m/s、Q/W逐点差<1.6e-19J。

| 实际路径（各一次） | 步数/时长 | 外墙钟 | 末态位置最大参考差 | 速度最大参考差 |
| --- | --- | --- | --- | --- |
| dt=2e-5s | 40/.0008s | .616505s | 4.27371e-11m | 6.39239e-7m/s |
| dt=1e-5s | 80/.0008s | .908472s | 1.06719e-11m | 1.59934e-7m/s |

两档都COMPLETE_MANUFACTURED_COMPONENT_NOT_AFSD、exit0；没有重试/调参。
参考外墙钟.769624/1.229920s；最初成本微探针2步.0097995s，实际均低于60s硬限。
组件内部含输出耗时.139036/.259226s，数组800bytes，只是属性大小，不是峰值RSS。
不拿该微型接触组件速度预测论文大规模AFSD成本或CFL。

粗/细末态Q=8.960054022142603e-5 / 8.960054022208344e-5J；
Wwall=.00014079985542973322 / .0001407998540721458J；
K=5.119931520830716e-5 / 5.119931385006234e-5J。
原记录全历史最大能量残差3.1564e-20/5.0640e-20J，独立逐子步读回结论见下节；
不是通过往热源填差获得。与ref的改善是本组件两档证据，不授予连续体收敛阶。

## 测试与证据路径

logs/rotating_contact_20261002/保存真实命令输出，COMMANDS.md可重跑。
kinematics_final.log：17PASS；component_final_gate.log：11PASS；old_contact_events.log：12PASS。
新11项含两步自由加速literal、法向反力/功、标签交换、ID重排、真实checkpoint继续一致、
失败回滚REJECTED、VTK只含解出的字段、独立RK4无接触调用、配置漂移拒绝。
首次8项里有1项测试要求浮点reaction_work严格等0，实际5.08e-24J而失败；
仅改该断言为0±1e-18J（比原合同更严），没有改物理/试跑参数，初日志完整保留。
同组重复测试不相加。历史大套件本轮NOT_RUN，原FAIL未删改。

outputs/rotating_contact_20261002/有两参考JSON、两目录history.json、
contact_intervals.jsonl、41/81帧真实VTK、last_accepted.npz和manifest来源/内容哈希。
该checkpoint仅支持制造组件重启，production_authorized=false，不是AFSD生产重启许可。
独立读回已实际执行，PASS_COMPONENT_READBACK_NOT_AFSD，exit0/外.365338s。
读回脚本不导入求解器/几何/参考程序，17项独立合成与篡改回归PASS；
最终新测试合计45PASS（17接口+11短组件+17读回），旧选定接触12PASS，重复不相加。
逐子步另核对Coulomb方向/容量和原舍入界，不能只看Q>=0；该代数检查没有改物理或
放宽事前误差界。所有实际raw Q为正，最小粗1.68022e-7/细8.20113e-8J。
重算全程最大能量残差粗1.02972e-19/细6.77864e-20J；与程序记录的微小差别来自
求和顺序，两者均保存。角动量z末值7.03999277149e-6/7.03999270361e-6 kg m2/s，
等于中点作用位置的累计转矩冲量，Wwall=20*deltaLz逐帧通过。
全部122个VTK帧与history实际位置/速度/ID逐一一致；未伪造场。
最终状态与全部参考差：readback.json；真实输出日志readback.log。

## 下一项，不退回旧小块反复试跑

把已验证旋转肩面运动/功接口接到**有真实内力的可变形SPH小块**，从短程等温
制造载荷开始，明确正常约束/载荷来源，保留原材料能量和空间近似风险。
先读solver.py的evaluate、contact_kick.domain_kick和mechanics/material现有路径；
不要复制点质量演化去代替金属塑性，不给自由材料指定轨迹。
旋转孔壁的面点重构/转置力传递、正常开离、轴向供料、热/FKS/有限应变另列缺项。
当前真实AFSD历史仍2us/191步、无持续沉积；三工况生产不授权。
