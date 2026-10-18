# 局部轴向驱动—可变形材料—刚底座—释放：十项交付

2026-10-02。分类：`ENGINEERING_APPROXIMATE_LOCAL_AXIAL_FEED_NOT_AFSD`。
后端仍为 PySPH `1cd3dea9cfaac7d9f1bda2a2491ce4b42a31dbe5`。
本轮是边界组件集成，不是新求解器、论文复现或沉积结果。没有修改旧物理源码、旧参数、旧FAIL。
事前参数/容限见 AXIAL_FEED_CONTRACT.md；实际命令见 ../logs/axial_feed_20261002/COMMANDS.md。

## 十项输出及实际状态

|项|有效输出|实现/证据|实际状态与边界|
|---|---|---|---|
|1|真实调用、坐标、控制范围、量纲及预算冻结|AXIAL_FEED_CONTRACT.md；configs/benchmarks/axial_feed.json|已落盘且在运行source hash中；制造问题，不新增论文事实|
|2|原始材料顶层z-only执行器，400ns自由释放|src/afsd_backend/pysph/axial_actuator.py 的 control/validate_interval|PASS：12接口回归之一；不锁x/y，不按标签选施力对象|
|3|反力/外功独立解析基准|scripts/probe_axial_actuator_reference.py；outputs/axial_feed_20261002/actuator_reference.json|PASS：实际旧kick/HostLeapfrog两步，不用被测函数制造参考|
|4|叠加底座法向载荷的执行器功基准|同一reference中的 Fn=2N 分支|PASS：正常材料功为负，刚壁功为0；反力与各项功没有混用|
|5|释放、标签交换/清零、排序不改变物理|tests/review/test_axial_actuator.py 与 test_axial_feed_block.py|PASS：观察标签与缓存actuator/boundary_role不定义受力；自由点不被清零|
|6|全部64点可变形的真3D局部棒—刚底座|scripts/probe_axial_feed.py 的 Block|PASS组件运行：取消本组件底部材料夹持；实际3D邻域/内力/应力，无二维复制|
|7|初始K/P/L/U与执行器/底座反力、功、弹簧储能分账|Block._callback/_validate_candidate/record；逐步journal|PASS功/冲量身份式；Es仍诊断，能量screen有WARNING，不授自由能认证|
|8|非零态与实际释放点checkpoint/原子回滚|scripts/check_axial_feed_checkpoint.py；两份 *_restart.json|PASS：400ns真实checkpoint续两步，x/v/rho/S/ep逐位一致；REJECTED不能restart|
|9|先测成本，再同物理两档完整短算例|AXIAL_FEED_COST_PREFLIGHT.json；dt5e9、dt2p5e9目录|COMPLETE：160/320步均至800ns；单次外8.254/14.456s|
|10|独立逐子步/VTK读回、接触事件及共同终时对照|scripts/inspect_axial_feed.py；comparison.json|读回PASS，482帧；161共同帧对照SENSITIVE_WARNING，不称收敛|

十项是十种可检查的产物，**不是十项科学认证**。短回归最终独立计数40PASS，旧选定31PASS。
本轮没有第三次调参运行来消除警告。原失败、早期日志和旧阶段敏感性结果均保留。

## 物理范围和真实调用顺序

64个粒子，4×4×4，2mm立方局部方棒，dx=0.5mm、h=0.65mm、总质量2.16e-5kg。
世界底座z=-0.01m，粒子中心z=-0.01+(k+0.5)dx；全部材料点可变形，没有固定底16点。
本探针没有基板，不能把它描述为“完整棒—基板沉积系统”。工具方孔宽8mm、保持开放、静止；
2mm棒远离孔壁，只发生BASE接触，孔壁的动力传递本轮NOT_RUN。

原始z0位于顶0.5mm层的16点，在0≤t<400ns仅规定vz=-0.00212m/s。
所有x/y及其余48点xyz自由；400ns以后全部xyz自由，既不清零速度也不回写位置。
初速顶层已经是规定速度，K0=1.213488e-11J，P0z=-1.1448e-8Ns；这属于初态而不是未记账的第0步外功。
指令总行程只有0.848nm。系统封闭、质量流率0，没有连续供料或可见塑性沉积。

实际顺序（源：probe_axial_feed.Block._callback、继承的 Block.advance、HostLeapfrog）：

```text
接受态
  → 半漂移（已有PySPH积分设施）
  → 3D NNPS → corrected梯度 L_before
  → 原 Jaumann/J2/FKS material.update
  → 中点应力 + adjoint内力
  → 实际geometry BASE间隙/法向力
  → 逐点contact kick：内力、法向、受约束mobility一起，记录反力
  → 速度/密度kick → 半漂移
  → candidate验收 → 接受态/checkpoint/VTK
                 ↘ 拒绝则保留REJECTED并回滚
```

不是先contact后直接覆盖材料速度。释放只改变control的活动轴；转折跨步会拒绝。
Jaumann/J2/FKS、核、校正空间算子/伴随内力、旧contact_kick均未修改。
293.15K为规定材料输入，thermal NOT_SOLVED；mu=0显式隔离本轮法向/执行器功，未替换生产摩擦。

## 独立参考、回归与本轮修复

独立字面参考：m=2kg，v0=(1,-2,-0.5)m/s，Fother=(4,-8,6)N，dt=0.25s。

|阶段|v_after (m/s)|Jact (Ns)|Wact (J)|Wother (J)|Wnormal (J)|deltaK (J)|
|---|---|---|---|---|---|---|
|顶层z驱动、无接触|(1.5,-3,-0.5)|(0,0,-1.5)|0.75|5.5|0|6.25|
|相同驱动、Fn=2N向+z|(1.5,-3,-0.5)|(0,0,-2)|1|5.5|-0.25|6.25|
|随后释放、无接触|(2,-4,0.25)|(0,0,0)|0|8.5625|0|8.5625|

坐标初值(0,0,2)m，两步为(0.3125,-0.625,1.875)、(0.75,-1.5,1.84375)m。
以上来自独立常力/受约束运动积分，只有算完后才与被测输出比较。

|测试证据（logs/axial_feed_20261002/）|实际结果|解释|
|---|---|---|
|actuator_unit_initial.log|12PASS|选择、单轴约束、释放、字面反力/功、参数域、初始账|
|block_unit_initial.log|15PASS + 1ERROR|保留首次排序后checkpoint初始角动量不完全一致错误|
|block_unit_final.log|16PASS|修复后首次通过，不能与下一行重复计数|
|block_unit_source_final.log|16PASS|最终运行源码的短组件测试|
|readback_unit_initial.log / readback_unit_complete.log|各12PASS，计12|独立字面算术/故意篡改反力、功、控制和间隙的反例|
|old_contact_tests.log|12PASS|旧接触事件回归，本轮未修改该算法|
|old_normal_tests.log|19PASS|旧法向/移动壁/开离回归|
|actuator_reference.log|PASS|两个独立单点驱动—释放路径实际运行|
|coarse_readback_initial.log / comparison.log|PASS读回|粗档先单独验收才启动细档；共同比较仍有敏感性警告|
|dt5e9_restart.log / dt2p5e9_restart.log|均PASS|两档实际400ns释放checkpoint各续两步|

最小修复只在**新增**probe_axial_feed.Block.__init__：初始K/P/L/U按persistent ID固定求和顺序。
首次错误来自排列改变浮点求和顺序、checkpoint严格等值检查失败；没有把L0设为0或放宽容限。
另外将candidate末态时间统一为(step+1)*dt，避免释放边界元数据与接受态时间因浮点加法不同。
最终源码在成本探针与两档运行前冻结，运行期间没有改物理或参数。
test_15短释放回归将非零短态放在释放边界，仅为事件单元测试；真正400ns路径由两次完整运行和实际重启另外证明。

## 成本、真实结果及敏感性

事前两步探针外0.484s，计算段0.0188/0.0415s；全流程线性预算约8.95/17.90s，保守上限60/90s。
实际单次内部170s、外部180s硬预算；短回归、读回、重启外部60s。
初始声学dt上界1.04907e-8s，法向接触上界7.16455e-8s；两档满足两者且不缩放模量/质量/时间。
ParticleArray占26880bytes，仅粒子数组，不是峰值进程RAM。峰值RSS未测；不据此推算论文全尺度。

|实测量|dt=5ns|dt=2.5ns|
|---|---|---|
|接受步数 / 共同终时|160 / 800ns|320 / 800ns|
|内部 / 外部墙钟|6.897 / 8.254s|12.624 / 14.456s|
|含输出内部步速|23.20步/s|25.35步/s|
|物理时间/外部墙钟|约9.69e-8s/s|约5.53e-8s/s|
|读回帧数 / 最小邻域rank|161 / 3|321 / 3|
|非驱动层最大速度|0.00431935m/s|0.00428258m/s|
|最大xz/yz应力|约8360Pa|约8021Pa|
|末态底座活动接触数 / ep_max|16 / 0|16 / 0|
|执行器累计对材料做功|9.421290366e-11J|9.419096092e-11J|
|材料法向接触功|-1.654016258e-11J|-1.644266553e-11J|
|末态弹簧能|1.654028381e-11J|1.644269752e-11J|
|法向时间离散残差累计|1.212330518e-16J|3.198688119e-17J|
|末态全局能量诊断残差|5.737143918e-12J|2.623561495e-12J|
|ENERGY_WARNING帧数|141|156|

刚底座静止，壁功为0，而材料法向功不为0。残差原样保存，不裁剪、不补成热。
两档400ns重启分别续至410ns/405ns，与各自已保存的接受轨迹x/v/rho/S/ep最大差全部0。

161个共同帧最大绝对分量差（不是向量/Frobenius范数）：位置8.51883e-12m、速度1.79850e-4m/s、全应力2886.95Pa、
密度8.35996e-5kg/m3，均在事前screen内。累计内力功1.72844e-12J、K 1.64395e-12J、
Es诊断2.07015e-12J、全局能量残差3.11358e-12J超过1e-12J的事前screen，
因此整体为 **SENSITIVE_WARNING**，不是时间步已收敛，更不是全轨迹连续体解验证。
能量残差随本次减半步长下降，不足以推断非线性材料闭循环问题消失。

底部16个ID首次活动接触在粗档[0,5ns]、细档[0,2.5ns]采样区间。
这只是该离散初值问题的活动性变化，不能冒充独立精确弹性波到达时刻。
读回从已保存i/j/g/B及独立算术检查功共轭、净内力/矩、冲量/功、法向势能和VTK一致；
没有导入被测contact/material函数生成参考，但也**没有独立解完整非线性材料轨迹**。

## 风险与下一轮边界

- PASS：新接口、字面反力/功、标签隔离/排序、原子回滚、两个短3D的账/读回、真实释放重启。
- WARNING：本轮能量screen/时间步敏感性；上一NORMAL_CONTACT阶段速度/应力警告继续保留。
- 旧FAIL：Jaumann闭循环和空间自由表面/材料精度反例不因本轮通过而改判；见旧科学报告。
- NOT_RUN：本轮第三dt、同几何空间加密、旋转方孔动力传递、连续质量入口、热、塑性流动、
  完整可变形棒+基板耦合与三工况生产。本轮ep实际为0，不能声称有塑性送料沉积。
- 定量论文验证BLOCKED：材料缺项、模型适用性/误差和实验测量算子/原始数据仍有限制。

下一工作单元不是再调同一弹性块的能量比值，而是将已分别验证的轴向边界与开放旋转方孔接触衔接：

1. 先检查 geometry.py/surfaces、rigid_contact_kinematics.py/contact_rows、contact_kick.py/domain_kick。
   有限face切换的gap/势能要一致；孔壁投影点速度若重构，力必须取同一映射的伴随，不能只改壁速。
2. 给旋转孔壁—执行器单个/角部多面接触建立独立功/反力参考，覆盖开孔无力、共转零滑、
   自由xy与轴向约束耦合、feature切换。失败只修最小接口，不更换本构或恢复raw算子。
3. 上述接口通过后再加一个新命名短局部组件，把控制接口接入开放孔；先测两步预算，
   保持真实材料自由度，源/参数/状态快照完整。不直接启动论文长程或三工况。

实现边界/源码清单：新增 axial_actuator.py；五个 scripts/*axial*.py（probe/reference/reader/restart/freeze）；
三个 tests/review/test_axial*.py；一份新config及本组合同/预算/报告/日志/输出。
旧文件仅HANDOFF、IMPLEMENTATION_PLAN、OUTPUT_SCHEMA更新交接；冻结清单核对其余历史文件完整性。
AXIAL_FEED_FREEZE.json只证明证据完整性与后端pin，不证明科学精度，不授权生产。
