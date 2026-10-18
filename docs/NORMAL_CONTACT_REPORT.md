# 法向承载—卸载—开离：十项有效输出

2026-10-02。当前层级：ENGINEERING_APPROXIMATE_COMPONENT；不是完整AFSD或论文复现。
沿用PySPH `1cd3dea9cfaac7d9f1bda2a2491ce4b42a31dbe5`。旧物理/配置/失败记录未改。
事前物理、容限、预算见NORMAL_CONTACT_CONTRACT.md；成本探针见NORMAL_CONTACT_COST_PREFLIGHT.json。

## 十项交付及实测状态

所有路径以项目根为基准；本轮数据根为`outputs/normal_contact_20261002/`，日志同后缀。

| 项 | 实际输出 | 状态与证据 |
|---|---|---|
| 1 | 真实调用顺序、初始预压能、法向符号与预算合同 | DONE；NORMAL_CONTACT_CONTRACT；不沿用旧顶层z夹具含义 |
| 2 | `normal_penalty.normal_state`有符号重叠/Fn/U，实际上游无阻尼弹簧一致性 | PASS；test_normal_penalty 01–05 |
| 3 | 静壁碰撞、压缩、回弹、开离的连续解析对照 | PASS；static_dt005_guarded/static_dt0025及static_pair.json，400/800步共同2s |
| 4 | 移动壁冲量与壁功基准 | PASS；moving_dt005/moving_dt0025及moving_pair.json；解析Jn=2Ns、Wwall=.4J |
| 5 | 工具抬升同步移动几何和壁速，底座不动、方孔开放 | PASS；moving_tool_rows/tool_state，单元06–10；非材料轨迹驱动 |
| 6 | `probe_normal_deformable_block.Block`解锁顶层法向，保留真3D粒子/邻域/内力 | PASS组件；64点，底16固定，其余48点xyz自由，xz/yz应力及内层z响应非零 |
| 7 | 法向材料功、壁功、势能、时间残差分账 | PASS代数；逐子区间独立复算；Rn未补热，Es不获自由能认证 |
| 8 | 标签交换/清零、ID排序、拒绝回滚与真实非零态重启 | PASS；15组件回归；两档实际预抬升checkpoint各续2步与原轨迹逐位一致 |
| 9 | 同几何/边界/材料的两档真实3D短程 | COMPLETE；block_dt5e9/block_dt2p5e9，120/240步共同600ns；外6.178/11.365s |
| 10 | 独立读回、全部VTK字段、按ID开离区间、同终时敏感性 | 读回PASS；362帧，121共同帧；速度/应力SENSITIVE_WARNING，不能记为全绿或收敛认证 |

最终新短回归43 PASS（19接口/单点、15三维组件、9独立参考/坏数据检测），
原选定26 PASS（12 contact、14旧变形块）。重复执行不重复计数。
四个完整单点基准、两个完整三维组件、两个实际checkpoint续算另列，不混算单元测试数。
既有材料闭循环/空间patch/零能模态等FAIL仍保留，本轮未重跑全部历史科学测试。

## 实际调用链和最小改动

`HostLeapfrog半漂移 → 3D NNPS/corrected_operators → material.update(L_before)`
`→ sigma_mid/伴随internal_force → moving_tool_rows与gap-derived Fn`
`→ 每粒子contact_kick（内力、法向、摩擦、底层约束在同一个kick）`
`→ 密度/速度kick → 半漂移 → 候选检查 → 接受/原子checkpoint`。
本构状态在回调中先更新；接触功使用实际子区间中点，不是旧L或整步端点平均假代替。
此组件等温，thermal NOT_SOLVED，没有热回调；未改material/mechanics/contact_kick/integrator旧源码。

新增helper、两个运行companion、两个独立reader、checkpoint检查器与三组回归。
没有用新求解器替换生产默认，没有GPU/并行性能开发，也没有增加依赖。

物理明确区别于旧fixture：kn=657500N/m；初始每顶点1N、重叠1/kn m，
初始16点总弹簧能1.216730038e-5J，初始材料应力零，是未平衡预压，不是论文初始化。
工具200ns后以+20m/s抬升到600ns，仅刚壁规定运动；自由点受内力/接触力推动。
没有改变质量、模量、几何尺度或时间单位。正常Fn会随变形减小，开离后不吸引。

## 独立参考与实测

单点m=1kg,k=4N/m，接触段delta=.5sin(2t)，u=-cos(2t)。
pi/4达到delta=.5m、U=.5J，pi/2开离，2s末静壁z=1.4292036732051034m、v=1m/s。
移动壁V=.2m/s：末z=1.8292036732051034m、v=1.2m/s，壁功.4J。
`inspect_normal_contact`只导入独立算术reader与通用VTK读取，不导入求解器生成参考。

| 静壁指标 | dt=.005s | dt=.0025s |
|---|---:|---:|
| 末位置绝对误差(m) | 6.015507e-6 | 1.095120e-6 |
| 末速度绝对误差(m/s) | 1.247521e-6 | 1.263068e-6 |
| 全轨迹最大物理K+U残差(J) | 1.249999e-5 | 3.124998e-6 |
| 开离采样区间(s) | [1.570,1.575] | [1.570,1.5725] |

两档均满足事前解析容限。末速度/能量误差没有单调改善，不能挑掉此结果或据两个dt声称阶数；
单边开离在步级采样，孤立弹簧的连续接触修正不变量也不是开离段的物理能量。
移动壁粗/细实际壁功.399999750496/.399999747386J，独立对照通过。

3D实际：所有16顶层ID在粗档[270,275]ns、细档[272.5,275]ns开离；
末态Fn和法向势能严格为零。全程邻域秩3，maxcond≈2.15227，
最大空间功共轭残差7.99e-15W；最大内层z速度粗/细.362354/.358365m/s。
实际ep全零：证明短时3D弹性动力学/接触，不证明塑性流动或沉积。

| 共同121帧最大差 | 实测 | 事前筛查 | 结论 |
|---|---:|---:|---|
| 位置(m) | 9.728908e-10 | 1e-7 | WITHIN_SCREEN |
| 速度(m/s) | .016479862 | .01 | SENSITIVE_WARNING |
| 全应力(Pa) | 194150.838 | 100000 | SENSITIVE_WARNING |
| 密度(kg/m3) | .00626337 | 1 | WITHIN_SCREEN |
| rawQ(J) | 1.051515e-10 | 1e-7 | WITHIN_SCREEN |
| 壁功(J) | 6.659613e-9 | 1e-7 | WITHIN_SCREEN |

比较已完成且有用，但不满足所有筛查。本轮不调容限、不再启动第三dt、不隐藏警告。
允许继续独立接口开发；不允许以此宣称连续体收敛/长期精度/生产许可。

## 能量与负数没有被藏掉

分别核算 `Rn=delta(U)+Wnormal_material-Wnormal_wall`，
`delta(K)+Q-Wwall-Wconstraint-Winternal-Wnormal_relative=0`。
粗/细累计Rn=-6.33135e-9/-1.20970e-9J；包含旧Es诊断的全局末残差
5.56615e-8/2.77394e-8J。没有触发事前能量screen，但Es仍NOT_CERTIFIED。
壁功约-1.129e-5J：抬升执行器从原先明确供给的初始弹簧能取回能量；负壁功本身不是负耗散。

细档原始trace有4个极小负Q，未abs/裁零，见saved_negative_friction_records.log。
最大负幅在step84/ID32，Q=-3.791011956846164e-23J，原逐项浮点界2.463291410915824e-22J；
step78/87/90的负值更小，也在各自浮点界内。独立NumPy逐项重算仅1项仍为负，
min=-3.791018785146777e-23J；`negative_raw_intervals`是重算计数，不是原始trace负数计数。
这不是过去明显负耗散的失败重现；旧失败不因此注销。

## 开发中保留的失败/修复

1. 源码审查发现单点构造器可接纳缺少accepted标志的直接state输入；补严格accepted/restart_allowed
   检查，加test19。未改任何动力学公式。之前point_cost_probe/static_dt005及日志均保留，
   原脚本副本probe_normal_contact_before_restart_guard.py保存原hash；最终采用guarded粗档。
   原/guarded静壁history及normal_intervals文件hash相同，不能重复计为两种物理结果。
2. 独立reader的手写1kg十进制测试fixture发生4.16e-17Ns减法舍入，超出组件审计1e-17Ns。
   首次9项中的1ERROR/3FAIL保留在reader_reference_tests.log。
   改为二进制可精确表示的独立有理数fixture：m1,dt1/8,Fn2.5,mu.2；
   Jt=-1/16,Jn=-5/16，deltaK=-3/256,Q=31/512,Wn=25/512。
   未修改reader容限、真实算例、力学；最终9PASS见reader_reference_tests_binary_fixture.log。

## 成本、限制及下一工作单元

单点完整运行内部.65–1.34s，外.99–1.69s；三维内部5.193/9.788s，外6.178/11.365s。
三维初始声学dt上界1.04907e-8s，接触上界7.16455e-8s；每步同时检查两者。
64点ParticleArray=26880bytes，不等于进程峰值RSS（NOT_MEASURED）。
600ns物理时间/外墙钟约9.71e-8、5.28e-8；不能外推为论文尺度成本。
单次外部限制单点60s、三维180s，均未到600s；运行日志COMMANDS.md可重跑到新目录。

限制：normal_snapshot/VTK的平面gap势能仅适用于这次远离孔/外缘的小块。
实际力仍由有限肩面contact_rows生成；原几何+10um位移上限保证不越孔边。
不得把此势能诊断泛化到孔壁/多面切换，需要按实际有限feature集合核算。
初始有弹簧储能、规定抬升制造程序、固定底层16点、等温且弹性，不是基板+方棒的AFSD。
论文PDF本轮未重读；没有从新图像/参数读取作新论文事实声明。

下一工作单元：局部轴向送料/执行器接口，先独立测约束反力、功及自由材料响应；
若进入旋转孔壁，先验证接触点速度重构与伴随力传递、有限feature势能，不直接拼生产运行。
本轮速度/应力敏感性作为后续物理精度门槛保留，不阻止做独立接口单测。
局部送料、孔壁力传递、热接入、有限应变改型、塑性沉积、空间细化、第三dt、三工况生产：NOT_RUN。
历史真实AFSD仍2us/191步，无持续沉积；不授权三工况生产。
