# 旋转肩面—三维可变形小块：五步推进报告

2026-10-02；后端仍PySPH 1cd3dea9cfaac7d9f1bda2a2491ce4b42a31dbe5。
工程近似方法组件，不是AFSD沉积或论文完成声明。旧物理、配置、失败、上游均未改。
新文件为伴随试验入口/独立读回/测试/配置/合同；不重搭后端或界面。

## 1. 冻结物理与真实时序

事前参数/容限见ROTATING_DEFORMABLE_CONTRACT.md与
configs/benchmarks/rotating_deformable_block.json；每项SI/来源/状态已记录。
64点(4x4x4)、2mm立方实体单元域、dx=.5mm、G=26.3GPa、rho0=2700kg/m3，
每点3.375e-7kg、总2.16e-5kg；未缩放质量/模量/时间。
底层16点xyz制造夹具；顶层16点只固定z，xy受接触力；内层32点xyz自由。
刚壁omega=20rad/s，开放方孔仍保留，本块位于孔外肩面下方；不是论文工况。
Fn=1N/顶点规定、mu=.2、法向gap=.2mm(检测overlap=.05mm)，非penalty/开离解。
所谓4MPa仅16N除以4mm²的名义接触载荷：顶层法向载荷被z夹具局部平衡，
**没有求解块体4MPa预压或工具下压**；本次内力传递由切向摩擦驱动。

实际调用映射(全部原设施复用)：

| 阶段 | 代码位置 | 独立检查 |
| --- | --- | --- |
| 半漂移 | integrator.HostLeapfrog / 固定上游LeapFrog | 保存x_before/x_mid/x_after |
| 真三维邻域与校正L | mechanics.corrected_operators | 非共面rank3、字面核导数、完整半径邻域、仿射A |
| Jaumann/J2/FKS | material.update，L_before | 非零初始应力客观旋转；旧材料回归 |
| 中点应力及内力 | pressure(rho_mid)、internal_force(...,B) | 独立配对力、总内力/矩、伴随功 |
| 旋转接触/约束 | contact_rows、原kick | 首步解析解、逐子区间signedQ与实际reaction |
| 速度/密度kick与半漂移 | 原HostLeapfrog | 保存状态逐步独立读回 |

材料在接触加速之前；不伪称按contact→material次序。另算真实各子区间中点速度的
dt加权平均，并用同一冻结i/j/g/B算L_midpoint，只诊断、不改变本构。
material_old_work和实际Winternal对应的midpoint_stress_work分开记录，不能强行相等。

## 2. 实际接入，未求字段不伪造

scripts/probe_rotating_deformable_block.py的Block提供advance/record/stable_dt/
checkpoint/from_checkpoint/vtk。源哈希02875c7fa0319e251d3362981f37b3d8147e83ae1f28b62497e592bd96aab11e。
真正调用原3D LinkedListNNPS、原校正SPH伴随内力与原材料更新，没有复制二维结果。
自由点无指定轨迹、无刚体拖拽；标签只输出，不入受力。底部夹具不等于全块刚化。
原FKS的b1量纲/低速正则化等问题保持APPROXIMATE；没有以JC/流体/新本构替换。
293.15K是规定输入，thermal NOT_SOLVED；Q原始值仅入耗散账，不变成已求温度。
checkpoint仅当前固定组件重启，生产false；候选失败全属性/累计账回滚，
REJECTED单独保存。被拒绝非有限诊断显式标记，不变零、不写入accepted末态。

## 3. 独立短回归与旧接口兼容

tests/review/test_rotating_deformable_block.py：最终14PASS。
包括完整SI/质量/角色、首步独立冲量和Q、三维仿射/非仿射伴随功、带非零应力旋转、
内层出平面响应、标签交换/清零、非零材料态checkpoint后续2步、重排按ID对应、
注入失败全状态回滚、VTK只含已求字段、声学界、拒绝REJECTED重启。
标签清零不是删除属性；不报告删除schema属性回归通过。
首步答案直接使用m=3.375e-7kg、Jt=.2dt*wall/|wall|、Q=Jt·wall-Jt²/(2m)，
不调用被测kick生成参考。没有当前自由表面块体的连续体解析轨迹。

旧选定兼容：test_contact_events.py 12PASS(反向/静动壁/多面/执行器/功等)，
test_mechanics_review.py默认14PASS(3D/FKS/客观率/ID等)。旧风险组本轮NOT_RUN，
以前raw Table1、Jaumann闭循环、边界精度FAIL全部原样保留，不据新PASS撤销。
初版13PASS；随后只补拒绝诊断日志完整性，重跑13PASS；另加材料/速度已更新后
acceptance失败的回滚回归，第14项也PASS。旧日志都保存，不重复计数。
独立reader合成/篡改回归12PASS，含signedQ/冲量/反力/舍入界/时间层、
邻居体积权重、独立首步值、参数漂移/假温度拒绝、无求解器导入。
本轮新26PASS(14+12)，选定旧26PASS(12+14)，不把41/81帧计为新增单元测试。

## 4. 实测成本与唯一粗档

所有命令见logs/rotating_deformable_20261002/COMMANDS.md，单线程，外限120s/内110s。
两步成本探针外.461538s、含输出内.0797125s，纯计算平均.02231035s/步；
声速6195.9647m/s，初始声学界1.04906989e-8s。这不是空间算子谱稳定证明。
计算dt=5e-9s/40步，细档2.5e-9s/80步，共同2e-7s；成本事前保守估计20.94/36.88s。
粒子属性26880bytes，非峰值进程RSS；峰值RSS NOT_MEASURED。未以此推算论文全尺度成本。

粗档实际一次完成40步/200ns，无重试：外2.070255s，含输出内1.4872354s。
末态速度最大.09793556m/s、密度相对变动2.2291971e-6；内层自由点w最大.018658m/s。
max|sigma_xz|=32978.3Pa，max|sigma_yz|=440242.9Pa；不是只靠三维点云宣称3D。
ep全0：这是已接入材料更新器的**三维弹性响应**，不是已展示屈服/塑性金属流动。
工具仅转4e-6rad，不是宏观旋转过程或沉积完成。

粗档Q=9.49019362848e-8J、Wwall=1.28199238143e-7J、K=2.41276992938e-8J，
Es=9.11211291404e-9J(诊断)，Pplastic=0；全局诊断残差=-5.74896503638e-11J。
虽在事前短程筛查内，旧有限应变材料能量问题仍未解。
当前末步1068邻居对、cond最大2.152132；实际材料功分裂完整保存，未补热。

## 5. 独立读回与同终时对照

粗档完整后实际运行独立reader，PASS_DEFORMABLE_READBACK_NOT_AFSD，外2.190215s。
保存prefine_gate后才唯一启动细档；没有调参、没有物理失败重试。
细档80步同到200ns，外3.714041s、含输出内2.9466191s；独立reader也PASS，外4.690264s。
比较命令只读既有结果、不重跑动力，外7.651396s，WITHIN_SCREEN。

| 共同41时刻的最大分量绝对差 | 实测 | 事前筛查 |
| --- | --- | --- |
| 位置 | 3.49470e-11m | 1e-7m |
| 速度 | 4.66523e-4m/s | 1e-2m/s |
| 全应力 | 4530.068Pa | 1e5Pa |
| 密度 | 1.23457e-4kg/m3 | 1kg/m3 |
| 累计rawQ | 9.23421e-11J | 1e-8J |
| 累计壁功 | 9.51349e-14J | 1e-8J |

这是同几何/边界/参数/源代码、同终时的两档时间筛查，不是实测收敛阶、
连续体精确解或相对实验误差。空间加密、第三时间档、完整AFSD仍NOT_RUN。

独立reader从文件重算字面核导数/完整半径邻域、rank3、moment逆、伴随力与功、
全部粒子与子区间冲量、Coulomb方向/容量、原signedQ/舍入界、实际dP/dL、
Wwall=omega*接触转矩冲量，以及全部122个VTK帧/持久ID/材料字段。
两档maxcond约2.152133，最大功共轭残差1.388e-16/2.776e-16W；
最大内力模0.215074/0.219225N、内层|w|max=.018658/.018753m/s；
rawQ无负值，最小子区间4.53233e-11/2.27202e-11J，未abs/裁零。
细档ep也全0；max|sigma_xz|=33386.3Pa、max|sigma_yz|=442523.3Pa。
两档末态全局Es诊断残差分别-5.74897e-11/-2.91053e-11J，保持原符号；
全程最大绝对残差5.97650e-11/3.02470e-11J。短程筛查内不等于材料能量获认证。
reader的checkpoint检查仅内容hash；真正非零态重启续算验证在第8/9项短回归。
它不导入求解器/本构/接触；不提供当前块体完整连续体解。

完整输出：outputs/rotating_deformable_20261002/cost_probe、dt5e9、dt2p5e9、
readback_dt5e9.json、readback_dt2p5e9.json、comparison.json及NOT_RUN.json。
新reader编写曾遇并行代理连接中断，接管完成后测试/数据读取均一次通过；
不是物理失败或恢复数据覆盖。真实全AFSD历史仍2us/191步，无持续沉积。

## 下一项与保留缺项

不再反复打磨四质点或本轮固定法向夹具。下一阶段重点是法向承载/开离：
先用独立压入—卸载接触功基准验证正常力与间隙/接触势能，再以命名companion
解除顶层z夹具接回小块，记录块体压缩/回弹及原始能量差，不能默认为已受法向预压。
范围主要涉及contact_kick.domain_kick、WallContact、几何与新适配层；原FAIL保留。
旋转孔壁面点重构须配转置力传递、局部轴向送料、热耦合、FKS材料细节、
有限应变自由能、空间分辨率和论文对照都仍缺项。禁止直接启动三工况生产。
