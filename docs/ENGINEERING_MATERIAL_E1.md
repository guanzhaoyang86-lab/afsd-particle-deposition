# E1：三维材料点多步塑性基准（2026-10-01）

## 预注册（运行前）

分类MATERIAL_POINT_BENCHMARK / SYNTHETIC_NOT_AA6061。
仅调用现有material.update，测试专用常屈服律Y=90e6Pa，不接生产配置；
G=30e9Pa、V=1e-6m³、T=293.15K，均为制造值或规定温度。
不是AA6061-T6标定，不是替换FKS的生产模型，不求解温度/接触/空间运动。
无质量/时间缩放：这是明确给定参数的独立制造问题，非改动AFSD物理。
所有运行参数、单位、来源、状态在configs/benchmarks/engineering_j2_cycle.json。

四个独立取向n=(1,0,0)、(0,1,0)、(0,0,1)、(1,2,2)/3；
A=1.5 n tensor n -0.5I，trace(A)=0，A:A=1.5。
规定材料点速度梯度L=e_dot*A，spin=0，等体积J=1。
e从0线性加载到0.004，再线性回程到0；每段1s，共同终时2s。
这是材料点加载输入，不是给自由粒子预设轨迹；F=exp(eA)只作为规定运动学记录。
不可由此称为三维空间SPH算例。

独立参考（不调用被测函数）：
q_signed=3G*e，首次屈服e_y=Y/(3G)=0.001；
加载饱和后q_signed=Y、ep=e-0.001。
峰值S=(2Y/3)A，ep=0.003；
卸载阶段q_signed=Y+3G*(e-0.004)，到e=0.002开始反向屈服；
此后q_signed=-Y，ep=0.003+(0.002-e)。
末态F=I，但有残余应力S=-(2Y/3)A、ep=0.005，这不是弹性闭循环。
累计塑性功=Y*ep*V=0.45J，末弹性剪切诊断=Y²V/(6G)=0.045J，
连续应力总功=0.495J。塑性功不是已解热量；温度恒定且热传导NOT_RUN。
返回算法的1.5G*(dep)^2*V单独记录，不加到塑性热/解析能量答案中。

时间网格：
- 每段N=16/32/64，对齐两次屈服转折，端点应力梯形功应为0.495J；
- 每段N=17/33/65，不对齐且N=1(mod4)。独立逐段梯形误差
  W_N=0.495-0.315/N² J；随步长减小趋于0.495，不用对齐采样伪称普遍精度。
  0.315=G*V*(.004)^2*(.28125+.375)，分别来自首次/反向屈服跨越区间。
所有网格到达同一终时，不因失败提前结束后直接比较。

预设容限：全应力张量相对Y最大误差1e-10；ep绝对1e-12；
每步塑性功与Y*dep*V/返回诊断公式误差<=1e-10J；
累积功/储能解析误差<=1e-10J；不对齐梯形功按独立W_N参考<=1e-10J。
未对齐相邻细化的连续功误差比<0.4，最细相对连续功误差<1e-3。
所有原始塑性/返回功非负（无abs/max/clipping修正），应力迹误差/Y<1e-10。
这些门槛针对合成材料点，不能直接套到实际自由表面/AFSD全局误差。

上游反例（不修上游、不把FAIL做xfail）：
- 固定hvi.VonMisesPlasticity2D.loop，纯xz剪切tau=2Y/sqrt(3)；
  独立3D q_trial=2Y，应径向返回q=Y，允许相对误差1e-10。
- basic.EnergyEquationWithStress.loop，m=rho_a=rho_b=1、p_a=p_b=2、VIJ=0；
  独立静止功率0，绝对容限1e-12（代码的SI制造输入）。
先执行真实loop再报告；不因类名或注释猜测。

本轮命令每入口timeout=60s、沿用隔离环境，不增加依赖。
旧边界RED、旧材料能量RED和闭循环诊断仍运行留档；
完整耦合回归/物理块加载/真实AFSD/热耦合/空间敏感性均NOT_RUN。
预注册/原AGENTS/原计划/HANDOFF保存在logs/engineering_e1_20261001/。

## 实际结果（2026-10-01，保留上述预注册）

最终入口test_engineering_material_e1.py：6 PASS，0.409s；
probe：6网格×4取向，全部到2s材料加载终时，PASS，内部计时0.416s。
这不是2s AFSD求解，不能外推沉积成本。
全历史应力张量最大相对误差1.68e-15，ep最大误差7.81e-18；
每点末态塑性功0.45J、剪切能诊断0.045J，残余应力/ep与独立分段解一致。

| 每段N | 应力梯形功/J（每点） | 连续功误差绝对值/J | 算法返回诊断/J（非热） |
|---|---:|---:|---:|
| 16 | 0.495 | 3.33e-16 | 0.05625 |
| 32 | 0.495 | 3.33e-16 | 0.028125 |
| 64 | 0.495 | 8.88e-16 | 0.0140625 |
| 17 | 0.493910034602 | 0.001089965398 | 0.051851211073 |
| 33 | 0.494710743802 | 0.000289256198 | 0.026983471074 |
| 65 | 0.494925443787 | 0.000074556213 | 0.013771597633 |

未对齐连续功误差比0.26538/0.25775；最细相对误差约0.0151%，仅为本基准数值误差，
不是论文误差。算法返回诊断随细化下降且单列，未加到物理热。
闭合的是规定等体积应变路径，塑性残余应力非零符合解析解；不和旧纯弹性闭循环混同。

| 实际入口 | 状态/证据（日志目录logs/engineering_e1_20261001） | 含义 |
|---|---|---|
| test_engineering_material_e1.py | 6 PASS；material_tests_final.log | 有限制造材料点范围 |
| probe_engineering_material.py | PASS；material_probe_final.log | 六档完整时序/配置和源码SHA |
| test_upstream_solid_limits.py | 2 FAIL/exit1；upstream_red.log | 上游限制真实重现，非xfail |
| test_spatial_probes.py --boundary-red | FAIL/exit1；old_boundary_red.log | 旧节点残差0.326157187666保持 |
| test_mechanics_review.py --hypoelastic-open-risk | FAIL/exit1；old_material_red.log | 旧Es非共轭相对缺陷0.9999985保持 |
| test_coupling_benchmarks.py --cycle | 1诊断PASS；old_elastic_cycle.log | 16/32/64/128步闭循环功均0.146111111111J；材料能量仍未通过 |
| 完整回归、受载块、空间细化、真实AFSD、热耦合、GPU/生产 | NOT_RUN | 本轮未启动 |

固定上游限制实测：
1. hvi.VonMisesPlasticity2D.loop：xz纯剪切q_trial=180MPa，返回仍180MPa，
   独立3D屈服参考90MPa，相对误差1.0。impact3d.py实际采用该类，不能照搬此塑性。
2. basic.EnergyEquationWithStress.loop：静止正压对原始d_ae=2，独立功率参考0；
   源码计算vijdotdwij却未乘入此项累加。原始2缺正确单位，不能当作有效能量率。
结论仅针对所测固定版本/路径，不意味着PySPH设施不可用。
现有AFSD没有因本轮而切换上述函数，上游/许可证未改。

初版5测试PASS；独立审查指出固定解析公式不能泛用N=18等网格。
新增拒绝未注册网格/制造参数变化的接口测试后6 PASS，不改物理/容限。
初版日志/输出保留；probe_before_oracle_guard.py快照对应初版报告source SHA。
最终6档数值与初版一致，没有因失败调参。

结论：进入E2独立受载小块设计/短基准，不直接进入大变形热耦合或AFSD。
现有更新器无需因放弃逐式复现而整套推倒；尚未选择生产替代材料模型。
Jaumann变体积能量、自由表面空间误差及真实接触耦合仍须量化。
目前层级=独立制造基准验证；AA6061标定、定性论文对照、定量复现均未达到。

## 重跑与输出接口

精确命令/exit见logs/engineering_e1_20261001/COMMANDS.md。
使用隔离.envs/pysph/python.exe -B；每入口subprocess.run(...,timeout=60)。
probe --output必须给新目录；最终outputs/engineering_e1_20261001_final。
无新增依赖、无src/生产参数变更。

NPZ是材料基准而非AFSD粒子schema：
time_s[s]、prescribed_log_strain[1]、prescribed_F[1]是规定加载；
stress_Pa[Pa]、plastic_strain[1]、equivalent_strain_rate_per_s[1/s]来自实际更新器；
cumulative_stress_work_J/cumulative_plastic_work_J/elastic_shear_diagnostic_J及增量[J]；
cumulative_algorithmic_return_J及增量[J]单列非物理热；
reference_stress_Pa/reference_plastic_strain是独立参考而非求解字段。
应力轴(time,orientation,3,3)；标量(time,orientation)，增量无初始行。
没有粒子ID/空间坐标/接触力/已解温度，因为不是粒子空间求解。
report.json含配置/源码与6个NPZ的字节SHA256、所有检查及未解项。
artifact_check.log：6个NPZ实际重读、独立末态应力/ep、全部字节hash PASS；
初末版六个NPZ字节完全一致，初版源码快照hash对应其报告。
最终冻结ENGINEERING_E1_FREEZE.json；文件完整性PASS不等于科学验收。


