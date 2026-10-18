# AFSD 独立科学审查 — 2026-09-24

当前更新：METHOD_CLASS_REVIEW.md。固定432接触锥类及439局部kick通过，
116项111PASS/5FAIL；边界固定域探针、本构非零闭循环仍不满足完整方法验收。
本轮无新AFSD试跑，以下为历史审查，不能据其旧439状态覆盖新证据。

最新限定修复见CONTACT454_REPAIR_REPORT.md：454已通过独立参考；
真实试跑接受26步后候选27/439双面锥失败，107项101PASS/6FAIL。
原材料/边界风险不变，未获生产授权。以下接触摘要亦为第一轮历史。

2026-09-30授权接触诊断修订见 `COUPLING_FAILURE_REVIEW.md`。
原首发已定位到32/485/BORE_Y_MINUS；原证据保留。新局部事件变体基准通过，
但真实固定问题第3候选在454/BORE_X_MINUS迭代失败，仍BLOCKED；没有生产许可。
以下为2026-09-24历史审查，不覆盖新报告，也不把旧冻结自动改成新源码。

## 结论与证据等级

**不允许进入三工况生产运行。当前级别：真正三维的 APPROXIMATE 方法开发/单元基准；
尚未获得稳定耦合试算、宏观定性对照或定量论文验证。**
原交接中“通过”的日志保留为历史证据，不能覆盖本轮新发现。
唯一后端仍为 PySPH `1cd3dea9cfaac7d9f1bda2a2491ce4b42a31dbe5`，
上游未修改、许可证保留；未重搭项目、切换后端或制作界面。

先冻结 `VALIDATION_CONTRACT.md`，再测试。预注册副本：
`logs/review/contract_v1_pretest.md`，SHA256
`e5a39ffec3d333c512afd1e25a3151f72743a6a0e5491d6a6858a6e88908e091`。
所有容限是数值门槛而非实验误差；失败后未放宽、未拟合论文面积比。
源代码没有顶层Git仓库，故用文件hash和 `logs/review/before_source/` 比较，不伪造项目commit。

## 独立检查矩阵

下列PASS严格限于所列基准；完整38项独立测试实际为35 PASS、3 FAIL，exit1。
三个FAIL不是隐藏skip：两个保留原始Table1离散的反例，一个是修正路径仍存在的本构能量反例。
总日志 `logs/review/independent_suite_final.log`；旧骨架19/19、后端21/21在本轮重跑通过，
分别为 `scaffold_final.log`、`legacy_backend_final.log`，它们不能推翻独立失败。

| 检查 | 状态与证据 | 修复 / 剩余风险 |
|---|---|---|
| 3D核及真实邻域 N01–02 | PASS：独立分段积分/差分、非共面暴力邻域、z向及三个剪切面；tests/review/test_mechanics_review.py | PySPH LinkedListNNPS(dim=3)，3D位置/速度/3×3应力；无2D厚度复制 |
| 材料点Jaumann非零应力 N03 | PASS：非零初始偏应力，独立轴角sin/cos参考，不调用被测旋转函数 | 只证明材料点及已修正空间旋转基准；不等于有限变形热力学正确 |
| 原空间梯度/功 N03/08 | FAIL（原式）；mechanics_raw_red.log | 刚体转动产生伪对称梯度，相对0.4208896；异密度应力功缺口9.76246% |
| 修正空间梯度/伴随力 | PASS限定基准；mechanics_corrected_green.log、mechanics_integrated_final.log | R04：完整性矩阵修正及功共轭力；全3D旋转、线/角动量、功检验通过。明确APPROXIMATE，非逐字Table1 |
| FKS/J2 N04 | PASS公式/数值；BLOCKED物理标定 | 独立Decimal/math式及独立二分返回映射，三向剪切；b1、参考率/温度未解决，不能称paper配置 |
| 压力体积能 | PASS：独立EOS原函数；同上 | K线性EOS是声明的近似；非P44全套压力系数 |
| 本构储能一致性 | FAIL：hypoelastic_open_risk.log；UnresolvedThermodynamics | 预应力各向同性压缩产生非零储剪切能变化但S:D=0；不是时间步误差。未加补偿能/假热源 |
| 导热、塑性耗散 N05 | PASS限定两体及材料点；contact_post_fix_corotation.log、独立suite | 热到冷、成对守恒；q_new*dep为热源，返回数值耗散另列；真实温变参数/热边界BLOCKED |
| 离散接触摩擦功 N06 | 原实现FAIL，单接触修复后PASS；contact_discrete_defect_pre_fix.log、contact_post_fix_corotation.log | R01：冲量上限时旧算法热2J而实际动能损失1J；改实际速度kick中点功。多力耦合反向滑移仍FAIL，见四档试算 |
| 方孔/深陷工具/完整接触 | 开孔与深陷拒绝PASS；完整边界FAIL/BLOCKED | R02深穿透fail-closed；外圆/顶面粒子半径重叠却零接触力，contact_unsupported_faces.log；静粘着、棱角、完整外表面未实现 |
| 物理边界/缩放 N08 | PASS代码审查及固定边界基准；tests/review/test_resolution_review.py | 基板/棒均变形；仅上部执行器x/z规定，y自由。无自由材料速度场、事后搬运、标签力、截温/隐藏质量时间缩放 |
| dt稳定性 | PASS显式热库限步回归；BLOCKED全耦合稳定性 | R05加接触导热行和限制；波速/罚接触限制已用真实物性，但修正梯度谱半径、拉伸失稳和粘滑仍无充分界 |
| 标签/ID/重启 N07 | PASS：交换及删除两种标签后物理逐位一致；非零应力/温度重排和重启 | 全部AA6061-T6，局部与全棒分开；ID int32溢出拒绝。无限入口动力学未实现，不能把append身份测试当入口验证 |
| 时空敏感性 N09–10 | 已运行，四档全部FAIL；同终时比较BLOCKED | 无档达到预定2µs；不比较不同终时、不宣称收敛，详细下表 |
| 观测算子/坐标 | PASS防误用9项；measurement_protocol.log | 仅R_count/R_volume代理；未知ROI拒绝、空分母null；实验面积重建BLOCKED，见VALIDATION_CONTRACT |
| Fig8/10/11/12 | 页图核对PASS；仿真形貌/轨迹对照NOT_RUN；定量验证BLOCKED | 未产生足够物理时长/沉积；无原始图像/XRCT、精确ROI或参考CSV |
| 三档长时间生产、分布式/加速性能 | NOT_RUN（未授权且科学门槛未过） | 没有长任务/付费任务；当前冻结只审查starved近似案例的离散敏感性 |

## 修复和独立依据

- R01 `contact.py:WallContact.evaluate`、`solver.py:evaluate`：
  力仍按原接触律；热改为 `-Ft·(v_actual_mid-v_wall)*dt`，
  其中中点使用全部力与执行器约束后的实际kick。恒壁孤立冲量满足
  ΔK=J·(v_old+v_new)/2；这是独立能量恒等式，不是被测函数生成参考。
  负耗散显式报错，绝不取abs/max(0)来“修复”。当前多接触/内力耦合可能反向滑移，
  需要一致约束/粘滑求解而非更大摩擦、更小温域或裁剪热量。
- R02 `geometry.py:assert_outside_tool`：工具实体内部不能因有限面mask为空而穿过。
  旋转方孔始终开放。该防护不是完整外表面/棱角接触；中心尚未进实体的球面接近
  外圆/顶面仍有盲区。未把工具改成封闭实心柱。
- R03 `particles.py`、`config_adapter.py`、`solver.py` 和配置：
  原细化把执行器体积减半、单位面积罚刚度加倍，证据
  `resolution_boundary_red.log`；改固定物理深度0.5mm、
  kn=0.05G dx²/(0.5mm)，故kn/dx²=2.63e12 N/m³。
  `resolution_boundary_green.log` 和最终suite通过；总质量两网格同为0.0002217375kg。
- R04 `mechanics.py:corrected_operators/internal_force`：
  V=m/rho，M_i=ΣV_j(r_ji⊗grad_i W)，B_i=M_i^-1，L_i=L_raw B_i；
  F_ij=V_i V_j(σ_i B_i^T+σ_j B_j^T)grad_i W。
  独立仿射/任意轴旋转检验和Σv·F+ΣVσ:L=0验证，不只检查力反对称。
  奇异/条件数>1e6邻域拒绝，不伪逆兜底；导热仍用原径向梯度以保持符号。
  原始operators保留为audit。修正不解决下面的本构自由能缺陷。
- R05 `stable_dt/evaluate`：显式热扩散行和加入每活跃面的hc dx²，
  dt<=0.25mCp/Σconductance，具有SI秒量纲；独立大hc探针回归捕获漏项。
- 输出/试验接口：`application.py` 新增 `--end-time`，实际最后dt取剩余时间
  （非时间缩放）；未到终时为INCOMPLETE/FAILED。运行分类改APPROXIMATE，
  保留契约hash、source/config hash及 `production_authorized=false`。
  原日志/初始失败命令均保存，没有覆写到“全绿”。

未修复反例：G=26.3GPa、m=1kg、rho=2700kg/m³、
S=diag(2,-1,-1)MPa、D=-I/s时，σ偏应力功S:D=0；当前constant-G
Jaumann保持S不变但体积减小，E_s=V S:S/(4G)，其导数=-3E_s=-0.063371356J/s。
独立小差分得到-0.063371261J/s，不能以小dt收敛消除；应重审本构自由能/体积耦合
及能量字段解释，不能在账本里添加同量“修正热”。

## PDF与来源复核

本轮重新渲染、查看用户PDF第5/6/14页及Fig8/10/11/12所在第9/10/11/12页，
另看第7页实验测量说明；图像在 `logs/review/paper/`。不是沿用上一模型的通过声明。
方程、FKS完整括号、b1量纲矛盾、侧带、Fig12的7.63mm疑点详见契约第5节。
`PARAMETER_SOURCES.md` 中P40/P41只取得落地页，非全文；本轮未取得新全文。
P44可读，但不能用它的旧FKS参数填充P0新标定。仍缺k(T)/Cp(T)、真实热分配、
摩擦粘滑、启动/停留/抬升程序；全部近似显式保留，不默认为零。

无二维结果沿厚度复制、全基板刚化、封死方孔、自由粒子预设旋转/径向漂移或侧偏
代码路径。旋转从接触力传给金属，热启动/初始过盈/预结合SPH界面是明示假设。
局部3–4mm标记和全棒诊断不能互换。AS/RS由(-omega*y)*Vx符号确定；
目标富集趋势未写入受力。当前正旋向配置的manifest正确，但application的说明字符串
硬写正旋向/corrected；反转旋向或选择raw_audit前须先修动态元数据，不能沿用该说明。

## 已执行的敏感性与成本

每档上限600s、预定同终时2.0e-6s、最多10000步；只改CFL或dx，
物理几何/温度/进给/转速/执行器厚度/单位面积接触刚度保持一致，h=1.3dx。
缩小横向几何为DEBUG来源，不能称论文尺寸；层高1.5mm分别为3/6个粒距。
所有运行均在 `contact.py` 负实际中点摩擦耗散处exit1。

| 输出目录 outputs/review/ | N / CFL / dx(mm) | 已接受步 / 终时(s) | 中位步耗时(s) / 总墙钟(s) | 峰值内存(MiB) |
|---|---|---|---|---|
| coarse_dt | 657 / .1 / .5 | 31 / 3.2519643e-7 | .0295674 / 1.08083 | 57.75 |
| dt_half | 657 / .05 / .5 | 1 / 5.2453475e-9 | .030446 / .05715 | 56.82 |
| dt_quarter | 657 / .025 / .5 | 2 / 5.2453421e-9 | .0258242 / .08326 | 57.63 |
| space_half | 5256 / .1 / .25 | 3 / 1.5735911e-8 | .2188855 / .88409 | 94.83 |

粗档dt约1.049e-8s；细档空间约5.245e-9s。较小dt没有消除反向滑移失败。
失败发生在不同终时，N09/N10的误差/下降趋势/精细两档比值全部null/BLOCKED，
绝不把“更早失败”称空间发散率或给出虚假的收敛阶。未继续调参。
所有已接受步的塑性应变仍为0，故这些尝试不证明AFSD塑性沉积或混合。

已记录前缀的|能量残差|/(0.01*交换尺度+1e-7J)最大值依次为
0.307598、1.807e-8、9.491e-7、0.00045246，均小于1；质量/有限值前缀也通过。
这不是完整N08通过：失败的部分kick没有合法能量状态，同终时仍未达到。

每目录有performance.json、diagnostics.jsonl、resolved_config/source_hashes、
manifest/series和checkpoint；对应日志 `logs/review/{coarse_dt,dt_half,dt_quarter,space_half}.log`。
为保留失败现场，未把部分完成步伪装为成功；本轮output-every=10000，因此四个有效
checkpoint均只到step0，31/1/2/3个接受步仅有标量诊断，失败中间态不写成有效场。
`readback_failed_runs.log` 的PASS仅指这些初态VTK/NPZ文件完整，运行状态仍FAILED。
`sensitivity_summary.json` 是只读汇总，不生成不存在的末态场或实验数据。

外推不是实测论文成本：按80,312粒子、名义dx=1mm、线性N步耗时和同CFL波速dt，
粗/细空间样本估计每1秒物理时间约4.43–4.79万墙钟小时（约5年单进程）；
半/四分之一CFL更贵，约9.85/16.72万小时。此层厚仅1.5粒距，
不满足本项目最小3粒距调试门槛；解析层厚的实际成本更高且未知。
初始化/IO/并行扩展、稳定接触、邻域退化都未进入模型；短失败步数样本不支持精确预算。
没有可可信执行的三档生产成本承诺或付费估价。

## 复跑与下一步

在项目根目录，用现有隔离环境（没有安装新依赖）：

~~~powershell
.\.envs\pysph\python.exe -B scripts\verify_scientific_freeze.py
.\.envs\pysph\python.exe -B scripts\run_scientific_tests.py
.\.envs\pysph\python.exe -B scripts\run_physics_tests.py --gate all
.\.envs\pysph\python.exe -B scripts\summarize_scientific_review.py
~~~

当前独立suite应exit1（保留3个科学反例），不是部署失败。冻结校验/旧后端测试应exit0。
汇总首次实际运行exit0；现有sensitivity_summary.json存在时会拒绝覆盖并非零退出，
这是证据保护，不应删除原文件来重跑。重检已有汇总用冻结校验即可。
复现耦合阻塞用下面命令；输出必须为未存在目录，不覆盖本轮记录：

~~~powershell
.\.envs\pysph\python.exe -B scripts\run_afsd_debug.py --case configs\debug\starved.json --output outputs\review\recheck_coarse --steps 10000 --end-time 0.000002 --output-every 1 --wall-seconds 600
~~~

该建议命令尚未执行，当前预期负摩擦耗散exit1；output-every=1只增加有效接受步保存，
不改变物理或提供生产授权。另三档将case换configs/review/*.json并用独立输出目录。

下一轮按顺序：
1. contact.py/solver.py：一致的多面法向/切向约束与实际中点功，独立组合内力、
   移动壁、执行器及粘滑基准；完整外圆/顶面/棱角及接触对力矩检验。不得裁负热。
2. material.py/solver.energy/output.py：建立与Jaumann/压力/体积变化一致的储能/热力学
   声明与积分；带预应力压缩、循环加载和旋转的独立能量基准。不得用账本补差。
3. mechanics.py/stable_dt：修正矩阵谱稳定界、自由表面/拉伸稳定化及孤立粒子处理，
   每项说明与论文差异；修复后重新通过同终时dt/空间门槛，再讨论更长物理时长。
4. 补材料/边界来源；实现守恒入口而非位置搬运；编译/加速现有唯一后端前后逐字段回归，
   再测真实吞吐。先修通用metadata，不让反向旋转输出错误AS/RS说明。
5. 精确实验ROI/原始图像和采样算子到位后才评Fig8；局部/全棒、校准/验证严格分开。
   可视化仅消费签名核对后的原始字段，不得改变冻结物理。

本轮没有获得“完整复现”标签；资料不足允许继续APPROXIMATE方法开发，但生产禁入。
