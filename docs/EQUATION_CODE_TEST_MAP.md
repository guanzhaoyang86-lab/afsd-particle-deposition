# 方程 / 边界 → 实现 → 测试（2026-09-24）

2026-09-30耦合补充（既有方程/物理参数不被静默替换）：

| 接口/数值关系 | 实现 | 独立证据 |
|---|---|---|
| 原实际kick完整时序 | solver.evaluate + 上游LeapFrog stage1/2 | COUPLING_FAILURE_REVIEW第1节；31步逐位一致观察 |
| 多力/多面+约束迁移率C_ab=R_a P R_b/m | contact_kick.kick/cone_solve，显式APPROXIMATE变体 | test_contact_events解析反向/多面执行器/壁功；新候选RED仍FAIL |
| 每接受子区间Q=-Jt dot u_mid、反力功 | contact_kick.interval/domain_kick | 解析标量解、独立向量极坐标积分、RK4转向 |
| 真实material.update闭弹性路径 | tests/review/test_coupling_benchmarks.py | 非零极限功0.146111111111J；不是新的free-energy模型 |
| 常应力边界/零模/畸变自由表面 | test_spatial_probes.py，mechanics.py未改 | 弱仿射PASS，节点牵引patch FAIL；单云仅6刚体零模 |
| REJECTED候选不得重启/末态输出 | output.restore/checkpoint + application | 拒绝重启/中断原子写入回归、失败run输出完整性 |

本轮独立审查覆盖：SCIENTIFIC_REVIEW.md / VALIDATION_CONTRACT.md；下面历史21项
后端测试不能授予生产资格。新增映射：

| 修复/假设 | 当前路径 | 独立证据 |
|---|---|---|
| 完整性/离散功共轭（明确偏离原Table1） | mechanics.corrected_operators/internal_force；solver.evaluate | tests/review/test_mechanics_review.py；原始2反例FAIL、修正基准PASS |
| 实际kick中点摩擦热 | contact.evaluate；solver.evaluate | test_contact_review.py；孤立功热PASS，多力耦合4档FAIL |
| 深陷工具拒绝、开方孔 | geometry.assert_outside_tool/contact_surfaces | test_contact_review.py；外圆/顶面近接触仍BLOCKED |
| 固定物理执行器/面罚刚度、热库限dt | particles/config_adapter/solver | test_resolution_review.py（2项PASS） |
| 预应力压缩自由能 | material.update；solver.energy | UnresolvedThermodynamics反例FAIL，未补假热 |
| 观测算子而非实验面积 | afsd_validation/measurement.py | test_measurement_review.py（9项PASS，仅代理） |
| 同终时敏感性/冻结 | application --end-time；summarize_scientific_review.py；verify_scientific_freeze.py | outputs/review/sensitivity_summary.json，N09/N10 BLOCKED |

实施前建立；每一关的实际结果写入 logs/physics/ 和 docs/PHYSICS_STATUS.md。
唯一后端 PySPH 1cd3dea9cfaac7d9f1bda2a2491ce4b42a31dbe5。

| 论文/来源 | 计划代码 src/afsd_backend/pysph/ | 必须通过的测试 |
|---|---|---|
| P0 p5 Eqs2–3 hyperbolic spline | kernel.py | 3D 归一化、分段导数、方向、零距 |
| P0 p5 Table1 质量/动量/梯度 | mechanics.py + PySPH LinkedListNNPS | 3D brute-force 邻域、仿射场、成对动量守恒 |
| P0 Table1 Jaumann/应力分解 | material.py | 非零初应力刚体旋转、三向剪切、零迹 |
| P0 p6 Eqs4–8、p14 Table2 FKS | material.py | 独立公式、1/c2 括号、温度/率有效域、3D J2 返回残差 |
| P44 pressure closure | mechanics.py | 密度/压力一致性、体积能导数 |
| P0 Table1 热传导 / P44 热源 | thermal.py | 两体热平衡、解析模态、塑性功热转换 |
| PySPH RigidBodyWallCollision / P44 接触 | contact.py | 3D 法向、滑动耗散、功热账本、分离、方孔转动 |
| P0 方棒/基板变形、平肩/旋转方孔刚体 | geometry.py | 开孔、层厚解析、仅执行器规定速度 |
| P0 leapfrog forward difference | integrator.py + upstream LeapFrogIntegrator | 常加速度、二阶位置收敛、耦合能量检查 |
| P0 p6 Fig5 局部3–4mm标签 | particles.py | 局部/全棒独立、材料不分支、ID重排/注入 |
| 工具平移系与 AS/RS | geometry.py | 旋向反转、平移变换 |
| checkpoint / output | output.py | 重启连续性、VTK/NPZ回读、不可变身份 |
| 小规模三维耦合 | application.py | dt/墙钟上限、质量/功/热/能量/密度/穿透/NaN 门限 |

上述映射在实施前建立。现有路径已落地；实际测试名/判据和通过范围见
tests/backend/test_gate_{a,b,c,d}.py 与 docs/PHYSICS_STATUS.md。
21项后端测试已运行通过，不代表整张表中所有目标都完成。
P44压力行实际采用PARAMETER_SOURCES.md所注明的bulk-EOS变体，非逐字P44声速系数。
接触为有限解析面/无记忆滑动DEBUG；三角网格/静粘滑未实现。
无限入口只测试append身份接口；当前有限棒无外部增质量。定性/定量论文验证未进行。
