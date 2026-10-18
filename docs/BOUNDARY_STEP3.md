# 边界第三步：固定边界的单内部节点扰动（2026-10-01）

## 预注册：实验设计与验收

本轮仅静态制造解/独立原型，不是粒子时间演化或AFSD轨迹。
src/、configs/、旧脚本/测试/失败/冻结不改，仅新增本轮文件并更新HANDOFF。
不换后端/本构，不做性能优化，不运行AFSD。沿用隔离Python/NumPy/PySPH。

固定 [-1,1]^3 m、125点、dx=0.5m、h=0.65m、rho=1kg/m³；
沿用STEP2体积V和几何边界面积E，质量、密度、边界位置/法向/牵引均不变。
只将原中心(0,0,0)m移动到(0.025,-0.015,0.010)m，即各轴(5%,-3%,2%)dx。
不换随机种子、不做幅值扫描、不试另一扰动直到通过。
这是固定积分权重能否用于不规则采样的反例，不是一个动态可接受物理状态。

### 独立推导：先查是否可能有解

取ell=1m，P=[1,x/ell,y/ell,z/ell]，d_kP的对应坐标列为1/ell。
三个4x4必要矩条件为

    C_k = P^T H d_kP + (d_kP)^T H P - P^T E_k P = 0.

由 Q+Q^T=E 和 QP=H dP 直接推出，与SPH核、邻域、求解器迭代数无关。
本轮P无量纲，C所有分量单位m²。容限max|C|<=1e-10m²，
与上一轮ell=1数值门槛一致，仅用于双精度代数恒等式，不是实验误差。

原中心体积0.125m³。边界完全不动，体积一阶矩改变为
mu=(0.003125,-0.001875,0.00125)m^4；质心误差mu/8。
故C的常数行/列变化为0，坐标块独立参考为

    (C_k)_{lm} = (mu_l delta_mk + delta_lk mu_m) / ell².

例如C_x的xx元素为0.00625m²，明显大于1e-10。
这是运行前的解析预期：冻结V/E且按上述位置取P时，不可能同时满足
严格的SBP和仿射再现；改进SVD求解精度不能消除此矛盾。
并非证明所有不规则SPH都不可能，只是否定这组冻结积分/边界权重的直接推广。

实际路径：固定fixture -> 独立矩检查 -> 原corrected矩阵提取/条件数诊断
-> 新独立checked_project入口。C失败必须在调用STEP2 project之前抛出明确异常。
不构造/保存“貌似通过”的最小二乘折衷算子，不放宽门槛。
原规则点云作为唯一对照，照旧调用STEP2 project；生产不导入新入口。

### 测试与判据

- 体积8m³、六面各4m²、原边界/质量/密度/h/V/E逐值不变：1e-12几何量/精确数组。
- 只移动中心点；解析mu、质心、完整C矩阵与直接组装比较：1e-12绝对（各自SI单位）。
- 原规则对照C<=1e-10、仿射/矩阵SBP/常应力节点残差<=1e-10。
- 扰动原corrected梯度仿射误差<=1e-10，实际cond(M)<1e6；
  不把条件良好误认为体积—边界相容。
- 使用mock哨兵保证不相容时project确实NOT_RUN，不能只写一条warning继续。
- 新回归的PASS只表示独立解析反例与拒绝行为被确认。扰动原型本身仍FAIL/BLOCKED。
- 原STEP2十项及旧boundary-red重跑；旧门槛不变。
- 扰动分支新算子谱/零模态、实际动态畸变、空间/时间敏感性、耦合、AFSD均NOT_RUN。

每脚本subprocess timeout=60s，BLAS/OMP单线程沿用上轮，仅用于确定性小矩阵。
落盘固定配置hash、源hash、原位置/新位置/V/E/C；不包含任何物理有效末态。
probe命令若矩不相容返回exit1，日志与报告标BLOCKED；不吞成“全部通过”。

## 实测结论：直接推广被否定，不改权重继续追绿

探针exit1，状态BLOCKED_BY_MOMENT_COMPATIBILITY；诊断正常完成。
控制点云PASS，扰动点云没有被接受；独立checked_project在投影前明确拒绝。
无新物理末态、无REJECTED物理candidate，只有标为REJECTED的静态采样夹具。
求解器源码、配置、STEP2算法、旧测试与失败全部不变。

| 固定125点/同一物理立方体 | 原规则对照 | 中心点扰动 |
|---|---:|---:|
| C元素绝对值最大值，m² | 0 | 0.006250000000000089 |
| 原M最大条件数 | 3.2550731723 | 3.2550731723 |
| 原corrected仿射误差 | 1.28608e-15 | 1.22125e-15 |
| NNPS边数 | 2557 | 2557 |
| 本轮投影构造 | RUN，节点patch 6.69342e-16 | NOT_RUN，必要条件不成立 |

解析预期0.00625与实测相符。完整三个4x4矩阵见report.json/NPZ及独立test的
字面常数参考；不是被测函数生成reference。控制算子SBP误差0、仿射误差9.72e-15。

CONFIRMED：这组固定H、几何节点对角E及扰动P不能同时满足严格分部积分和仿射再现，
即使扩大邻域、改核或提高线性求解精度也不能解除这个必要条件矛盾。
不据此否定所有不规则SPH/边界重构方案，不把静态采样变化当作真实物理运动。

RULED_OUT（仅本反例）：此次拒绝不是邻域病态、质量变化、几何边界移动、
浮点舍入或SVD迭代不够；求解器都尚未进入投影构造。
相容性缺口约为容限的6.25e7倍；没有用伪逆/最小二乘残差掩盖它。
扰动分支零模态检查因没有合格算子而NOT_RUN，不声称其无零模态。

实跑清单：

| 命令入口 | 状态 | 日志 |
|---|---|---|
| probe_boundary_perturbation.py | BLOCKED，exit1；0.146434s | probe.log |
| test_boundary_closure_step3.py | 6PASS，exit0；0.150s | diagnostic_tests.log |
| test_boundary_closure_step3.py（新增文件回读回归后） | 最终7PASS，exit0；0.342s | diagnostic_tests_final.log |
| test_boundary_closure_step2.py | 10PASS，exit0；0.358s | step2_regression.log |
| test_spatial_probes.py --boundary-red | 原FAIL，exit1；残差0.3261571877 | original_boundary_RED.log |
| 完整回归/新材料测试/耦合/真实AFSD/时空敏感性/生产 | NOT_RUN | 不存在新通过声明 |

本轮未反复修测试输入或调参数。每入口实际耗时不足1s，60s超时包装。
并发运行的只是独立小测试命令，不是求解器并行优化。
原559项冻结在修改前PASS，结束的新冻结只允许HANDOFF变化。

收尾文件回读发现新增脚本的配置SHA256不等于落盘字节hash：Windows默认将LF
转换为CRLF。artifact_readback.log保留FAIL，artifact_hash_diagnosis.log证实
LF规范化后的hash与报告一致、21项源码hash当时均一致，数值内容没有不同。
只将新脚本fixture_config.json写入指定newline=LF，新增实际main写入/回读回归。
修前脚本及6测试源码快照保留，初始输出不覆盖。最终7PASS；最终探针仍exit1/
BLOCKED，0.141005s，artifact_readback_final.log PASS。
这项记录修复不是边界/物理修复，矩阵、几何、材料、权重和容限均未变。

### 下一轮范围与路线判据

不再扫描更多扰动、增加投影迭代或直接改生产mechanics.py。
先给一个质量一致的“物理体积—数值求积—边界迹”接口设计：
物理V_i=m_i/rho_i、质量/密度更新不能为了patch被偷偷改写；
若引入其他求积权重或重构边界迹，必须分别命名、说明功与边界牵引含义，
并以独立几何/制造解验证，不能从-F_internal反推牵引。
该设计与可行性小基准是下一工作单元，本轮未实施。

如果该约束下没有可推广方案，应明确提交方法取舍，而不是继续逐粒子修复。
材料Cauchy/Jaumann闭循环非零功0.146111J仍是另一独立方法门禁。
完整复现时间判断另见COMPLETION_OUTLOOK_20261001.md；不因此扩大本轮范围。

## 文件和重跑

新增scripts/probe_boundary_perturbation.py的moment_defect/checked_project/run_case/evaluate，
新增tests/review/test_boundary_closure_step3.py的7项回归。仅独立脚本导入，不接生产。
outputs/boundary_step3_20261001_final/含配置及其SHA256、源hash、JSON、
control_sampling_fixture.npz及REJECTED_sampling_fixture.npz。
REJECTED文件没有candidate_D_per_m，更不包含有效物理末态。
初始outputs/boundary_step3_20261001/保留已知换行hash不一致，不用作最终hash证据。
logs/boundary_step3_20261001/含预注册、命令日志、原HANDOFF、run_ledger.json。
新BOUNDARY_STEP3_FREEZE.json验证文件完整性，不是方法/论文验收。

项目根，沿用OPENBLAS_NUM_THREADS/MKL_NUM_THREADS/OMP_NUM_THREADS=1。
以下是60s subprocess包装实际执行的内部命令；再次运行使用新输出目录：

    .envs/pysph/python.exe -B scripts/probe_boundary_perturbation.py --output outputs/boundary_step3_recheck_01
    .envs/pysph/python.exe -B tests/review/test_boundary_closure_step3.py
    .envs/pysph/python.exe -B tests/review/test_boundary_closure_step2.py
    .envs/pysph/python.exe -B tests/review/test_spatial_probes.py --boundary-red

预期exit分别1、0、0、1；两个exit1均有明确方法状态，不应改成“测试全部通过”。


