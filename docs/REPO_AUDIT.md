# REPO_AUDIT — 实测后端审计与选型

## 2026-09-24 后续实施增量（当前状态）

下文仓库比较/构建结论保留为准备阶段历史，没有重新选择后端。
PySPH固定commit不变、上游tracked源码干净。现已新增项目host适配层的真实三维
hyperbolic核、完整J2/FKS/Jaumann、bulk压力闭合、导热/功热、有限面刚体接触、
执行器、持久ID和VTK/NPZ输出；实际复用上游粒子、3D邻域、leapfrog、壁接触和检查点。

已有代码/实测通过/仍需开发分别见PHYSICS_STATUS.md四层表、EQUATION_CODE_TEST_MAP.md
及IMPLEMENTATION_PLAN.md顶部。已实测19骨架+21后端测试；三档657粒子微秒DEBUG耦合
和重启/告警/回读通过。未实现论文完整三角面粘滑，未解析P40/P41参数，尚无沉积尺度
或实验验证，不能把本页历史的“缺FKS”读作当前未编码，也不能把新增代码读作论文已复现。

实际证据：logs/physics/20260924_physics_verified/commands_and_status.json；
outputs/debug/20260924_physics_verified/。当前CPU主机求值路径无JIT/无GPU。

## 以下为准备阶段审计记录

审计日期：2026-09-24。范围是本工作目录内固定的上游源码、当前执行环境和本次真实命令；不是论文复现认证。README、文件名和三维绘图都不计作能力验证。

## 1. 结论

下一阶段只选 **PySPH** 作为开发基线：`pypr/pysph` 固定在 commit `1cd3dea9cfaac7d9f1bda2a2491ce4b42a31dbe5`，BSD 3-Clause。项目内 Conda 环境一次安装成功，官方通用例、二维固体例和真正 `dim=3` 的固体例均完成有界 smoke run，输出也被重新读取验证。

这只是“有一个能运行的开源基线”。PySPH 现有的 FKS、本项目所需完整三维塑性、金属温度/导热、变形体—刚体三角面接触和论文 leapfrog 路径均未达到 AFSD 要求。不得称为论文复现。

WeldForm CPU 按两轮修复限额停在缺少 GSL 的构建错误；WeldFormGPU 编译通过，但真实固体例在第一次 GPU 邻域搜索处失败。两者均保留固定源码、许可证、日志和源码审计，仅作参考，不与 PySPH 并行开发。

## 2. 当前执行环境门槛

`ENVIRONMENT.md` 记录的是本次 Codex 进程实际可见的原生 Windows 环境，不外推为其他电脑或容器：Windows 11 build 26200，Intel Core Ultra 9 275HX（24C/24T），约 64 GiB RAM，RTX 5090 Laptop GPU（24,463 MiB，CC 12.0），driver 595.97，`nvcc` 12.8，MSVC 19.44、MinGW g++ 8.1、CMake 4.1.0-rc1、默认 Python 3.13.7。GPU 硬件门槛通过只授权评估 WeldFormGPU，不代表上游可运行。

PySPH 依赖安装在 `.envs/pysph`，没有升级系统 Python 或系统包。固定环境清单见 `logs/upstream/pysph/12_conda_list_explicit.stdout.log` 和 `13_pip_freeze.stdout.log`。

## 3. 候选的构建/运行结果

| 候选 | 固定版本 / 许可证 | 构建或安装 | 实际案例 | 最终状态 |
|---|---|---|---|---|
| WeldForm CPU | `19b14e20fb0940b2991cb141ef9e1b80cb796e0b` / GPL-3.0-or-later | 1/2：CMake 4 policy 失败；2/2：配置成功，构建因 `External/matvec.h:34` 缺 `gsl/gsl_math.h`，exit 2 | 固体、热学、接触输入均存在，但无可执行文件，全部 `NOT_RUN_BUILD_BLOCKED` | `BLOCKED_BUILD` |
| WeldFormGPU | `ed567833a7490bf5df829817f5e5214929045a9e` / GPL-3.0 | 1/2：CMake policy 失败；2/2：VS2022 + CUDA 12.8 编译 exit 0，二进制含 `sm_75` 与 `sm_120` | 原样 `Compression.json` 创建 11,316 粒子后，首次 nsearch 报 CUB `invalid argument` / `invalid device ordinal`；无时间步和输出，236.8 s 后终止 | `BUILD_PASS / BLOCKED_RUNTIME_NEIGHBOR_SEARCH` |
| PySPH | `1cd3dea9cfaac7d9f1bda2a2491ce4b42a31dbe5` / BSD 3-Clause | 项目内 Python 3.11 环境、requirements 和 `setup.py develop` 均首轮 exit 0 | 通用 2 步、固体 2D 2 步、固体 3D 1 步均 exit 0；NPZ 可读 | **`SELECTED_WORKING_BASELINE`** |

逐字命令、退出码和原始输出分别在：

- `logs/upstream/weldform/00_exact_commands_and_status.log`
- `logs/upstream/weldformgpu/00_exact_commands_and_status.log`
- `logs/upstream/pysph/00_exact_commands_and_status.log`

详细报告见 `docs/WELDFORM_BASELINE.md`、`docs/WELDFORMGPU_EVALUATION.md` 和 `docs/PYSPH_BASELINE.md`。来源、commit、工作树状态和许可证哈希见 `third_party/*.UPSTREAM.md`；许可证原文保留在各上游仓库中。

## 4. WeldForm CPU 重点源码审查

活动调用链是：

```text
Source/WeldForm.cpp
  -> Domain::SolveDiffUpdateLeapFrog (Source/SolverLeapfrog.cpp)
  -> MainNeighbourSearch (Source/Neighbour.cpp)
  -> CalcRateTensors (Source/InteractionAlt.cpp)
  -> Particle::CalcStressStrain (Source/Particle.cpp)
  -> CalcPlasticWorkHeat / ThermalCalcs (Source/Thermal.cpp)
  -> WriteXDMF (Source/Output.cpp)
```

`Source/Material.h/.cpp` 有 Johnson–Cook、GMT、Hollomon和基类，但没有 FKS；基类多个方法返回零，`getRefStrainRate()` 还有无返回路径。`Source/Domain.h` 声明的能力多于主程序实际分派，且直接包含多个 `.cpp`，因此以 `WeldForm.cpp` 与 `SolverLeapfrog.cpp` 的实际调用为准。`Source/Thermal.cpp` 的独立热循环和力学耦合入口都存在，但活动导热累加有共享写风险，接触热只更新固体侧；这都没有运行通过。

静态源码确认了 z 向 cell 遍历、完整 3×3 率/应力字段、J2 更新、Wang 接触和 Leapfrog 代码存在；实测只到编译早期 GSL 阻塞。另发现：`integrationMethod` 不控制最终分派，`Mech-Fraser` 分派链可能继续报 `Invalid solver`，旧案例 `const1/const2` 不被当前材料解析器读取，`Output.cpp` 的 `ShearS` 数据集误传 `Sigma`，`ID/ID_orig` 不是持久唯一粒子 ID。详见 CPU 报告，不能把这些“已有代码”写成“验证通过”。

## 5. WeldFormGPU 能力边界

源码存在 27-cell 三维邻域、六分量应力/Jaumann 更新、J2/多种材料分支、Fraser 积分和接触 kernel；但实际运行在第一次 `find_neighbors()` 即失败，没有到达一个应力或积分时间步。热目标与热参数解析被注释，故标 `NOT_AVAILABLE_AS_ACTIVE_ENTRYPOINT`；接触源码和 JSON 存在，但共享邻域前置阻塞，标 `FOUND / NOT_RUN_BLOCKED`。材料模型枚举和部分 `Ep` 路径还存在未初始化风险。

因此 GPU 编译成功不是 GPU 求解成功，更不是三维物理验证。

## 6. PySPH 实测基线

项目隔离环境的核心版本：Python 3.11.16、PySPH `1.0b3.dev0`（上述 commit）、NumPy 2.4.6、Cython 3.3.0、Compyle 0.9.1。主要真实命令为：

```text
conda create --prefix .envs\pysph python=3.11 pip setuptools wheel -y
conda run --prefix .envs\pysph python -m pip install -r third_party\PySPH\requirements.txt
conda run --prefix ..\..\.envs\pysph python setup.py develop
conda run --prefix .envs\pysph python third_party\PySPH\pysph\examples\elliptical_drop_simple.py --max-steps 2 --pfreq 1 --directory outputs\upstream\pysph\elliptical_drop_smoke
conda run --prefix ..\..\..\..\.envs\pysph python ..\..\..\..\third_party\PySPH\pysph\examples\solid_mech\rings.py --max-steps 2 --pfreq 1
conda run --prefix ..\..\..\..\.envs\pysph python ..\..\..\..\third_party\PySPH\pysph\examples\solid_mech\impact3d.py --max-steps 1 --pfreq 1
```

结果：

- 通用官方例：5,025 粒子、2 步、exit 0。
- `rings.py`：17,592 粒子、二维固体、2 步、exit 0。
- `impact3d.py`：32,805 个 plate 粒子 + 4,157 个 projectile 粒子，`Gaussian(dim=3)`、`Solver(dim=3)`、`VelocityGradient3D`、`HookesDeviatoricStressRate`、`EPECIntegrator(SolidMechStep)` 的实际路径完成 1 步，exit 0。
- `validate_outputs.py` 重新读取三个末态 NPZ，核对时间、数组名和粒子数，exit 0。

第一次通用例命令使用不支持的 `--output-dir`，CLI 在求解前 exit 2；依据上游参数定义改为 `--directory` 后通过。失败日志保留，且这不是第二轮依赖/构建修复。`pytest` 不在 requirements 环境中，`python -m pytest --version` exit 1；没有为此追加安装。因此上游单独的 3D NNPS 单测本次 `NOT_RUN_MISSING_TEST_DEPENDENCY`，但 `impact3d` 已真实经过三维 NNPS/方程执行。PySPH 默认把运行时生成的 Cython 模块缓存到 `<LOCAL_USER_PROFILE>\.pysph`；依赖本身仍在项目环境，后续运行器要显式记录/隔离该缓存策略。

## 7. 选中后端的能力矩阵

| 能力 | 固定 PySPH 中已有代码 | 本阶段实际验证通过 | AFSD 仍需开发/验证 |
|---|---|---|---|
| 真三维邻域 | `pysph/base/*nnps*.pyx` 支持 `dim=3`；上游测试含三维坐标与 brute-force 对照 | `impact3d.py` 36,962 粒子完成 1 步；不是绘图证据 | 增加本项目非共面邻域金标准、排序后 ID 同步和重复性测试；上游 NNPS pytest 本次未跑 |
| 三维率/应力 | `pysph/sph/solid_mech/basic.py` 的 `HookesDeviatoricStressRate` 使用九个速度梯度和六个应力分量 | `impact3d.py` 的 `VelocityGradient3D`、动量、能量和应力率活动路径完成 1 步 | 做仿射 patch、单轴/三向剪切、预应力刚体旋转客观性和能量测试；一时间步不证明数值正确 |
| 塑性 | `pysph/sph/solid_mech/hvi.py` 有 `VonMisesPlasticity2D`；3D 示例也调用它 | 仅证明该调用没有在首步崩溃 | 该类的屈服不变量只用 `s00/s01/s11`，不能算完整 3D 塑性；FKS 不存在。需新建完整 3D 返回映射、FKS 状态和材料点测试 |
| 热传导 | 固体能量方程会更新比内能；gas-dynamics 有人工热传导相关代码 | `impact3d` 能量方程首步运行 | 金属温度、k(T)、Cp(T)、塑性热、实体间导热、摩擦/接触热的活动固体接口 **NOT_AVAILABLE**；需独立实现并做守恒测试 |
| 刚体接触 | `rigid_body.py` 有刚体/DEM 球接触；固体 impact 通过普通 SPH 粒子相互作用 | 没有 AFSD 型刚体接触实测 | 变形固体—旋转/平移三角网格工具/底座接触、摩擦和热分配 **NOT_AVAILABLE**；不能把 DEM 球接触或两种可变形粒子相互作用当替代 |
| 积分器 | `PECIntegrator`、`EPECIntegrator`、`SolidMechStep` 存在 | rings 的 PEC 和 impact3d 的 EPEC 有界运行通过 | 论文要求的 forward-difference leapfrog 路径未验证/未提供；需定义状态时序、稳定步长、收敛和能量漂移测试 |
| 输出/身份 | 上游写 NPZ/info/log/profile，ParticleArray 可添加自定义属性 | 三组 NPZ 实际生成并可重新读取 | 上游输出尚不符合本项目 schema；新增不可变 `particle_id`、初始位置/材料/标签和原子 manifest/VTK 转换，禁止以数组下标充当 ID |

## 8. 本项目已建立并实测的接口骨架

- `configs/base.json` 与 `configs/cases/{starved,copper-core,overfed}.json`：计算值统一 SI；每个参数有 `value/unit/source/status`；未决值为 `null`，不填零。三个 feed 速度分别为 0.00106、0.00212、0.00424 m/s。
- `schemas/afsd-case.schema.json`、`run-manifest.schema.json`、`particle-record.schema.json` 与 `src/afsd_scaffold/`：配置合并/拒绝未决项、输出契约、legacy VTK 读写。
- `tests/run_tests.py`：19 项通过，覆盖三工况、非有限/非 SI 拒绝、运行门禁、manifest 来源要求、ID/标签和 VTK round-trip；最终日志为 `logs/tests/20260924_final_scaffold.log`。
- `outputs/geometry_fixture/AFSD_GEOMETRY_FIXTURE.vtk`：118 点静态三维 fixture，manifest 明确 `GEOMETRY_FIXTURE`、`physics_claims_allowed=false`。它只验证几何标签与 I/O，不是仿真结果。

## 9. 未决阻塞和选型约束

论文侧阻塞仍包括：FKS 参数/参考量、k(T)/Cp(T)、完整工具/基板/仿真进料几何、启动/路径程序、摩擦/热接触/热分配、EOS、h/spacing、稳定化/CFL/缩放、示踪采样/坐标/测量算子和 Fig.8 原始数据。详见 `docs/PAPER_SPEC.md` 与 `docs/ASSUMPTIONS.md`。任何缺失项都不得默认为零。

下一阶段只在 PySPH 适配层开发并先做小型模块测试。WeldForm/GPU 的源码可以作为公式/接触实现参考，但不得形成第二套活动后端。只有完整 3D 塑性、热、刚体接触、积分和输出门槛分别通过后，才允许创建 coarse `DEBUG` AFSD case；参数门禁未清除前不允许 `PAPER_COMPARISON`。
