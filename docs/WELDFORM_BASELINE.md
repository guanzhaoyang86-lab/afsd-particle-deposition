# WeldForm CPU 基线审计

审计日期：2026-09-24

## 结论

后端候选固定为 WeldForm CPU，commit `19b14e20fb0940b2991cb141ef9e1b80cb796e0b`。最终状态为 **`BLOCKED_BUILD`**：第 2/2 轮配置成功、`LSDynaReader` 构建到 60% 并完成，但 `WeldForm.cpp` 编译因缺少 `gsl/gsl_math.h` 退出 2；没有生成 `WeldForm.exe`。因此固体、热学、接触三类输入均只是上游源码中可找到，全部为 `NOT_RUN_BUILD_BLOCKED`，没有任何求解结果，也没有复现论文。

本报告严格区分：

- **已有代码**：只能证明固定 commit 中存在相关实现和调用语句。
- **实际验证通过**：仅限命令确实返回成功且有日志；本次没有任何 WeldForm 求解器运行通过。
- **仍需开发/验证**：进入 AFSD 基线前必须完成的适配、修复和数值测试。

权威逐字命令与退出码见 `logs/upstream/weldform/00_exact_commands_and_status.log`，原始日志见同目录 `01`–`15`。上游固定和许可证说明见 `third_party/WeldForm.UPSTREAM.md`。

## 上游、许可证与依赖固定

| 项目 | 固定值 | 实测状态 |
|---|---|---|
| WeldForm | `https://github.com/luchete80/WeldForm.git` @ `19b14e20fb0940b2991cb141ef9e1b80cb796e0b` | detached HEAD；工作树无受控文件修改 |
| WeldForm 许可证 | `third_party/WeldForm/LICENSE.txt`，GNU GPL v3 或更高版本 | SHA-256 `5C2212175AAFF4054EB204302B14085C6EEF8FD0F713770B4DCE0DF7B7DB3271` |
| LSDynaReader 子模块 | `10eda875f2324c1953539efe32a60368c238a1a9` | 用等价 HTTPS URL 初始化；最终 `git submodule status` 前缀为空格 |
| BLITZ fork | `luchete80/blitz-0.9-cmake` @ `da019d3911d99c4a9848e878d09e707429222e9c` | 源码及 `LICENSE`/`COPYING` 保留在 `.deps/weldform/src` |
| HDF5 1.12.0 | `HDFGroup/hdf5` @ `4aa40365623cf35643b262781f92578220c42bbb` | 静态 C/HL 库隔离安装成功；许可证保留 |

最初的 shallow clone 包装器在捕获 stdout 和 `$LASTEXITCODE` 前失败，所以 clone 的数值退出码诚实记录为 `UNKNOWN`；随后独立验证了 remote、commit、干净工作树、许可证和子模块，不能倒推或伪造 clone exit code。

依赖没有安装到系统环境。`SPH` 指向 `third_party/WeldForm`，`PKG` 指向 `.deps/weldform`。WeldForm README 要求“BLITZ++ 0.9”，但其链接 fork 的 CMake 元数据报告 1.0.3，且缺少 `blitz/config.cmake.h.in`；老 Autotools 的 `config.guess` 又不识别当前 MINGW64。最终隔离 BLITZ 安装使用固定源码中的预生成头，并实际编译 `src/globals.cpp` 得到 `libblitz.a`，不能把它表述成标准上游包构建通过。

## 两轮构建证据

本次严格停止在两轮 WeldForm 构建，没有为补 GSL 发起第 3 轮。

| 轮次 | 实际命令摘要 | 退出码 | 结果 |
|---|---|---:|---|
| 1/2 | `cmake -S third_party/WeldForm -B build/weldform-attempt1 -G "MinGW Makefiles"` | 1 | 上游 `cmake_minimum_required(VERSION 2.8)` 与 CMake 4.1 不兼容，配置失败 |
| 2/2 配置 | 同上，使用隔离依赖并加 `-DCMAKE_POLICY_VERSION_MINIMUM=3.5 -DCMAKE_BUILD_TYPE=Release` | 0 | GNU C/C++ 8.1.0 配置成功；配置阶段已警告 GSL/GSL CBLAS 缺失 |
| 2/2 构建 | `cmake --build build/weldform-attempt2 --target WeldForm -j 8` | 2 | `lsdynareader` 构建完成；`External/matvec.h:34` 的 `#include <gsl/gsl_math.h>` 致命失败 |

最终检查未找到任何 `WeldForm*.exe`。`cmake --build ... --target help` 只盘点目标，没有触发第 3 次构建；目标中有 `WeldForm`、`lsdynareader` 和 `test`，但 `ctest --test-dir build/weldform-attempt2 -N` 报告 `Total Tests: 0`。

## 三类上游案例

以下命令是构建成功后应执行的真实 CLI 形式；由于可执行文件不存在，本阶段**没有执行**这些命令。

| 类别 | 输入 | 预定命令 | 状态 | 静态风险 |
|---|---|---|---|---|
| 最小固体 | `examples/input/Compression.json` | `build/weldform-attempt2/WeldForm.exe third_party/WeldForm/examples/input/Compression.json` | `AVAILABLE_SOURCE_INPUT_ONLY; NOT_RUN_BUILD_BLOCKED` | 使用当前解析器需要的 `const`，选择 `Mech-LeapFrog`；仍需真实 smoke run |
| 热学 | `examples/input/SolidBlock_Heat_Conv.json` | `build/weldform-attempt2/WeldForm.exe third_party/WeldForm/examples/input/SolidBlock_Heat_Conv.json` | `AVAILABLE_SOURCE_INPUT_ONLY; NOT_RUN_BUILD_BLOCKED` | 选择 `Thermal`，含导热参数和对流边界；未验证温升或能量守恒 |
| 刚体接触 | `examples/input/Contact_Compression.json` | `build/weldform-attempt2/WeldForm.exe third_party/WeldForm/examples/input/Contact_Compression.json` | `AVAILABLE_SOURCE_INPUT_ONLY; NOT_RUN_BUILD_BLOCKED` | 使用旧 `const1/const2` 字段且选择有分派缺陷的 `Mech-Fraser`；不能作为现成通过案例 |

三类输入都存在，因此没有一类标为 `NOT_AVAILABLE`。不过 `examples/CMakeLists.txt:6-123` 的案例可执行列表全部被注释，实际默认构建入口只有根目标 `WeldForm` 加 JSON。早先发现的 `Compression_fast_short.json` 更短，但仍使用当前解析器不读取的 `const1/const2`，所以最终固体 smoke 输入改选 `Compression.json`。

## 实际调用路径

活动入口不是 README 或绘图，而是下面的源码链：

```text
Source/WeldForm.cpp
  -> JSON 解析、Material_ 派生类、3D 几何和 TriMesh 刚体
  -> solver 字符串分派
  -> Domain::SolveDiffUpdateLeapFrog (Source/SolverLeapfrog.cpp)
       -> MainNeighbourSearch (Source/Neighbour.cpp)
       -> CalcAccel / 接触力
       -> 速度、密度和位置更新
       -> CalcRateTensors (Source/InteractionAlt.cpp)
       -> Particle::CalcStressStrain (Source/Particle.cpp)
       -> CalcPlasticWorkHeat -> ThermalCalcs (Source/Thermal.cpp)
       -> WriteXDMF (Source/Output.cpp)
```

包含关系同样是活动证据：`WeldForm.cpp:21` 包含 `Domain.h`，`WeldForm.cpp:25,31` 包含 `InteractionAlt.cpp` 和 `SolverLeapfrog.cpp`；`Domain.h:544-552` 直接包含 `Domain.cpp`、`Neighbour.cpp`、`Output.cpp`、`Thermal.cpp`。这说明相关代码进入同一翻译单元，但由于 GSL 头在早期中断，仍不等于这些路径编译或运行验证通过。

## 能力审查：已有 / 实测 / 待开发

| 能力 | 上游已有代码（固定 commit 的静态证据） | 本阶段实际验证 | AFSD 仍需开发或验证 |
|---|---|---|---|
| 真三维邻域 | `WeldForm.cpp:175,420-449` 默认 3D，仅薄于 `h` 时降为 2D；`Neighbour.cpp:56-89,109-115` 分配 xyz cell/HOC；`346-354,418-472` 遍历 z 层和上方 9 个 cells | **未运行**。不能用 3D 几何或 VTK 证明求解器 3D | 建立非共面粒子邻域金标准，断言 z 邻居、成对唯一性和多线程一致性；案例必须断言 `Dimension==3` |
| 真三维率张量/应力 | `InteractionAlt.cpp:263-361` 计算 xz、yz、zz 应变率和三轴旋转率；`Particle.h:112-135` 保存完整 `Mat3_t`；`Particle.cpp:960-1059` 有 3D Jaumann 项、J2 等效量和塑性张量更新 | **未运行、未编译到目标完成** | 做 3D patch test、刚体旋转客观率测试、单轴/纯剪切解析回归；审查 yield 更新时序和单位 |
| 塑性材料 | `Material.h:96-101` 基类默认返回 0；`117-224` 声明 Johnson–Cook、GMT、Hollomon；`Material.cpp:6-161` 有相应流动应力/切线；Bilinear 由 `Ep` 路径处理；`CalcStressStrain` 做 J2 限制，`CalcPlasticWorkHeat` 在 `Particle.cpp:922-939` | **未运行**。只能确认源码存在 | 论文所需 FKS **不存在**；须新增清晰材料接口和 FKS 后再做。先对现有 Bilinear/JC/GMT 做曲线、屈服返回、温度/应变率边界和塑性功测试 |
| 热传导/热源 | `Domain.h:208-211,249-254` 声明热求解；`Thermal.cpp:53-219` 汇总粒子导热、对流、源项、塑性热和接触项；`567-718` 是独立热循环；`721-748` 被力学 LeapFrog 每步调用 | 热输入存在，**未运行** | `WeldForm.cpp:325-329` 的 `k_T/cp_T` 未初始化后即尝试读取，缺项不能默认为零；配置适配器必须按启用物理显式校验。`Thermal.cpp:60-138` 多线程对共享 `temp/axiterm` 做无显式同步累加，需确定性/TSAN 检验并改成 reduction；做解析导热、对流和全局能量守恒测试 |
| 刚体接触 | `WeldForm.cpp:581-734` 创建 Plane/File `TriMesh` 并加接触粒子；`Domain.cpp:169` 默认 Wang；`SolverLeapfrog.cpp:234-239` 调 Wang/Fraser/Seo/LS；`Contact_Wang.cpp:15-345` 计算法向、摩擦、热生成和接触导热 | 接触输入存在，**未运行** | `WeldForm.cpp:695` 的算法读取被注释，`Fraser` 选项明确回退 Wang；必须让配置实际控制算法。`Contact_Wang.cpp:323-328` 只更新固体侧 `q_cont_conv`，刚体侧 `tot_cont_heat_cond` 累加被注释，需做双侧能量平衡。增加穿透、反力、摩擦功和刚体热交换测试 |
| 时间积分器 | 活动入口包括 `SolverLeapfrog.cpp`、`SolverFraser.cpp`、`SolverKickDrift.cpp`；LeapFrog 在 `124-335` 完成邻域、力、接触、速度/密度/位置、应力和热更新 | **未运行** | `WeldForm.cpp:217` 读取 `integrationMethod` 到 `dom.Scheme`，但 `929-940` 实际按 `solver` 字符串选积分器。两个连续 `if` 使 `Mech-Fraser` 返回后落入第二链的 `Invalid solver`。统一成单一分派表并做一阶/二阶收敛与能量漂移测试。`SolverVerlet.cpp` 虽存在且 `Domain.h` 声明，但活动入口未包含/分派 |
| 输出与粒子身份 | `Output.cpp:75-270` 写 HDF5/XDMF 的位置、速度、应力、温度等；`Particle.h:67-68` 有 `ID/ID_orig` | **未运行**，没有输出文件 | `ID` 注释和构造器表明它是区域 tag，不是唯一粒子 ID；`ID_orig` 仅在接触求解器 `SolverLeapfrog.cpp:79-83` 中复制 tag。需新增不可变 `particle_id` 和初始坐标。另有明确错误：`Output.cpp:250-251` 创建 `ShearS` 数据集却传入 `Sigma` 数组。输出必须通过项目 schema/VTK round-trip 测试 |
| 上游测试 | 配置生成了 `test` 目标 | `ctest -N` 实测为 0 tests | 所有物理能力都必须由本项目添加小型、可判定的测试，不能以案例能启动或图像“看起来对”代替 |

## Material.h、Domain.h、Thermal.cpp 的重点判断

### `Source/Material.h` / `Material.cpp`

材料多态接口存在，但基类多个虚函数直接返回 `0.0`，`getRefStrainRate()` 的基类实现甚至没有返回值。当前 JSON 入口只接受 `Bilinear`、`Hollomon`、`JohnsonCook`、`GMT`（`WeldForm.cpp:244-290`），没有论文所需 FKS。Johnson–Cook 仅在应变率高于参考值时启用对数倍率，温度归一化没有在函数内夹紧；GMT 才按输入范围夹紧应变、应变率和温度。它们都属于“代码存在”，不是数值正确性证明。

当前解析器只读 `Materials[0].const`（`WeldForm.cpp:252`）。若旧输入只提供 `const1/const2`，长度 10 且预置为零的 `c` 不会得到预期参数；这就是不选择 `Compression_fast_short.json` 做最终 smoke 的原因。AFSD 配置层应执行带单位、来源和 `resolved/pending` 状态的显式转换，缺失值报错，不向上游静默填零。

### `Source/Domain.h` / 实现包含链

`Domain.h` 同时承载邻域、接触、热、积分器和输出接口，且直接包含多个 `.cpp`，耦合度很高。它声明的能力比活动主程序多，例如 Verlet；审计因此以 `WeldForm.cpp` 的真实分派和 `SolverLeapfrog.cpp` 的真实调用为准，不能把声明当成完成。`Domain.cpp:78-175` 默认 `Dimension=2`、`Scheme=0`、热/接触关闭、接触算法 Wang；主程序随后覆盖部分状态。下一阶段适配器必须在生成输入和读取启动日志时双重检查最终状态。

### `Source/Thermal.cpp`

独立热求解确实调用邻域搜索、初始 `CalcTempIncPP()` 和逐步 `ThermalCalcs()`，耦合力学路径也在塑性功之后调用 `ThermalCalcs()`。但 `Thermal.cpp:3` 直接包含名为 `ThermalTest.cpp` 的实现，且活跃导热累加存在并行共享写风险；接触导热又缺少刚体侧累加。这些问题必须先经解析解、线程确定性和能量账本验证，不能因为 README 列出 “Thermal Solver” 就标为已验证。

## 仅供未来重新评估 WeldForm 时的修复清单

本项目最终已选择 PySPH 作为唯一下一阶段后端；当前权威计划是
`docs/IMPLEMENTATION_PLAN.md`。以下内容仅记录若未来明确决定重新开启
这个已阻塞候选时所需的工作，不是当前任务，也不得与 PySPH 并行开发。

优先保持单后端，不并行开发 GPU/PySPH。建议顺序如下：

1. **解除唯一构建阻塞**：在 `.deps/weldform` 隔离固定 GSL 2.5，保存许可证和 commit/tag；只进行下一阶段第 1 次受控构建。若 MinGW GSL 仍阻塞，记录新的完整编译命令和首个致命错误，不无限修环境。
2. **建立后端适配器**：新增 `src/afsd_scaffold/backends/weldform.py`，只负责从三工况 schema 生成单位一致的 WeldForm JSON、执行单一二进制、捕获 manifest/log；所有 `pending` 参数拒绝启动，显示单位转换不得进入求解输入。
3. **清理入口和校验**：修改或维护小型 patch queue 覆盖 `Source/WeldForm.cpp` 的必填参数检查、`const` schema、接触算法读取和互斥 solver 分派；不要在此阶段实现 FKS。
4. **最小 smoke**：先运行 `Compression.json` 并要求 exit 0、有限值、至少两个时间帧；然后分别运行热学和接触输入。原始 stdout/stderr、输入快照、二进制/上游 commit 必须进入 run manifest。
5. **源码级物理测试**：新增 `tests/test_weldform_3d_patch.py`、`test_weldform_material_curves.py`、`test_weldform_thermal_balance.py`、`test_weldform_contact_balance.py`、`test_weldform_integrator.py`。每项先覆盖现有上游行为，再决定最小修复。
6. **输出适配**：修复 `Source/Output.cpp` 的 `ShearS` 指针错误；在 `Particle.h/.cpp` 增加不可变唯一 ID 和初始位置，或在受验证的适配层稳定映射。把 HDF5/XDMF 转为 `schemas/particle-record.schema.json` 所需字段，并用现有 `tests/test_output_contract.py` 与独立 VTK fixture 做 round-trip。
7. **热/接触并行安全**：重构 `Source/Thermal.cpp` 的 pair accumulation 为线程私有数组加 reduction；补全 `Source/Contact_Wang.cpp` 的双侧接触热账本，再做 1/2/多线程等价性测试。
8. **后续物理开发门槛**：只有上述基线测试通过后，才在 `Source/Material.h/.cpp` 设计 FKS；论文参数仍为 pending 时不得以零或任意常数代替。

## 阻塞与选择门槛

- 当前硬阻塞：隔离 GSL 2.5 缺失，导致 CPU 主程序无法编译。
- 即使构建解除，接触样例分派、旧 JSON 字段、热并行累加、单侧接触热、非唯一 ID 和输出数组错误仍阻止其成为 AFSD 可验证基线。
- WeldFormGPU 只有在硬件/CUDA 工具链与项目门槛都满足时才值得评估；它不能用来掩盖 CPU 基线的接口与验证缺口。
- PySPH 只在 CPU WeldForm 经限定修复后仍有实测阻塞时作为后备候选，不能与 WeldForm 同时开发。
