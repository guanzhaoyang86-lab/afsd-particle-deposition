# WeldFormGPU 候选评估

评估日期：2026-09-24（UTC-05:00）  
最终状态：**BUILD_PASS / BLOCKED_RUNTIME_NEIGHBOR_SEARCH**  
结论：本次没有任何 WeldFormGPU 物理案例完成一个时间步；不能把编译通过写成基线运行通过，也不能据此声称 AFSD 已复现。

## 范围与硬件门槛

仅因 `ENVIRONMENT.md` 已实测当前执行进程可见 RTX 5090 Laptop（24,463 MiB、compute capability 12.0）、driver 595.97、`nvcc` 12.8，才进入 GPU 候选评估。没有把其他用户电脑或假想容器硬件用于结论；没有系统级安装，也没有修改 WeldForm CPU 后端。

本候选严格限制为两轮有证据的构建修复。第二轮成功后只做相同配置的 clean rebuild 留全量日志；这不是第三轮修复。

## 上游固定与许可证

- 上游：`https://github.com/luchete80/WeldFormGPU.git`
- 固定 commit：`ed567833a7490bf5df829817f5e5214929045a9e`
- `LSDynaReader` 子模块：`c352a4d893c6f22b93c523c02c31f1a533ba04f2`
- 许可证：仓库内 `License.txt`，首行为 `GNU GENERAL PUBLIC LICENSE`，SHA-256 为 `E5CDBABBE1AE60E8EACD40D03061826291E2F65447B8466F1BE8DA8C512495A2`
- 最终 `git status --porcelain=v1 --untracked-files=all` 为空，上游工作树未改；构建产物均在 `build/weldformgpu-attempt*`。

上游 `.gitmodules` 使用 SSH URL。首次 `git clone --recurse-submodules` 主仓库已下载，但子模块因 SSH host-key 校验退出 1；随后只在单条命令中用 `url.https://github.com/.insteadOf=git@github.com:` 获取子模块，未修改 `.gitmodules`。详见 `logs/upstream/weldformgpu/01_clone_attempt.log`、`02_provenance.log` 和 `third_party/WeldFormGPU.UPSTREAM.md`。

## 构建实测

上游 README 只声称 VS2019 与 CUDA 11.4 已测试；根 `make.bat` 还硬编码 VS2019 路径。实际使用本机已存在的 VS2022 Community 17.14 / MSVC 19.44、CMake 4.1.0-rc1、NMake 和 CUDA 12.8。

### 第 1 轮：失败

```text
cmd.exe /d /s /c call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat" && cmake -S third_party\WeldFormGPU -B build\weldformgpu-attempt1 -G "NMake Makefiles"
```

退出码 1。CMake 在编译前拒绝上游 `cmake_minimum_required(VERSION 2.8.12)`，明确提示 CMake 4 已移除 `<3.5` compatibility。证据：`03_attempt1_configure.log`。

### 第 2 轮：配置与构建通过

```text
cmd.exe /d /s /c call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat" && cmake -S third_party\WeldFormGPU -B build\weldformgpu-attempt2 -G "NMake Makefiles" -DCMAKE_POLICY_VERSION_MINIMUM=3.5 -DCMAKE_CUDA_FLAGS="-gencode=arch=compute_120,code=sm_120" && cmake --build build\weldformgpu-attempt2 --target WeldFormGPU
```

退出码 0。policy 参数采用第 1 轮错误原文建议。`src/CMakeLists.txt:175` 将目标架构硬编码为 75；命令行另加 `sm_120`，所以没有改上游源码，又避免只生成不能在本机 CC 12.0 上运行的 `sm_75`。实际 `flags.make` 同时包含 `compute_120/sm_120` 与 `compute_75/sm_75`，`cuobjdump --list-elf` 也实测可执行文件同时含两类 cubin。证据：`03_attempt2_configure.*.log`、`05_attempt2_clean_rebuild.*.log`、`07_architecture.*.log`。

clean rebuild 从 17:36:52 到 17:38:37，退出码 0，生成 `build/weldformgpu-attempt2/src/WeldFormGPU.exe`。完整 stdout 约 175 KiB、stderr 约 17 KiB；存在大量警告，包括 C++14 已被 CUDA/Thrust/CUB 弃用、未初始化变量、缺少返回路径和 host/device 声明不一致。它们不是本次链接失败，但属于后续可信度风险。

`04_attempt2_clean_rebuild.*.log` 是日志采集器第一次传递带引号内联命令时的 Windows quoting 错误，发生在构建前，已保留且明确不计后端构建轮次。

## 上游固体案例运行

选择仓库原样 `examples/input/Compression.json`，未缩短时间、未改参数，并在隔离目录 `logs/upstream/weldformgpu/runs/compression/` 运行：

```text
build\weldformgpu-attempt2\src\WeldFormGPU.exe third_party\WeldFormGPU\examples\input\Compression.json
```

进程从 17:40:09 到 17:44:06。stdout 证明可执行文件实际看到约 25.65 GB GPU 总内存，建立了 11,316 个三维圆柱固体粒子，分配/复制设备数组、材料与边界 ID，并进入 `Initializing nb search data`。随后 stderr 记录：

```text
CUDA error 1 ... cub/util_device.cuh, 101: invalid argument
CUDA error 101 ... cub/util_device.cuh, 375: invalid device ordinal
CUDA error 101 ... cub/device/dispatch/dispatch_for.cuh, 188: invalid device ordinal
```

进程在错误后仍不退出，GPU 进程可见但 CPU 累计仅约 0.19 s，没有完成时间步，也没有生成 `out_*.csv` 或 `test.csv`。按有界 smoke 策略在 236.8 s 强制终止，因此日志中的 `EXIT_CODE=-1` 是终止结果，不是假称程序自行优雅返回。状态为 **BLOCKED_RUNTIME_NEIGHBOR_SEARCH**。完整证据：`06_compression.stdout.log`、`06_compression.stderr.log`、`06_compression.meta.log`。

CUB 的第一处报错在 `CurrentDevice()` 调用 `cudaGetDevice`；后续错误可能是前序 CUDA 错误传播，也可能是邻域库与当前运行时的兼容缺陷。本轮没有证据足以进一步归因，且两轮修复预算已经用完，故不继续试错。

## 工况/入口可用性

| 类别 | 上游资产 | 实际状态 |
|---|---|---|
| 最小固体 | `examples/input/Compression.json` | **ATTEMPTED / BLOCKED_RUNTIME**；构建通过，但在首次邻域搜索中止前无时间步、无结果文件。 |
| 热学 | `main_thermal.cu`、`cuda/Thermal.cu`，以及标为 `Mech-Thermal` 的 `MetalCut_mm.json`、`MetalCut_3D_AA_red_rac10um_fraser.json`、`Traction_JC.json` | **NOT_AVAILABLE_AS_ACTIVE_ENTRYPOINT**。CMake 中热学可执行目标被注释；主程序对 `solver`、热物性和初温的读取均被注释，`thermal_solver` 默认 false。JSON 名称不能证明热学实际启用。 |
| 刚体接触 | `Contact_Compression.json`、`Contact_Compression_2surf.json`、`cuda/Contact.cu` | **FOUND / NOT_RUN_BLOCKED_BY_SHARED_NEIGHBOR_INIT**。主 Fraser 循环存在接触 kernel 调用，但上游 README 明示 contact “in progress”，且共享邻域初始化已在无接触案例失败；重复跑接触输入不会绕过该阻塞。 |

因此没有热学或接触 PASS。热学缺少可从当前构建/JSON 激活的入口，明确记 `NOT_AVAILABLE`；接触资产存在，故不能误写成缺失，只能记 `NOT_RUN/BLOCKED`。

## 源码能力与实证分层

| 能力 | 已有代码 | 本次实际验证 | 仍需开发/修复 |
|---|---|---|---|
| 真三维邻域 | `cuNSearchKernels.cu:142-163` 遍历 x/y/z 的 27 个相邻网格，并用 `dx²+dy²+dz²` 判距；不是三维绘图替代。 | **未通过**；正是在 `SolverFraser.cu:109` 的首次 `find_neighbors()` 阶段出现 CUB 错误。 | 最小 cuNSearch 隔离测试、每次 CUDA 调用后的同步错误检查、修复当前驱动/toolkit/CC12 路径。 |
| 三维应力更新 | `SolverFraser.cu:204-218` 调用应变率与应力 kernel；`Mechanical.cu:458-551` 使用六分量张量、Jaumann 旋转项、J2 与三向剪切项。 | **未到达**；没有一个时间步。 | 单元测试刚体旋转客观性、单轴弹性、三维剪切和能量一致性。 |
| 塑性 | `Mechanical.cu:503-535` 有 J2 回归、塑性应变、Bilinear/Hollomon/Johnson-Cook 分支。 | **未到达/未验证**。 | `WeldFormGPU.cu:752-792` 创建 `Material_` 后没有给 Bilinear/Johnson-Cook 设置 `Material_model`；`Mechanical.cu:513-525` 的 `Ep` 在部分分支未初始化。必须先消除未定义行为，再做材料点测试。 |
| 热传导 | `Thermal.cu:68-99,101-142` 有导热/对流与 leapfrog 温度更新；Fraser 循环在 `:285` 调 `ThermalCalcs`。 | **未启用**；`thermal_solver` 默认 false，主入口配置代码被注释。 | 恢复有单位的热参数/初温入口；把 cuNSearch 邻域显式交给热核。目前 `neib_part/neib_offs` 在所审 Fraser 路径只分配，未见由 cuNSearch 填充。 |
| 刚体接触 | `SolverFraser.cu:233-257,291-318` 会调用表面、接触邻域、接触力和网格更新 kernel；`Contact.cu` 有三角面法向/穿透/摩擦计算。 | **未运行**。 | 邻域修复后建立单平面无摩擦/有摩擦 benchmark；再验证多表面。不能把 README/devlog 当 PASS。 |
| 积分器 | Fraser 实现存在。其他 solver 源文件也存在。 | 只验证到 Fraser 初始化前。 | `WeldFormGPU.cu:204-210,875-882` 的配置选择被注释，`:914` 无条件调用 Fraser；注释还写明 KickDrift 不工作。需显式、可测试的 solver 枚举及拒绝未知值。 |

## 候选决定与下一步

WeldFormGPU 当前不是“可工作的上游基线”。可以保留为源码参考和后续 GPU 修复候选，但在以下门槛通过前不应接 AFSD 三工况：

1. 用独立 cuNSearch 小程序复现并修复首次邻域调用，要求进程自然退出 0，且至少检查一组已知三维邻居。
2. 初始化并验证材料模型枚举及全部本构状态，增加弹性、Bilinear、Johnson-Cook 材料点测试；消除未初始化读。
3. 建立真实热入口，将温度/导热/热容与当前邻域数据连通，再做解析导热 benchmark。
4. 在无接触固体基线通过后，按单平面接触、多表面接触顺序测试；所有 CUDA kernel 加统一错误检查。
5. 把硬编码 `CUDA_ARCHITECTURES 75`、求解器选择和 Debug/Release 模式变成固定、可审计的构建配置。

由于 CPU WeldForm 已另有实测阻塞，而本 GPU 候选也为 runtime blocked，后续应按用户顺序进入单一 PySPH 候选评估；不要并行开发多套后端。

完整命令索引见 `logs/upstream/weldformgpu/00_exact_commands_and_status.log`。
