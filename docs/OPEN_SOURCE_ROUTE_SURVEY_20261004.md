# 开源过程求解路线调查（2026-10-04）

状态：DESK_RESEARCH / NOT_RUN。核对公开仓库、部分源文件、官方文档和研究论文；没有安装、编译或运行候选后端，没有新增物理接受时间。本报告不替代上游版本冻结、代码审计和独立复跑。

## 结论和本项目基线

有可借鉴的完整相邻工艺实现，也有针对接触或慢过程求解的开源模块；本次没有找到可直接满足全部约束的 AFSD 开源求解器。不能将“底层软件开源”当作“作者的专用求解器已公开”，也不能将同一仓库的不同求解分支功能相加。

当前交接已更新到 1.3503936329109039 微秒。500 ns 根求解成功，但与同起点的 250+250 ns 分支相比，末速最大差约 0.293670 m/s、应力差约 2.51499 MPa；两者新增塑性和肩面中心穿越均为零。依据本地 CONTACT_MATERIAL_SEARCH_REPORT.md；本次未重新运行。改用隐式软件只能改变求解能力，不能自动证明大步长精度。

优先级分两个目的：

- 获取相邻工艺的完整实现经验：先审计 FE-SPH-GPU 和 FS-SPH-GPU。
- 在当前不启动 GPU、不做质量/时间/模量缩放的阶段取得可承担的输运：优先评估独立初始化的热态黏塑性慢过程；OpenFOAM/rheoTool 有相应基础及 AFSD 研究先例，但仍需集成。
- 保留当前 PySPH 模型改进多面停滑：Siconos 最直接针对该数学子问题；它不能单独解决过程时间跨度与有限几何有效域。
- 隐式 MPM 为备选平台，不能仅因有 MPM、隐式和热塑性三个功能名就迁移。

这些是项目适配判断，不是已测得的性能排名。

## 1. FE-SPH-GPU：摩擦挤压，最接近“受力驱动穿孔”

[仓库](https://github.com/SPH-SSMP/FE-SPH-GPU)公开 C++/CUDA 求解器、构建说明和 PID 送料控制；GPL-3.0，README 指向 CUDA 11.0。研究对象是变形坯料、刚性容器与模具，不是 AFSD 的旋转方孔及可变形基板。

[2026 论文](https://doi.org/10.1016/j.aime.2026.100184)给出热力 SPH 与力控挤压，采用速度缩放 10、质量缩放 10^5，并省略初始加热阶段。作者报告 120 s 工艺在 A100 上约需 5 h；这是其算例结果，不是本项目预估。作者用动能/内能比支持缩放，不等于本项目的误差已受控。

建议借用/核查：模具孔口接触、摩擦热、有限几何、入口驱动与力反馈、输出及实验对照方法。不能原样继承缩放、热态初始化或材料参数。当前仅核对仓库说明与论文，尚未完整审计 CUDA 源码，也未建立当前 Windows 环境可构建性。

## 2. FS-SPH-GPU：摩擦堆焊，最接近“材料铺展和沉积”

[仓库](https://github.com/SPH-SSMP/FS-SPH-GPU)公开摩擦堆焊 SPH 实现（GPL-3.0，CUDA/CMake），与 FE 项目同源于 [iwf_mfree_gpu_3d](https://github.com/iwf-inspire/iwf_mfree_gpu_3d)。

[2025 论文，尤其第 4.2–4.3 节](https://link.springer.com/article/10.1007/s40571-025-01048-2)展示材料沉积及实验比较，但工艺没有 AFSD 空心工具。模型使用速度缩放 10、部分棒材/基板刚体切换及加入沉积层的判据；约 35 万粒子的算例优化后仍约需 A100 运行 50 h。其公开实现有参考价值，不能把其沉积判据或刚体切换当作本项目被动标签，也不能照搬为全域可变形模型。

建议借用/核查：有限变形材料流、棒材与基板接触、热交换、过程验证；每个数值稳定项和状态切换都需列入模型范围。未复跑、未审计全部源码。

## 3. OpenFOAM / rheoTool：热态黏塑性慢过程候选

这条路线有 AFSD 应用证据：[Kincaid 等 2023 论文的 ORNL 记录](https://impact.ornl.gov/en/publications/a-finite-volume-framework-for-the-simulation-of-additive-friction/)报告用 OpenFOAM 黏塑性流动模型模拟单层沉积，比较三种工况温度并定性比较沉积形貌。[作者博士论文](https://ir.ua.edu/items/a24a18ea-10e4-427e-9e3b-5b5c70a4d0be)说明组合现有求解器并加入定制代码。本次未找到该专用求解器的可核实公开仓库。

[rheoTool](https://github.com/fppimenta/rheoTool)提供非牛顿流动、热耦合及动网格基础（GPL-3.0）。源码核对的重要边界：

- [rheoHeatFoam.C](https://github.com/fppimenta/rheoTool/blob/master/of90/src/solvers/rheoHeatFoam/rheoHeatFoam.C)实际包含动量/压力、本构、热更新与被动标量。
- [rheoInterFoam.C](https://github.com/fppimenta/rheoTool/blob/master/of90/src/solvers/rheoInterFoam/rheoInterFoam.C)明确为等温两相 VOF 求解器。不能将它与前者合称已有完整非等温自由表面 AFSD。

拟议模型是将热塑性金属作为黏塑性连续介质，从力平衡、质量守恒和温度依赖流动应力求速度，再积分被动示踪轨迹；轨迹不预设。它可以另立慢过程范围，避免必须解析原预压初态的全部弹性瞬态；是否适用需先检查惯性、松弛和热时间尺度。后续源码审计更正：当前初温已为723.15 K，不能称为冷启动；慢过程的关键是独立力学平衡初始化，不是单纯升温。

主要研发缺口：非等温自由表面耦合、旋转方孔边界、工具相对滑动/摩擦热、低率本构处理、基板的变形区域与材料历史输运、入口质量和完整功热账。必须使基板参与力学，不能只把它作为导热壁面。转为热态独立初值，不能清零当前动态 checkpoint 速度后宣称连续。

这是当前工程近似目标的优先备选研究路线，不是已证明最快的后端。

## 4. Siconos / PETSc：针对当前耦合求解的模块

[Siconos 仓库](https://github.com/siconos/siconos)及[官方三维摩擦接触文档](https://nonsmooth.gricad-pages.univ-grenoble-alpes.fr/siconos/users_guide/problems_and_solvers/friction_contact.html)提供多接触 Coulomb 问题、非光滑 Newton、Alart–Curnier/Fischer–Burmeister、投影与 NSGS 等算法。

本项目可先从真实候选中构造接触相对速度—反力响应，联合处理多面粘着/滑动，再与材料和温度外迭代。必须先检查离散方程匹配：文档的标准约化问题有其矩阵性质假设，当前割线不自动满足；原罚接触不能无说明变成刚性互补接触。用同一力学问题比较残差、执行器功、摩擦耗散及成本，才可宣称改善。

[PETSc SNES 官方文档](https://petsc.org/release/manual/snes/)有 Newton、线搜索、预条件及带盒约束的变分不等式求解。盒约束 VI 不是现成的 Coulomb 摩擦锥模型；安装 petsc4py 或替换求根 API 不会自动消除当前耦合难点。

这些模块只针对求解，不提供穿孔几何、金属本构、热学或 AFSD 验证；当前均 NOT_RUN。

## 5. Kratos MPM：具备隐式架构，但热塑性分支存在缺口

[官方 MPMApplication README](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/MPMApplication/README.md)列出三维大变形、隐式 Newmark/Bossak 静态/准静态/动态、移动边界及刚体耦合。然而同页明确限定 Johnson–Cook 热塑性模型“just for explicit MPM”。

因此不能据此声称已具备隐式金属热塑性全过程。需要证明目标本构的一致切线、热方程、接触和入口能在同一隐式分支工作。对本项目属于更换离散平台的工程量；粒子化和隐式本身不保证足够大的准确时间步。

## 6. 其他已核对的候选

| 候选 | 确认的内容 | 对本项目的判断 |
| --- | --- | --- |
| [MOOSE](https://mooseframework.inl.gov/modules/solid_mechanics/) / [接触模块](https://mooseframework.inl.gov/modules/contact/) | 固体力学与接触、热传递模块可耦合；有三维摩擦接触例子 | 适合独立热态压缩/挤压基准；本次未核实可支持持续大流动的重网格/历史转移完整链 |
| [Uintah ImpMPM](https://github.com/Uintah/Uintah/blob/master/src/CCA/Components/MPM/ImpMPM.cc) | 隐式 MPM 源文件实际创建隐式导热与热接触；读取刚体接触配置 | 不能从显式分支摩擦功能推断隐式分支已有所需多面 Coulomb 热塑性耦合；暂列备选 |
| [Newton SolverImplicitMPM](https://newton-physics.github.io/newton/latest/api/_generated/newton.solvers.SolverImplicitMPM.html) | 隐式 MPM、弹塑性/黏性等参数 | 本次未确认金属温度依赖本构和热反馈完整链；不启动 GPU 评估 |
| [AFSD-Code](https://github.com/SSahaNU0923/AFSD-Code) | AFSD 热源与残余应力，附教程 | 并非目标材料输运求解器；见下述源码证据 |

[AFSD-Code 主程序](https://github.com/SSahaNU0923/AFSD-Code/blob/main/Tutorial/AFSD_new_model.py)读取既有网格并使用 element_birth/node_birth 激活材料、求热应变与弹塑性位移；可视化含位移放大，机械计算温度输入含 cp.clip。它有自己的热应力用途，但不能把激活得到的几何当成受力产生的料流，也不宜直接继承其处理到本项目。

## 下一步的有界验收建议（尚未执行）

1. 从 FE/FS 上游只读冻结目标 commit、许可和核心文件，建立本构/接触/缩放/材料切换差异表；提取可复用方法，保留当前 PySPH 为唯一活动后端。
2. 对慢过程备选先做最小三维热态方孔挤出基准，保留可变形受料基板和自由出口。热态是独立初值；初期可用明确标记的等温基准检验力学，随后必须接热反馈才能称目标过程。
3. 先测成本，再做每次不超过 600 s 的运行；成功条件是接受状态中持续正向穿越质量、非零新增塑性及孔下铺展，入口/出口质量与功热账可读回。
4. 至少有时间步减半和空间分辨率检查；数值收敛不等于实验验证。不得因刚体穿透、删守卫、材料切换或预定激活而获得假穿越。
5. 如果单位物理输运的成本仍不可承担，记录该路线失败；不能以更多组件测试替代过程交付。

本次改变：新增本报告和 HANDOFF 调查索引。后端、参数、守卫、checkpoint 和物理输出均未改变。
