# 完整沉积所需的下一项具体集成判断

**延长当前 342 粒模拟不能把它变成完整沉积模型。** 当前小域可以继续给出有限料棒的
出口穿越和局部塑性响应，但 4 mm 基板初始即不覆盖 5.4 mm 肩径，有限 eligible 库存
没有补料；时间步仍出现局部应力比较失败。下面仅为实际源码核对所得的可实施改进，
未修改冻结文件，未运行初始化、NNPS、物理候选或新测试，不声称入口已经接入。

## 先建立更大基板的独立热态分支，不能给当前 checkpoint 补一圈粒子

可直接复用 `src/afsd_backend/pysph/particles.py::create_geometry`：它按 dx 的整倍数
生成方形基板及方棒，`make_array` 已生成完整物理字段、持久 ID、初始位置和被动标签。
当前 `scripts/run_hot_quasistatic_probe.py::parameters` 固定读取
`outputs/contact_material_search_20261004/root/resolved_config.json`，仅覆盖孔宽为
2.504 mm。各活动 controller 又比较参数和来源 SHA；新域需要另一个有来源的参数
入口，不能改该冻结 JSON 让旧 checkpoint 冒充同一模型。

独立初始化可以复用 `quasistatic_process.py::independent_initial_state`、
`discrete_process.py::independent_initial_state` 的无载热态、底面固定、棒顶仅规定 z
进给、八面接触 history 以及初始出口资格逻辑，再为新域建立明确的 visible-support
初始资格。`visible_quasistatic.py::explicit_transition` 当前针对有旧分类/父 SHA 的
显式转接；不能伪造旧过程父态来获得新域资格。

按原 dx=0.5 mm、基板厚 1.5 mm、原 150 粒料棒和 25 个棒顶 z 约束计算：

| 方形基板宽 | 总粒子数 | 自由位移 DOF | 当前 dense D 单份 | 初始肩边到基板边 |
|---|---:|---:|---:|---:|
| 4 mm（当前） | 342 | 809 | 24.09 MiB | -0.7 mm |
| 6 mm（较小覆盖候选） | 582 | 1289 | 69.77 MiB | 0.3 mm |
| 8 mm（较宽域对照候选） | 918 | 1961 | 173.60 MiB | 1.3 mm |

这些是按构造公式计算的点数/矩阵字节，**没有执行扩大域案例**。6 mm 只解决初始
几何覆盖；当前核邻域半径为 2h=1.3 mm，0.3 mm 余量不能证明边界效应足够小。
8 mm 也只是域对照候选，不是已验证的无限板替代。当前工具从方板中心向 +x 走，
静态未变形几何必须满足 `W/2 >= R + traverse*t + buffer`；因此域宽还必须由拟做
短轨迹长度确定，不能只看初始覆盖。

最小独立验证范围应是：相同 dx 的 6/8 mm 两个短域，先核对初始无交叠、真实 NNPS、
修正矩阵秩/条件和固定自由模态；然后一次小载荷完整候选及能量/热/接触账、同父
粗细短窗比较，并比较肩区反力、温度和塑性响应的域宽敏感性。真实运行前分别测
完整单根成本、峰值内存，满足原 450 s 预测与 600 s 外部限制后才能延长。
首个更大域资格不是三工况生产，也不能继承 342 点已有接受时间。

## 有限工具可以复用；单体距离势能组件不能直接替换当前多面力

当前真实链已经使用 `quasistatic_contact.py::QuasistaticContact._patches/geometry`、
`discrete_contact.py::ChainRuleDiscreteContact.evaluate` 以及
`scripts/batched_contact.py::BatchedChainRuleContact`。其有限肩面、顶面、开放方孔、
外圆柱和底面进入真实受力；`visible_quasistatic.py::evaluate` 将同一过滤图用于
密度、梯度、力和热，`visible_support_path.py` 与 `contact_shell_finite_enclosure.py`
检查整步路径。新域应复用这条已经明确标注的力学路线。

较早的 `finite_tool_geometry.py::FiniteOpenSquareTool.query` 是独立几何/势能组件，
文件明确写着未集成、其单体最近距离势能不同于逐面弹簧之和；它可作几何独立参考，
不能作为同物理公式的替换。现有 TOP/SHOULDER 白名单也不能因新域而直接放开
OUTER_CYLINDER/BASE 新接触；出现相关事件仍需独立有限步功和真实同父事件比较。

## 持续入口有可复用协议，缺的是物理接受事务

`inlet_protocol.py` 已有明确可复用的接口：

- `initialize_registry`：持久 ID high-watermark、事件号和累计来源账；
- `propose_insertion`：复制父 PA 后追加新粒子、保留旧粒子所有 history、分配不复用
  的新 ID，逐字段核对单位/来源/初始条件；
- `restore_registry`：PA/registry 同时重启、密封摘要校验；
- 来源账已算质量、动量、显热、动能和体积，但事件明确为
  `INLET_PROTOCOL_CANDIDATE_NOT_ACCEPTED_NOT_INTEGRATED`。

`tests/review/test_inlet_protocol.py` 具有上述账、历史纯性、排序、两次插入、重启、
标签无力学影响、重复事件和 ID 溢出等制造检查。此次只核对测试源码，没有重跑。
`particles.py::append_feedstock` 及 `tests/backend/test_gate_c.py` 的 inflow 用例只检验
追加和身份跟随，不能承担过程入口的物理接受。

缺失接入位置是明确的：

1. **registry 随物理接受更新。** registry 的 `particle_state_sha256` 密封全部 PA
   物理字段；真实一步后位置/应力等变化会使原 `restore_registry` 拒绝。需要新增
   有版本的接受事务层，保持 allocator/event/source 累计历史，同时与新的真实
   accepted PA 原子发布。不能每次按 max(ID) 重建 registry 或抹掉旧来源账。
2. **入口力学与开口功。** 当前协议把 pressure-flow work、elastic energy 和
   mechanical reaction work 显式记为 `NOT_COMPUTED`。需规定入料断面的应力/密度/
   温度/history，补齐入口牵引功、对流储能及体积源的总账，避免将压力功和执行器
   反力重复计入。当前 QS 不含惯性储能；协议中的动量/动能账不能自动变成它已
   求解惯性，必须明确在该近似中的地位，不能将缺口补成热。
3. **全局状态形状与执行器。** `contact_history` 是 `(N,8,3)`、`free` 是 `(N,3)`；
   必须按持久 ID 延展并保留旧行。当前 actuator 由初始化高度一次标定，
   `QuasistaticTrial.prescribed` 永久给这些点规定 z 增量；反复 append 会留下旧
   执行器层并新增受控层。需要明确入口控制区域的接入/释放规则及反力功，不能让
   已离开入口的材料继续被隐含驱动。新粒子标签仍只观察。
4. **几何/邻域/材料接受。** 每个入口候选必须做无交叠、孔口有效域、真实 NNPS/
   visible graph/修正秩、接触 history 初值和完整热材料求值。插入使拓扑改变，
   必须重建切线，不对插入前后做 Broyden secant。未完成求值的 scratch 字段仍
   不能标成真实结果。
5. **观察器与读回。** `interval_acceptance.py::observe_square_bore_crossings` 明确
   要求 ID/component/mass 集合不变；现有 trajectory/局部比较也依赖固定 N 和
   ID 顺序。需要有版本的出生事件、出生时上方资格及可变粒子集读回；新增库存
   不是出口穿越。不存在的出生前帧不能填成伪轨迹。细/粗入口时刻和 ID 分配还
   必须一致，才能做同父时间比较。

最小集成验证应先在独立的三维短棒入口基准中做**两次**插入，中间至少经历一次真实
机械/热接受，并在两次之间重启一次；核对旧 history、ID、不重复来源账和边界功。
再与同域内预先显式放置较长料棒的有限参考比较出口断面响应及来源收支；入口处给定
边界条件，内部材料运动仍由力学产生。通过后才在上述更大基板新分支接入口，保留原
342 点分支为独立证据。当前代码没有可直接启用的成熟持续 inlet runner。

## 扩大域前最实际的数值成本改进位置

`quasistatic_tangent.py::gradient_map` 分配 dense `(9N,3N)`，虽有邻域局部填项，
并非稀疏存储；`assemble_frozen_force_tangent` 用局部块乘法减少计算后仍分配
dense `(3N,3N)`。`secant_equilibrium.py::_full_matrix` 另取自由矩阵及其完整逆矩阵，
`visibility_equilibrium.py::solve` 每次方向仍 `np.linalg.solve`；
`implicit_thermal.py::advance_implicit_thermal` 也组 dense `(N,N)` 并直接求解。
本次在实际 PySPH 隔离环境检查 `find_spec('scipy')` 为 None；项目活动源码未找到
已实现的 CSR/稀疏解法，不能说换一个配置即可启用。

6/8 mm 的自由矩阵三次复杂度相对当前约 4.04/14.24 倍，**仅为运算量量级，不是
实际耗时预测**。342 点的 46.5–47.5 ms 粗细实测外部耗时为 51.419+146.225=
197.644 s，不能按点数线性外推。旧 profile 的 residual 0.148 s、tangent+linear
0.774 s 来自更早状态，只能说明其历史位置，不能代替当前扩大域测量。

若新域单根成本超预算，最小数值改进应是单独适配热矩阵及局部力学矩阵的稀疏存储/
求解，保留 PySPH 邻域和全部物理公式。先用冻结的同一实际图/材料/接触状态对比
矩阵作用、自由轴方向、线性后向误差、完整残差及原能量账，再跑小型实际短窗；不能
仅比较求解速度。稀疏 Broyden 还需避免 rank-one 更新将矩阵重新变密，现有完整逆
条件监控也不能直接删掉，必须另定可核验的数值资格。这是尚未实现的局部求解器适配，
不是已经拥有的新后端，也不需要启动 GPU 项目。

以上改进均应添加新参数来源、初始化/入口/求解器适配文件和独立输出，不就地改冻结
源码或现有接受历史。当前再得到若干毫秒的穿越响应，也仍须按其有限小域范围报告。
