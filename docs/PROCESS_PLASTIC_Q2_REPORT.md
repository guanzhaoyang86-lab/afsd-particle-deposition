# Q2：三维准静态塑性压缩—回程（2026-10-01）

## 结论和范围

四个独立小基准均到达规定的0.2s加载—回程终点，并保存真实三维位移、应力、
塑性应变、逐粒子塑性功和执行器反力。PySPH仍固定为
1cd3dea9cfaac7d9f1bda2a2491ce4b42a31dbe5；未改src、旧配置、旧测试或旧失败证据。
分类是ENGINEERING_APPROXIMATE_QUASISTATIC_SMALL_STRAIN_PLASTIC_BENCHMARK。
不是AFSD推进了0.2s，不是自由卸载后的残余形状，也不是论文复现完成。

本轮按用户允许粗精度的工程近似目标推进下一机制，不把空间误差继续调到全绿。
但保留有限值、试探回滚、力平衡和功账等基本检查；无质量/时间/模量缩放。
惯性、热、接触、进料和有限应变密度本轮均NOT_RUN。没有GPU或生产运行。

运行前定义见PROCESS_PLASTIC_Q2_CONTRACT.md及configs/benchmarks/process_plastic_q2.json。
所有完整算例均一次运行，外部timeout60s，未改变参数重试；命令/预算和日志在
logs/process_plastic_q2_20261001/COMMANDS.md。输出在outputs/process_plastic_q2_20261001。

## 实际调用与最小实现

    main -> PlasticBlock(QuasiStaticBlock -> E2 Block/PySPH ParticleArray/3D NNPS)
         -> inherited advance(t)
           -> solve_candidate(t)
             -> trial(u,t): accepted S/ep -> material.update(grad(delta_u)/physical_dt)
                           -> spherical stress K*trace(grad u)
                           -> paired mechanics.internal_force -> free residual
             -> tangent: 9-component local material finite differences
                       -> global assembled Kff -> bounded Newton/line search -> trial
           -> commit: recompute trial, check equilibrium/domain/virtual work
                    -> accept position/S/ep/plastic work/reaction once
         -> accepted history / upstream diagnostic NPZ / VTK / report+hashes

新scripts/probe_process_plastic.py继承已有Q1设施，不调用生产HostLeapfrog.advance，
也不把原显式积分器偷偷改为准静态。它复用原material.update和corrected梯度的伴随力，
固定参考构形的小应变近似明确命名。K的号与Fint=-L^T V sigma成对；不强行对称或pinv。
有体积、x/y/z独立自由度和完整3x3梯度/应力，不复制二维结果。

边长10mm；27/125点分别45/275自由DOF。E70GPa、nu0、rho2700kg/m3、h/dx1.3，
体积1e-6m3、总质量0.0027kg不变。三零坐标面的各自法向roller固定；X=L仅规定ux，
0至0.1s压到-40um，再到0.2s回原位。切向和内部自由度由平衡求解，无参考轨迹驱动。
这是双向位移执行器，允许拉力；尚不是只压不拉的刚体接触。

## 实测结果（不是实验误差表）

| 算例 | 粒子/增量数 | 完整进程墙钟 s | 终反力 N | 终体积平均ep | 累计塑性功 J |
| --- | ---: | ---: | ---: | ---: | ---: |
| j2_n3_i8 | 27/16 | 1.211684 | 9000.000000 | 0.004142857143 | 0.372857142857 |
| j2_n3_i16 | 27/32 | 1.916840 | 9000.000000 | 0.004142857143 | 0.372857142857 |
| j2_n5_i8 | 125/16 | 2.717260 | 7009.843543 | 0.003126008462 | 0.281340761536 |
| fks_n3_i8 | 27/16 | 1.856490 | 151.971287 | 0.007934901581 | 0.012050195460 |

全部DISCRETE_PLASTIC_PATH_COMPLETE、exit0；共同终时0.2s是本制造基准时间。
advance分别0.836427、1.523989、2.315397、1.459907s；不是步长稳定性或AFSD成本测量。
27点/125点selected_array_bytes分别235764/4913500；不含临时切线/求解缓冲或进程RSS，
峰值内存NOT_RUN，不能当完整内存承诺。过程运行前均已估成本，没有超过预算的启动。

| 算例 | Wact J | Wstress-Es-P J | 单独return诊断 J | Wact-Wstress J |
| --- | ---: | ---: | ---: | ---: |
| j2_n3_i8 | 0.427500000112 | -0.003214285711 | 0.103928571440 | 1.080e-10 |
| j2_n3_i16 | 0.430000000001 | -0.000714285714 | 0.053303571429 | 8.578e-13 |
| j2_n5_i8 | 0.324017191660 | -0.002268556743 | 0.103239054339 | 7.660e-10 |
| fks_n3_i8 | 0.011009812234 | -0.001056879819 | 0.206705846765 | -3.084e-11 |

Wact是实际反力—位移端点梯形积分，Wstress是同时间层应力功。P是原更新器的
后向端点塑性功V*q_new*dep。return=V*1.5G*dep^2独立输出，不是物理摩擦热，
不是总能量残差，更不加入温度。Es是原小应变诊断，不宣称已解决有限变形自由能。

### 独立参考与敏感性解释

J2制造材料Y90MPa，在单轴应力路径下独立分段公式给出：
屈服轴向应变Y/E=0.001285714286；压缩峰ep=0.002714285714；回程在
epsilon_x=-0.001428571429反向屈服；终sigma_x=+90MPa、ep=0.004142857143。
连续塑性功0.372857142857J，终Es=0.057857142857J，连续总功0.430714285714J。
测试参考由标量公式给定，不调用被测更新器构造期望值。

独立分段应力在8/16均匀节点上的端点梯形积分恰为0.4275/0.430J，解释n3的
上述功差来自跨屈服拐点的时间积分；不是凭空漏掉一份必须补入的热。
同终时步长减半，功差绝对值0.0032143降到0.0007143J。仅支持这个固定主轴小路径，
不消除旧Jaumann闭循环非零功反例，也不证明任意旋转/非比例路径正确。

空间检查保持几何、边界、温度、E/rho/Y、执行器和0.2s程序不变，仅n3到n5。
n5末反力相对均匀连续体参考低约22.11%，ep分布也不均匀，最大0.0103630。
这是已有边界/离散精度风险，本轮没有证明空间收敛；n3接近解析值不能外推为更细网格可靠。
保持记录并继续下一独立机制，不通过放宽旧精度阈值改判历史FAIL。

FKS明确使用原默认式，包括1/c2括号、b1无量纲DEBUG解释及1/s率下限。
规定T723.15K，E/nu仍是制造值，未标定AA6061。接受更新432/432次触发率下限；
末态最大实际偏应变率约0.04/s。其PASS只说明原FKS可进入整体平衡/事务接口，
并由独立公式读回P；不证明低率规律、温度演化、应变率敏感性或热力学闭合。
FKS本轮未做时间/空间敏感性，保持NOT_RUN，不能由J2的功差解释直接替它背书。

## 实际测试与失败记录

- 核心15项PASS：logs/process_plastic_q2_20261001/q2_tests_fixture_fixed.log，0.515s、exit0。
  包括独立J2反向屈服/零应力卸载、FKS公式/非法域、完整三维约束、切线方向差分、
  实际塑性步、physical_dt、试探/拒绝回滚、提交重核、过期状态和标签不影响力学。
- 初轮核心14PASS/1FAIL保留于q2_tests_initial.log。原测试错误假定n3均匀仿射场
  一定不平衡，实际已平衡。只将测试改为中心自由DOF局部扰动并先核实残差；
  没改driver、配置或容限让测试通过。
- 四个实际算例及其独立读回全部完成。逐粒子塑性功独立公式缺口均小于2e-17J；
  ID、源/数据hash、接受时间、执行器功、NPZ/VTK一致性实际检查。
- 读回回归还检查文件篡改、重hash后伪造塑性账、BLOCKED不能晋升COMPLETE、重复ID、
  未解温度字段。初版6项PASS保存于readback_tests.log；最终修订结果见下节。
- 完整旧116项本轮NOT_RUN，旧111PASS/5FAIL记录未改；旧raw、客观率/闭循环、
  接触失败和空间反例都保留。Q1旧13项本轮NOT_RUN，不冒用为本轮新测试数。
- 三维自由无载释放、单边接触、旋转多面接触/摩擦功、热、真实AFSD及生产均NOT_RUN。

### 读回最小修复

独立审查发现初版reader没有检查report.records与history帧数/step/time一一对应；
截短ledger有可能跳过材料功核验。只在新inspect_process_plastic.py补齐对应检查，
新增RED→GREEN回归并保留日志；不改任何物理结果。
readback_tests_red.log实际8项/3处FAIL（含两个subtest）、0.310s、exit1；
readback_tests_fixed.log实际8/8 PASS、0.306s、exit0。最终新增共15+8=23项PASS，
不重复计入早期6项，也不抹去历史失败。旧完整读回JSON仍保留；最终reader对四案例
重检全部PASS_COMPLETE_READBACK、exit0，另存*_readback_v2.json及同名日志。

## 剩余范围和下一步

允许进入下一小单元Q3：在独立小块上增加**刚性压头单边法向接触**，先替换本轮
能拉能压的双向端面约束；建立基于整体切线和真实边界约束的接触柔度/耦合。
先检验单面压入—脱离、反力/间隙/接触功及拒绝回滚，再做两个面的耦合检查。
不直接照搬显式kick中的单粒子P/m；不同时加入旋转摩擦、热、入口和整段AFSD。
下一轮预期改新scripts/probe_process_contact.py、小型配置/独立测试；复用Q2的
trial/tangent/commit接口，不先改src/material.py、mechanics.py、生产solver默认路径。

Q3未实现/未运行，须先写单边条件、独立参考与预算；这些小算例的秒级成本不证明
完整AFSD可行。固定参考邻域、有限变形本构自由能、移动接触和热耦合仍是结构缺项。
早期ValueError可只有REJECTED清单而无完整candidate数组；只有携candidate异常才保存
REJECTED_candidate.npz。当前四案例均无拒绝，但不能宣称所有异常均可完整恢复。
internal_torque_relative使用参考位置X，不是当前构形大转动角动量验证。
诊断checkpoint禁止生产restart，内含未求解scratch字段，不得包装成完整物理末态。

真实AFSD保持历史2us/191步、3个自由塑性粒子、无供料穿越/沉积；未在本轮推进。
生产三工况/长跑未授权。后续可视化只能读取冻结字段，不能改变参数或补造热/速度场。
