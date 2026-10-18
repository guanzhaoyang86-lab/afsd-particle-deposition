# 独立审查覆盖通知（2026-09-24）

当前权威为SCIENTIFIC_REVIEW.md和VALIDATION_CONTRACT.md：不允许三工况生产。
新增独立38项测试35通过3失败；修复后的四档耦合敏感性试算均负中点摩擦耗散失败。
下文是本轮修复前的历史实现报告/日志，不删除，也不能当作当前源码的通过证据。

# 实现与实测状态（2026-09-24）

结论：已实现并实际运行 **PySPH上的真实三维热弹塑性 DEBUG 耦合探针**。
不是动画；不是论文尺度AFSD、稳定沉积或论文复现。未启动长时间/付费任务。

唯一后端：PySPH 1cd3dea9cfaac7d9f1bda2a2491ce4b42a31dbe5，BSD-3-Clause，
上游已跟踪源码干净。没有再比较/开发WeldForm。复用真实 ParticleArray、dim=3
LinkedListNNPS、LeapFrogIntegrator.one_timestep、LeapFrogStep.stage1/2、
RigidBodyWallCollision.loop、NPZ dump/load；新材料/热/应力等方程在项目主机适配层，
**不是已编译的PySPH Equation流水线，也不是GPU实现**。

## 可重跑入口和证据

在项目根目录运行（PowerShell 7；使用已隔离的.envs/pysph）：

```powershell
pwsh -NoProfile -File scripts/run_physics_checks.ps1 -RunId fresh_run_id
```

RunId不得重用。脚本运行原19项骨架回归、21项数值/身份测试、三工况各200步、
半dt的400步、重启20步、墙钟上限和故意密度告警、全时序回读。
每次求解默认600s上限，最多10000粒子；墙钟在步边界检查，结束输出有额外开销，
不是可中断任意原生调用的硬实时计时器。当前试跑最长约11秒，没有分钟级长任务。

最终证据根：

- logs/physics/20260924_physics_verified/commands_and_status.json：精确解释器、argv、预期/实际退出码。
- 同目录各项.log：真实输出。脚本整体exit0/CHECKS_PASS。
- outputs/debug/20260924_physics_verified/：与当前源码哈希一致的正式DEBUG证据。
- 早期probe/final/release等目录保留迭代历史，不能当当前可重启版本；早版初始导数字段为未评估的零占位，最终verified版本已经删除这些初始字段并记录原因。

单独运行/继续：

```powershell
.\.envs\pysph\python.exe -B scripts/run_physics_tests.py --gate all
.\.envs\pysph\python.exe -B scripts/run_afsd_debug.py --case configs/debug/starved.json --output outputs/debug/new_starved --steps 200
.\.envs\pysph\python.exe -B scripts/run_afsd_debug.py --case configs/debug/starved.json --output outputs/debug/new_restart --restart outputs/debug/20260924_physics_verified/starved/checkpoint.npz --steps 20
.\.envs\pysph\python.exe -B scripts/validate_physics_outputs.py outputs/debug/new_starved
```

NPZ使用上游pickle格式，只允许可信本地生成文件。源码/配置改变后拒绝旧checkpoint，不能绕过后声称连续重启。

## 四层验证状态

| 层次 | 本次事实 | 尚不能声称 |
|---|---|---|
| 代码实现 | 3D核/邻域/应力、FKS-J2、Jaumann、线性EOS、导热/功热、有限面接触、执行器、时序/ID/重启 | 完整三角面粘滑、真实启动、温变热物性、无限供料已实现 |
| 测试通过 | 21项后端测试+19项原回归；三DEBUG短算例；故障检测与回读 | 全工况鲁棒、网格/时间收敛或论文参数已验证 |
| 定性对照 | NOT_RUN；只观察短时接触、塑性和热变化 | 与Fig.8/11/12材料混合/沉积分布一致 |
| 定量验证 | BLOCKED；缺材料/边界来源、尺度和数据；当前塑性时间步敏感 | 已复现论文/实验误差范围内成立 |

## 关卡测试顺序

A先实测通过后才增加B；B通过后增加C；C通过后增加D及真实小算例。
最初日志gate_a_01、gate_b_01、gate_c_01、gate_d_01均exit0。
最终21项：A8、B5、C5、D3，详见tests/backend/和最终physics_gates.log。

主要数值判据：核3D归一化1e-8、分段导数1e-8；真实非共面NNPS对暴力金标准；
仿射内部梯度相对2.5%；成对合力1e-11；预应力任意轴旋转相对2e-12；
三剪切面返回残差相对1e-10；体积/偏量分离及弹性卸载；FKS独立标量式/域检查；
导热两体解析衰减、守恒和3D Fourier空间细化；三维接触功热等式1e-10；
非零9分量应力随ID重排；重启逐字段逐位相等；标签改99/77仍力学热学逐位相同；
leapfrog常加速度精确、谐振子位置误差缩减比>3.9；NaN显式失败。

这些接触测试不等于已通过完整静粘滑/三角网格边缘测试；热Fourier局部周期模板不等于
任意自由表面导热都已收敛。塑性热转换与耗散账本实现了，但耦合能量还需要进一步收敛。

## 真实小算例与性能（最终verified目录）

三案均657粒子（432可变形基板+225可变形方棒）；15,879无向邻居对量级，
dx=.5mm，层高1.5mm=3dx。共同300rpm、走速2.12mm/s；方孔/肩面/底座均刚性。
几何缩小、预热、过盈及所有参数偏离见PARAMETER_SOURCES.md D01–D10。

| 工况 | 步数 | 物理时间(us) | dt约(ns) | 墙钟(s) | 最终能量残差 | 最大塑性应变 |
|---|---:|---:|---:|---:|---:|---:|
| starved 1.06mm/s | 200 | 2.09794 | 10.49 | 5.5874 | 2.1231% | 3.24325e-5 |
| copper-core 2.12mm/s | 200 | 2.09794 | 10.49 | 5.5535 | 2.1094% | 6.04718e-5 |
| overfed 4.24mm/s | 200 | 2.09794 | 10.49 | 5.6730 | 2.1087% | 4.40918e-5 |
| starved半dt | 400 | 2.09795 | 5.245 | 10.9839 | 0.8820% | 5.46696e-5 |

starved最终质量0.0002217375kg，无外部新增质量；最大相对密度变化1.41233e-4，
穿透1.97868e-6m（约0.396%dx），最大速度0.69767m/s。没有NaN；所有D10门限通过。
塑性、摩擦、接触热三类实际非零，分量逐步保存；没有温度/屈服裁剪或残差补偿。
初始过盈释放主导这一微秒过程，不能把三案的微小差异解释成供料沉积规律。

半dt降低能量误差，但最大塑性应变变化约69%，速度也变化：**尚未时间步收敛**。
未用放宽容限把失败包装为通过；故意将密度门限设1e-12的独立故障试验第1步exit1，
manifest=FAILED，保留第0步有效checkpoint。wall_cap=.05s试验第3步停，exit2/INCOMPLETE；
实际约.097s含结束保存，符合步边界软停止的已声明行为。重启从200到220步exit0。

主机方程中位步耗时约.0198–.0203s，含诊断/输出的循环约35–36步/s；进程峰值工作集
约63.4–63.7MB（Windows工作集，不是GPU显存或完整提交内存）。wall_seconds从初始快照后
计时到循环/最后保存结束，不含Python导入/几何构建/初始快照；performance字段中的
end_to_end_steps_per_second仅指这个循环口径。

按N线性放大到80,312粒子，并以论文名义dx=1mm把dt翻倍的**乐观外推**约
3.2万小时/1物理秒（dt-half约6.5万小时）。未实测大规模，也未启动该计算；
接触比例、邻域、内存及并行效率变化会改变成本。名义1mm只有1.5层间距，不满足
本DEBUG至少3dx的最低要求，不能拿该外推当可用论文配置。

## 明确阻塞与下一关

1. 先检查当前一阶材料/热分裂和接触切换引起的能量/塑性步长误差；dt/2、dt/4、解析接触振子和更严格机械功共轭测试，未收敛前不延长物理时长。
2. 补P40/P41可信全文及材料适用域；D08 floor和热参数必须有范围敏感性，不能升级paper标签。
3. 三角面邻域/边角、静粘滑历史、有限体热分配、变形体接触/分离/结合、真实夹持和入口守恒仍缺。
4. 几何真实尺度、预热/驻留/抬升程序、h/人工稳定化/边界梯度校正仍需来源和独立测试。
5. 将已验证标量方程移入同一PySPH编译Equation路径并做逐步数值等价测试，重测成本；不要切换/开发第二后端。
6. 以上通过后才规划更长的3D定性试验；实验ROI、截面积指标和定量数据仍未具备。

