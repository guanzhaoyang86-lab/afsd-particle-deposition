# Q6：双切向圆形摩擦锥受限基准报告

2026-10-01。结论：本轮实现及独立小基准通过，可进入下一步“J2塑性＋双切向摩擦”
的最小耦合开发；**仍不允许AFSD三工况生产运行**。生产src、旧配置、旧tests及旧FAIL
未改。PySPH固定1cd3dea9cfaac7d9f1bda2a2491ce4b42a31dbe5，没有更换后端。
本轮没有重新跑失败的真实AFSD，也没有将它伪报为通过。

## 实际增加的能力与范围

新增companion scripts/probe_process_friction_cone.py；Q5只在Y方向有摩擦、Z自由，
Q6在YZ平面联合解 ||tau||<=mu*lambda 的圆形界，不是逐分量矩形限幅。
接触力方向来自解出的材料相对壁运动，不能直接用壁速替代。共享自由度及当前法向
反力通过整体切线耦合，既不使用独立粒子质量限幅，也不作自由粒子轨迹赋值。

算例是1cm立方体、27三维材料点/9接触/54自由DOF；参考点三维NNPS/算子来自已有
设施，未沿厚度复制二维解。X=0整面三分量夹持，X+刚壁先压入40um再保持，
切向幅值20um沿全局(Y,Z)=(.6,.8)，0→正峰→负峰→0。E=70GPa、nu=0、rho=2700，
制造弹性yield=1e12Pa、mu=.1、规定温度293.15K，全部是已有基准/显式制造参数，
不是AA6061拟合、FKS或论文实验条件。层厚、旋转方孔、供料并未用这块体近似替换。
此companion求准静态平衡，不调用旧动态kick；.2s是制造加载程序时间，不能据此声称
真实AFSD从2us前进到了.2s。源方程/PDF解释本轮未改。

实际调用顺序（与旧动力学求解器分开）：

```text
ConeBlock.advance：阻止跨加载拐点
  → advance_cone：accepted态纯探针/粘滑首界二分
  → ConeBlock.branch：约束法向和候选粘着DOF
    → PlasticBlock.trial：三维增量梯度 → 原材料更新 → 配对internal_force
    → ConeBlock.state：当前lambda/tau/slip、自然映射及力平衡残差
    → 实际material tangent + cone_rows → 整体Newton/有界线搜索
  → ConeBlock.commit：从accepted重新trial；接触定律/约束/逐接触功验收
    → 一次原子提交材料和接触 → accepted时序/VTK/诊断checkpoint
```

温度与惯性方程没有参与；摩擦耗散没有注入残差热。tau=-Fint_yz会让接触处机械
平衡恒等成立，因此提交另验自然映射根、圆界和滑移方向；不能单凭力平衡放行。
盘投影只用于求根残差，实际力没有事后投影裁剪。kappa=E*L是数值单位尺度，
独立测试改变.1/1/10后根不变；没有改变质量、模量或物理时间。

## 独立参考与实际测试

参考推导、单位和运行前容限见同前缀CONTRACT及REFERENCES，非调用被测函数生成。
关键字面参考：单位弹簧、单位圆容量、e=(.6,.8)，壁0→3e→0，终u=-e，
Q=3J、Wwall=3.5J；两次反向0→3e→-3e→0，终u=e、Q=7J、Wwall=7.5J。
两个接触共享自由度的基准总Q=4.5J、Wwall=5.25J、储能.75J；不能用两个孤立弹簧代替。
各向异性K=diag(2,1)、壁(1.5,2.4)时端点u=(.3,.8)、tau=(.6,.8)、slip=(-1.2,-1.6)，
它验证方向不是壁方向，但不是连续任意曲滑路径的积分参考。

| 项目 | 实际状态 | 证据（logs/process_friction_cone_q6_20261001/） |
| --- | --- | --- |
| 15项独立核心：斜滑/反向、共享自由度、法向容量、非共线、基底、kappa、标签/回滚 | PASS | q6_tests_all_final.log |
| 5项guard：接触定律、完整Jacobian、mask、拐点/预算、近界误差传播 | PASS | 同上 |
| 9项I/O对抗测试 | PASS，仅schema保护 | 同上、readback_tests.log |
| 2项比较保护：不同物理/事件冒充请求帧/部分终态 | PASS | 同上、comparison_tests.log |
| 新测试合计31项 | PASS，完整进程1.209s | q6_tests_all_final.log |
| 原Q5三组测试23项 | PASS，不重复计新测试 | q5_compatibility.log |
| 三组真实输出独立读回 | PASS_COMPLETE_READBACK | readbacks_final.log及outputs下*_readback.json |
| 原Q4严格ep尾差、原材料闭循环/空间失败 | 本轮NOT_RUN；旧FAIL保留 | 父Q5/Q4归档及原日志 |
| 物理刚体旋转且非零初应力、热、J2+摩擦、空间n5、动态AFSD | 本轮NOT_RUN | 未用基底旋转或I/O夹具冒充 |

9项reader测试使用显式SCHEMA_FIXTURE_NOT_PHYSICAL的合成I/O；真实求解证据来自下面
三个独立运行目录。reader不导入被测接触函数，独立重建壁程序、向量滑移/根/功与Es，
但没有独立重建材料点梯度或求解连续体，因此不等于独立空间/本构验证。

## 真实尝试与最小修复（失败不删除）

1. core_initial.log：15项中14PASS/1ERROR。共享双接触反向界后探针t=5.0001s，
   实际tau≈(-.6,-.8)N，slip≈(4.5e-5,6e-5)m，自然根2.22e-16N。
   旧方向检查除以||slip||后得1.18445e-12N，超过1e-12N检查线。
   core_shared_failure_diagnostic.log保留accepted/candidate和完整分量。
   CONFIRMED：方向检查放大位移相减的舍入误差；不是参考力或摩擦根错误。
   最小修复在feasible_cone：检查 || ||s||tau+C*s ||，容限为原ft*||s||+C*st。
   没改ft/st/参考答案/物性，不abs(Q)，不裁负耗散。近零滑移方向精度依然受st限制。
   contract附录明确这是测试后的误差度量修正，原ERROR依然在案。
2. guards_initial.log：4方法3个失败条目。两个非bool mask(nan/int2)原会强转true；
   release_inventory现在只接受bool。另一个是新增guard自身选.075s却尚未超圆的
   测试设定错误，改为原程序既定峰值.1s；不改模型或加载参数。
3. 力平衡但超圆候选、反向/垂直错误滑移力、过期候选、损坏事件/ID、隐匿未求字段
   均有拒绝回归。纯trial不改变accepted状态；标签修改不改变力学trial。

## 三组实测、成本和敏感性

全部输出在outputs/process_friction_cone_q6_20261001/，均首次真实运行即COMPLETE：

| 目录 | 每相请求增量 | 接受区间/事件 | 共同终时 | 完整进程墙钟 | 总摩擦耗散J |
| --- | ---: | ---: | ---: | ---: | ---: |
| oblique_i4 | 4 | 24/8 | .2s | 15.902s | .00109144760190150 |
| oblique_i8 | 8 | 40/8 | .2s | 16.692s | .00109146692254245 |
| basis37_i4 | 4 | 24/8 | .2s | 15.719s | .00109144760190150 |

三组均由外部subprocess.run(timeout=60)限制全进程，含初始化/输出；未启动600s以上
任务。预算和完整复跑包装见日志COMMANDS.md。selected arrays235764bytes，不是
进程RSS或峰值内存。M最大条件数3.386；没有pinv。准静态程序无声学CFL稳定步意义，
这些墙钟不能外推真实AFSD论文尺度生产成本，后者仍BLOCKED/未测。

最大自然接触根8.4482e-5N，小于既定9e-5N；这不等于机器精度，不能用极小的机械
平衡残差代替它。独立逐接触总功收支最大约4.16e-16J。原始最小单接触区间Q约
-4.88e-19J照存；两点累计Q亦有负舍入尾差，未裁成0。三组最终ep=0。
末态仍受40um法向压缩；弹性储能约.560118J非闭循环残差。切向受载残余不是塑性变形。

comparison.json只在17个保存的共同REQUESTED时点比较，不插值事件帧、不比较失败终时。
4→8终态最大位移差4.762e-10m、应力分量差1811Pa、单接触切向力差.00784N；
总摩擦耗散差1.932e-8J、壁功差1.610e-8J。8个匹配事件的最大时刻差6.774e-7s。
这是两档增量敏感性测量，不给收敛阶或空间精度结论；没有为减少差值继续调参。
37°基底终态位移差1.097e-20m、壁功差0，满足预定2e-9m/1e-6J。
它保持同一全局边界，只检验坐标表示，不证明旋转工具或Jaumann非零应力客观性。

## 文件、归档与下一步

新源码：probe_process_friction_cone.py、inspect_process_friction_cone.py、
compare_process_friction_cone.py、freeze_process_friction_cone_q6.py；新配置
configs/benchmarks/process_friction_cone_q6.json；4份test_process_friction_cone_q6*.py；
本报告/合同/独立参考、时序NPZ/VTK/诊断checkpoint、输出读回/比较和所有日志。
唯一允许修改的旧文件：HANDOFF、IMPLEMENTATION_PLAN、OUTPUT_SCHEMA。
归档通过PROCESS_FRICTION_CONE_Q6_FREEZE.json和preservation_report保存源/配置/输出
哈希及旧失败；实际最终检查结果以freeze_final.log为准，归档本身不授予物理PASS。

下一Q7新增小范围塑性-摩擦companion：复用Q2实际J2 trial/tangent与Q6双切向根，
统一原子commit和塑性/摩擦功分账。先独立弹塑性材料点/串联摩擦基准，后同几何小块，
单次预算仍<=60s。明确是常屈服J2制造变体，不静默冒充FKS。先不加热、供料、转工具。
旧ep严格尾差、有限应变闭循环、空间自由面风险不可被本轮PASS覆盖。
当前级别：工程近似方法组件＋独立基准；论文定性/定量对照均未新增。
真实AFSD仍2us/191步、无持续沉积；生产三工况不放行。
