# Q7R：有界求解恢复和接触事件诊断（2026-10-02）

本轮实现了命名变体Q7R_BOUNDED_NEWTON14，未重构/替换生产求解器。
旧12次Newton失败保持原样；新变体只增加到14次，不改变物性、力/滑移/功容限。
实际4增量程序仍 **BLOCKED**，有效末态仍制造时间.05s/step6，
并未推进AFSD沉积。通过原Newton停止点后，现有事件定位报cone event not bracketed。
不能把“求根恢复”称为“完整耦合修复”。

## 改动和真实调用链

- 新probe_process_plastic_friction_recovery.py：RecoveryBlock只改qcfg的Newton预算；
  临时构造器注入复用旧main/输出，finally恢复；没有复制/改变接触或材料算法。
- 新独立reader验证原16来源＋新增入口/配置来源，严格比较原配置，
  只允许12→14、变体元数据和4/8请求网格；再执行原逐粒子数据/功账审查。
- 原probe_process_friction_cone.py、probe_process_plastic_friction.py、src、
  配置/测试/失败和输出均未改。

实际链：advance_cone -> 小增量模态探测 -> branch目标求根 -> feasible_cone
-> 必要时strict二分定位 -> Q7 commit重新计算材料/接触、检查全部条件
-> 单次原子提交ep/P/R/Q -> accepted journal。
分支root合格不能直接commit；生产restart仍禁止。

## 实际执行清单

| 项目 | 状态 | 证据（logs/process_plastic_friction_q7r_20261002/） |
| --- | --- | --- |
| 原归档1429文件校验 | PASS | freeze_before.log |
| 新回归6项＋元数据测试12项 | 18 PASS，进程7.651s | recovery_tests_final.log |
| 原Q7解析/事务/标签/读回 | 21 PASS，13.563s | q7_compatibility.log |
| 原Q6兼容 | 31 PASS，1.293s | q6_compatibility.log |
| 新4增量完整程序 | BLOCKED，进程10.198s | n3_i4.log |
| 新输出独立审查 | PASS_PARTIAL_READBACK，exit2，.334s | n3_i4_readback.log |
| 8增量/空间加密/热/FKS/转工具/动态AFSD/生产 | NOT_RUN | outputs对应n3_i8/NOT_RUN.json |

初轮5个新核心测试PASS记录为tests_initial.log，最终不重复累加。
12项元数据测试是SCHEMA_FIXTURE_NOT_PHYSICAL，不是12个力学基准。
历史Q4严格ep FAIL和闭循环/空间FAIL未删除、未改判，本轮未重跑全历史套件。
每次子进程timeout60s，单线程，无依赖安装/GPU/长跑。

## 新回归实际证明什么

1. 原12次的19个REJECTED候选数组和全迭代trace精确重现。
2. 原同一接受态上14预算第13次root=1.637588335e-10；
   独立norm(tau)-.1*N识别ID21/feature1超圆5.691232926 N。
   该数值是冻结重放参考，不是独立连续体解析解。
3. 直接调用实际Q7 commit会拒绝这个候选；所有已接受材料/接触/粒子签名不变。
4. 新旧构造器除预算及元数据外完全一致；非法预算/强制类型转换被拒绝。
5. 实际短预压advance(.025)走现有探测/二分/提交，所有branch都收到14预算；
   有真实CONE_LIMIT，ID/标签、法向/支撑约束、原耗散/接触判据有效。
6. 原独立J2/双切向分段参考由旧Q7测试继续通过；未重新生成“参考答案”。

## 实际小块和成本

仍27粒子、9接触、54自由自由度；预估16请求区间＋事件、20–50s，
实际因事件失败提前停止于10.198s，不是外timeout。
有效记录0至.05s共7帧/6区间，与旧Q7全部30个history字段逐字节相等；
records/events亦相等。没有.05之后的接受事件。

已接受终态：max ep=.0037164738；塑性功P=.2396461244J；
摩擦Q=.00528460208J；Ws=.3012567643J；
Es=.06589598987J；Ws-Es-P=-.00428534998J仍告警，未转为热。
接触累计功账独立读回最大缺口4.35e-17J。
质量.0027kg；condM=3.385745；选定数组235764bytes不是进程峰值RSS。
准静态制造时间不是显式稳定dt、真实AFSD速度或论文尺度成本预测。

## 首个新阻塞：纯事件重放已完成

保存的REJECTED候选是目标.0625s、base step6，不是二分最后时刻。
异常iterations为空，不能据此声称二分里哪个节点/哪次比较失败。
因此仅补一次事前登记、禁止commit的事件重放，不是第二次完整物理运行。
进程2.571s、exit0；原37个分支全部返回root，37次feasible比较与独立向量公式相符。
19个目标候选数组逐字节重现；全部接受态指纹前后一致。
详细来源、独立prefix/功检查见同前缀RECOVERY_AUDIT.md；
完整trace在outputs/process_plastic_friction_q7r_20261002/event_diagnostic/diagnostic.json。

**CONFIRMED：无可行下界，不是二分宽度不足。**
在.05000125s模态探测中，耦合求解后的ID21余量为-6.095486e-5N；
普通9e-5N容限判可行，所以它仍为强制stick。目标.0625s余量则-5.691232926N，
必须被拒绝。34个strict中点全部False，ID21每次都是阻挡项。
最近中点.0500000000007276s余量-6.224624e-8N，普通容限可行但strict失败；
末宽度7.275985e-13s已小于1e-12s，lower仍None。ID24/26/27在近端也有负尾。
以上上下界由实际trace和未改源码推得，原异常没有直接保存这些局部变量。
无released-flow阻挡，未到commit；加Newton/二分次数不能解决已观测的缺失下界。

**SUPPORTED：应区别边界上立即起滑与区间内部事件。**
这个判断有全部已测中点支持，但不是对连续时间所有点的数学证明。
不可把strict一律放宽，也不可接受5.691N真实越界候选。
该事件修复及其独立新基准本轮NOT_RUN；原材料能量风险仍独立存在。

## 复现层级与下一步边界

代码实现和回归通过：有界预算接口、原子拒绝、来源/参数门禁、短预压事件路径。
完整三维塑性摩擦小块加载程序：BLOCKED；定量连续体/论文验证：未达。
真实AFSD仍原2us/191步，无持续沉积；不允许三工况生产。
下一次修改限于Q7 companion的初始接触分支分类和事件入口（参考
probe_process_friction_cone.py的solve_cone_modes/advance_cone，保留旧路径）。
新增测试入口应为tests/review/test_process_plastic_friction_onset.py，
先用独立切向弹簧k=1N/m、容量1N、初始tau=1N：
刚壁增量+dw时立即滑动，tau=1N、du=0、Q=(1N)*dw；
反向小增量-dw时先stick，tau=1N-k*dw、du=-dw、Q=0。
再加入耦合多接触和明确标记的容限内带符号起始扰动，不能用被测分支生成参考。
新接口须区分BOUNDARY_ONSET和INTERIOR_EVENT，保留相同最终commit验收；
没有独立基准证据前不实施此提案。
下一次修改必须基于事件trace和独立“边界起始/多接触方向变化”小基准，
保留当前Q7R失败，不能继续加Newton次数或把越圆候选当有效态。
命令见logs/process_plastic_friction_q7r_20261002/COMMANDS.md。
