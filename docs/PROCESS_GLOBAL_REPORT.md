# 联合求解实际进展与失败边界（2026-10-03）

项目仍未结束。此次把已验证的结构响应和隐式热接到真实五场同史候选中，
并运行了一次固定500 ns联合求根；没有把最后未收敛候选发布为物理末态。
不启动GPU/长时生产，不修改物性、旧失败、原始接受阈值或自由材料运动。

## 已实现且实际测试的接口

- `stage_preconditioner.py::StagePreconditioner`：保留第一版直接密度候选修正及16独立回归。
- `stage_continuity_preconditioner.py::ContinuityStagePreconditioner`：凝聚压力—连续性反馈；
  10独立回归，包含实际27点旧算子的冻结弹性map及字面压缩解(.32 m/s,-1.2 kg/m³)。
- `stage_structural_root.py::solve_stage`：只用结构预条件器生成未接受候选；
  原始五场残差控制接受/增长保护，10独立回归，含零/微小proposal不能假PASS。
- `process_thermal_trial.py::ProcessThermalTrial, acceptance`：旧NNPS/B、FKS/Jaumann/J2、
  事件接触/执行器、feature guard不变；具名BE热反馈，逐事件摩擦证书、逐粒子/总热账。
  16回归PASS；首跑15 PASS+1 ERROR保留，修复的是制造云初始穿透，不是放宽门槛。
- `process_global_state.py`与`run_process_global.py`：新变体明确分类、严格父checkpoint、
  source/config/parent绑定成本、候选日志、事务提交与真实字段输出。5存储/预算回归PASS。

统一预检103项PASS：57新回归+46既有回归；不是103个新物理算例。
入口`check_process_global.py`，证据`logs/process_global_20261003/preflight_initial/report.json`。
初始旧freeze3464项PASS，旧代码/FAIL/轨迹保留。独立根与运行器审查纠正了
成本漏计日志IO和接受记录发布后才做拒绝性检查的问题，均在真实试跑前修正。

## 真实调用顺序

```text
同一accepted父历史
  → 候选 xmid/rhomid/Tmid + 独立事件平均速度
  → 原三维 NNPS / corrected梯度
  → 原 Jaumann/J2/FKS 材料更新
  → 原伴随内部力
  → 原多面事件接触与执行器反力（实际vbar≠简单端点平均）
  → 原有限feature guard
  → 各事件真实Q证书 / BE热图
  → 原始五场map残差
  → 结构-连续性Schur生成下一候选 / Anderson
  → 仅原始残差及完整物理检查都满足才可提交
```

结构预条件器的反力和修正不是实际外部冲量；BE不是旧显式热的同一时间离散。
没有把线性K称为完整FKS/接触切线，没有把未认证Es诊断当自由能。

## 一次500 ns运行的真实结果

父状态为`outputs/nonlinear_continuation_20261003/continuation/stage_004/accepted.npz`，
342点、t=850.393632910904 ns，完整塑性/应力/温度/ID历史保持。
目标t=1350.393632910904 ns，但**未达到接受态**。

| 项目 | 实测/状态 |
|---|---|
| 成本探针 | build .7774535 s；eval+proposal+完整日志7.4581122 s；外8.9168804 s |
| 事前32评估预算 | 389.155566 s，含1.5倍与30 s开销，低于450 s启动门 |
| 实际尝试 | 32评估，内265.000049 s、外265.5593475 s，均<600 s |
| 峰值工作集 | 199172096 B（约190 MiB） |
| 终止原因 | `StageFailure('structural root evaluation limit')` |
| 最终速度raw残差 | 7.255309545e-8 m/s；原门1e-9，FAIL |
| 最终事件均速raw残差 | 3.627653478e-8 m/s；原门1e-9，FAIL |
| 位置raw残差 | 1.813826866e-14 m；原门1e-12，PASS该代数项 |
| 密度相对raw残差 | 1.498962336e-14；原门1e-10，PASS该代数项 |
| 温度raw残差 | 1.023181539e-12 K；原门5e-12，PASS该代数项 |
| 新接受物理步 | **0**；accepted_history为空；失败candidate明确REJECTED |

实际命令/日志：`logs/process_global_20261003/COMMANDS.md`、`root_attempt.log`。
完整失败：`outputs/process_global_20261003/root_attempt/REJECTED_candidate.json.gz`，
32个数值候选及最后完整trial均保留。没有随后换dt、扩迭代数或放宽门槛重跑。

## 能支持的结论与不能支持的结论

原始速度缺陷从2.3182 m/s降低到7.2553e-8 m/s，说明结构/连续性候选修正有效；
这不是接受，更不是过程完成。此次没有新跨肩面质量，也不从REJECTED候选推导沉积结论。
现有有效轨迹仍850.393633 ns，之前实际中心穿越质量为0。

独立尺度审查`PROCESS_GLOBAL_STEP_AUDIT`已显示500 ns最高弹性模态omega*dt约10，
中点法强相位误差。即使该root收敛，也不能拿它直接做长过程。
原法向势能/Jaumann总能量差当前为诊断，不是通过宏步准确性验收。

## 下一任务与结束门

最终只读复核已完成，不调用被测trial/acceptance生成参考：

- `PROCESS_GLOBAL_HEAT_REVIEW`：14008接触子区间与8131热对局部账通过，
  candidate节点热残差3.224e-16 J；候选动量仍FAIL(2.449e-14 N s)，没有接受。
- `PROCESS_GLOBAL_NONLINEAR_DIAGNOSIS`：99.5357%的均速残差平方位于本kick有停滑事件的粒子，
  72.6664%位于多面接触。末候选所有点都在弹性侧，离屈服至少2.916 MPa，
  距率地板至少1.29475/s；没有feature增删，线性后向误差1.284e-16。
  这排除了末候选塑性/率地板/热账作为当前直接阻塞，支持接触响应漏项，但未做因果Jacobian隔离。
  最后残差仍下降；不要将预算耗尽误写为稳定停滞。

下一最小数值任务依据最终只读诊断，而非按函数名猜测；双面事件接触的局部响应
尚未进入结构预条件器。应先用独立反向/多面/执行器基准验证事件平均响应与端点响应，
再考虑仅改候选求解的联合Schur，不改真实摩擦冲量/热或原残差门。
独立反向参考为B=dvbar/dF=19/54 s/kg、C=dv1/dF=7/9 s/kg，B不等于C/2。
对未测试的切线原因保持SUPPORTED/POSSIBLE，不冒称CONFIRMED因果隔离。

并行于这项数值求解工作，过程路线必须解决：保留惯性但控制欠解析高频的具名积分变体及算法耗散账；
有限孔口/肩缘事件功；接受history绑定入口/完整源账。它们不能由更多相同短窗代替。
需要省略启动/去惯性则另给初始化与授权，不能清零旧checkpoint速度。

当前宏步同终时准确性、真实过孔及持续塑性输运、动态入口接受事务、
实际重算重启等价、三工况/论文比较仍NOT_RUN或未完成。阶段性数值进展不授予完整复现/生产许可。
项目结束仍按`PROCESS_COMPLETION_CONTRACT.md`，不按测试数量或自定百分比。
