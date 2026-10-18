# 连续推进与真实输运：本轮结果（2026-10-03）

## 结论

**原来整体Picard不能接受下一步的阻塞已跨过；真实过程尺度材料输运尚未完成。**
同一342点三维问题、同一原史/物性/空间/接触/热模型、同一接受容限，新增有界Anderson求根，
真实完成4个连续50ns阶段，每阶段都从已写出并严格重载的checkpoint继续。
从旧step62的650.393633ns到850.393633ns，本分支增加200ns；不是从旧末态700ns再加200ns。
外进程237.386985s，无第二次真实推进尝试。没有GPU、生产三工况、缩放或放松物理守卫。

但自由棒最大净位移仅36.174nm、肩面上下中心穿越质量都为0；仅第一阶段PID269新增塑性，
后三阶段无新增塑性。不能称持续塑性料流、沉积、定量复现或“有效75%”。

## 来源、改动和实际调用

PySPH pin仍为`1cd3dea9cfaac7d9f1bda2a2491ce4b42a31dbe5`。
事前合同：NONLINEAR_CONTINUATION_CONTRACT.md；其几何简写的更正见NONLINEAR_CONTINUATION_CORRECTION.md。
旧UPDATED_TRIAL_FREEZE的3274项在开工前全部核验通过；旧物理和失败不重写，仅更新3个交接索引。

| 新文件 | 接口 / 作用 |
|---|---|
| src/afsd_backend/pysph/stage_anderson.py | solve_stage；尺度化type-II多割线、薄QR与显式秩保护，最多6历史，仅混合未接受候选 |
| src/afsd_backend/pysph/continuation_state.py | commit/state/save/restore；接受历史累加一次、真实字段VTK、严格同源NPZ重载 |
| src/afsd_backend/pysph/transport_observer.py | observe；按持久ID独立观测lab/工具平移系位移、塑性与固定肩面中心穿越 |
| scripts/run_nonlinear_continuation.py | 有界成本probe与最多4阶段driver，每候选留档，先验收/落盘再发布接受历史 |
| scripts/assess_continuation_result.py | 独立数组读回，不调用被测接触/材料/热/acceptance/observer生成参考 |
| scripts/assess_nonlinear_continuation_cost.py | 只读几何与固定步成本尺度；首版错误结果保留且判INVALID |

真实顺序（保留旧UpdatedTrial，未重新编写物理）：

```text
accepted history
  -> candidate midpoint NNPS / corrected B / event-mean velocity gradient
  -> material from SAME S0 / ep0; pressure and conjugate internal force
  -> coupled contact + axial constraint impulses / event-resolved mean velocity
  -> conduction + plastic/friction/contact heat; five-field map/residual
  -> Anderson candidate iteration (NO history commit)
  -> original full acceptance
  -> clone commit -> transport observer -> ledger -> NPZ/VTK
  -> strict checkpoint reload -> publish accepted_history
  -> NEXT stage uses reloaded particle state
```

不存在单独在摩擦计算后偷偷覆盖自由速度的新增路径。接受算法不把混合/拟合残差当物理收敛；
必须实际evaluate最后candidate，并通过原五场、空间功、热、接触、约束和feature检查。
旧updated_trial.py、stage_iteration.py及材料/空间/接触/热源码不变。

## 根因证据，而非按函数名猜测

NONLINEAR_CONTINUATION_AUDIT.md和outputs/nonlinear_continuation_audit_20261003/assessment.json给出独立推导。

- CONFIRMED：旧32次Picard停滞，最大速度缺陷5.0427446e-4m/s在无直接接触的PID268。
- CONFIRMED：旧REJECTED候选内邻点269出现塑性更新，应变率0.99134/s接近已声明FKS低率floor=1/s；
  这是拒绝候选的分支诊断，不是原已接受态的塑性增量。
- SUPPORTED：固定几何/ep的独立方向计算中，原0.75更新跨折角至1.10054/s，Y变动约26.3%；
  这是整体率依赖反馈证据，不是调整FKS参数的理由。
- RULED_OUT（旧停止原因）：并非负摩擦、NaN、温域/密度/穿透/空间条件数或超时触发。
- POSSIBLE/NOT_RUN：全Jacobian谱、根唯一性、多分支。不能凭本轮收敛证明这些性质。

## 实际测试与独立参照

5组最终不同回归合计40项PASS；不是整库重新跑40遍，旧FAIL不删除。

| 实际入口 tests/review/ | 最终状态 | 独立参照/主要覆盖 |
|---|---:|---|
| test_stage_anderson.py | 14 PASS | g(u)=-2u+.002的根1/1500；二维根221/146、152/73；折角g(z)=1.1-9max(z-1,0)根1.01；秩亏/预算/域失败/严格JSON |
| test_transport_observer.py | 10 PASS | 手写位置/质量穿越、ID重排、被动标签、工具平移系、反向穿越、守恒拒绝 |
| test_continuation_state.py | 7 PASS | 真实27点零平衡trial接受；制造增量仅测试存储一次；NPZ所有有效属性与ID/来源拒绝 |
| test_nonlinear_continuation_runner.py | 6 PASS | 27点制造零平衡4阶段、预算边界、累计状态量不相加、发布故障不算孤立checkpoint接受 |
| test_nonlinear_feasibility.py | 3 PASS | 跨肩面初始几何的正gap=.25mm、guard触及任一侧、无上游组拒绝负时间 |

制造四阶段是driver接口测试，不能当真实AFSD。求根器首轮2 ERROR记录保留；两轮修复分别为
回退不能丢最佳锚点/fixture形状，以及非有限QR诊断使用null+状态而非非标准JSON数。
真实运行前另经审查修正driver状态量累计、异常evaluate计数与发布事务，相关回归首跑通过。
未改变数值物理容限追绿。

## 唯一真实运行及独立读回

probe完整evaluate含压缩证据4.044927s，按事前公式选总预算69次，保守预计448.649913s。
实际61次；内540s/外570s保护均未触发。峰值工作集135417856B（约129.14MiB）。
配置SHA：`f8cb21fea743323fcbd65dbe69c17e958a23eb8a1db9d3e463157751682050ae`。

| 阶段 | 终时/ns | evaluate | 原容限归一根残差最大值 | 新塑性功/J |
|---|---:|---:|---:|---:|
| 1 | 700.393633 | 16 | .864020 | 6.4845985868e-9 |
| 2 | 750.393633 | 15 | .151771 | 0 |
| 3 | 800.393633 | 15 | .581994 | 0 |
| 4 | 850.393633 | 15 | .441446 | 0 |

归一残差<=1只表示达原求根精度，不是论文误差比例。
独立读回首跑PASS：43源/97产物、4接受态、23411实际接触子区间，外2.944432s。
核对NPZ/VTK/ID/质量/工具位姿、动量、真实带符号摩擦Q、热与空间功、累计账和穿越。
最大子区间动量差1.32349e-23N·s；最大mcΔT热账差1.32649e-15J；总raw摩擦Q=2.88148e-8J。
约1e-26J负舍入原样保留，未取abs或截为0。
Es仍非认证自由能：四阶段总能量诊断残差约[-4.51414e-9,2.43253e-12,-3.58889e-12,-7.40633e-12]J，
不补成热，不因热账PASS掩盖旧闭合弹性循环问题。

## 为什么不能继续靠相同步长堆到料流

有效只读报告：outputs/nonlinear_continuation_20261003/process_feasibility_corrected.json。
首版process_feasibility.json错误用全棒有符号min，产生负时间/负成本；**逻辑FAIL/INVALID**，
原文件/日志保留，不依赖它判断物理解。纠正记录与3条回归均落盘。

1. 原始方棒150点，75肩上/75肩下；后者是预置初态，不是本轮沉积。
2. 全部棒中心距固定肩面最近250um，冻结最大总位移只有50um。
   这意味着在当前debug守卫内，两侧都不可能出现中心穿越；只重复短步无法改变这个事实。
3. 肩上最近距离250um / 轴向规定速度2.12mm/s = 0.117925s，仅是执行器长度时间尺度，
   不是自由材料真实穿越时间预测。按当前50ns步、59.1181s/步粗线性外推约2358491步、1614天。
   该4.42年数字不是已验证长时成本；它只说明当前路线不能直接堆步成为过程求解器。

## 下一真正有效单元和禁止捷径

不再以增加50ns短步数量作为主要进度。下一单元必须同时解决**可用过程时间步与允许材料通过的物理域**：

1. 保留本轮为同史五场耦合基线。在updated_trial/stage_anderson/driver旁单独设计宏步准确性与成本合同；
   先联合惯性/滑移反向/热/客观率基准，再允许一次预估<=600s的同终时比较。
   不能因隐式根收敛直接放大dt，不能从当前预压瞬态直接去惯性，不能用旧制造Q7当真实料流。
2. geometry.py、rigid_contact_kinematics.py、coupled_boundary.py需要有限孔口/肩缘feature及进出接触能量测试；
   位移安全域只能在覆盖所经feature后重新声明，不能直接删除50um守卫。
3. particles.py/output.py及observer需要持续入口协议：新ID、质量/动量/焓外流入账、原始标签与checkpoint一致；
   当前observer明确拒绝插粒，有限150点棒不是已实现持续入口。
4. 新同变体外部restart入口必须核公开接受历史；再运行不间断/重新求解同终时对照。
   内部每步实际重载已经运行，但独立两条动力轨迹等价NOT_RUN。

宏步/时间步与空间敏感性、持续入口、有限棱边过程、持续塑化/沉积、三工况和论文定量验证均NOT_RUN或BLOCKED。
材料来源/低率正则、热启动、预结合和Jaumann自由能风险沿用旧近似标签，不在本轮偷偷替换。
本轮达到：**代码实现 + 短回归PASS + 真实4步耦合接受 + 独立读回**；没有定性AFSD图形对照或定量论文验证。
命令见logs/nonlinear_continuation_20261003/COMMANDS.md；新FREEZE只证明文件完整性，不授予生产资格。
