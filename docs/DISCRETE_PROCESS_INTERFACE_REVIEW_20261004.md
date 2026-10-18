# DG 准静态变体的提交与能量接口审查

2026-10-04；只读代码和已保存输出审查，没有修改冻结源，没有重新运行 AFSD。
审查对象：`discrete_process.py`、`discrete_contact.py`、
`interval_acceptance.py`、`run_discrete_contact_process.py` 及其调用的原QS接口。
范围为 v1 源和 v1 macro 实际输出；不是全项目物理验证。

## 已确认的实际控制

1. `DiscreteContactTrial` 只接收新分类的独立历史；累计能量记录逐前缀验证，
   并要求记录数等于物理接受步数，旧模型不能通过分类适配自行恢复为新模型。
2. `accepted_copy` 从提供的位移与 dt 完整重新计算 `fresh`，随后用同一个
   `fresh` 计算原始能量并执行 `require_interval_energy`，之后才调用原提交函数。
   因而累计门槛确实进入提交路径，不是仅输出到报告。
3. 原提交函数另检扫掠、自由力、当步原始能量、穿透、质量，复制粒子后才更新
   材料与切向簧历史；新穿孔观察随后处理初始肩上棒粒子与当时实际旋转方孔交点。
   这些观察字段没有进入力学或热学。
4. `energy_diagnostics` 使用新接触模块实际返回的 `wall_work_J` 和当前真实
   力反力功，储能由同一候选材料/接触状态计算。兼容字段
   `normal_endpoint_work_error_J` 在新模块内是 DG 力的恒等式误差；旧端点力
   的假想误差另存 `normal_endpoint_quadrature_error_J`。
   未见用后者抵扣全局缺口或把材料 BE / 切向算法损失转成热。
5. 热过程仍严格拒绝温域外、负物理热源和热代数不平衡，未将 DG 变体改成另一套
   FKS 参数或材料热源。原材料 BE、密度 BE 和有限应变 hypoelastic 限制仍存在。
6. 先输出不可恢复的 `UNPUBLISHED_candidate.npz` 并读回粒子字段和原状态关键量；
   源和参数冻结核对完成后，`accepted.npz` 使用临时文件原子替换发布。

## 真实 v1 macro 结果

`outputs/discrete_contact_process_20261004/v1/macro/manifest.json`：

```
status = REJECTED_NO_ACCEPTED_STATE
error = ArithmeticError('tangential spring work identity failed')
accepted_steps = 0
accepted_time_s = 0
wall_seconds = 10.628151900000375
```

异常发生于 `newton_direction -> trial.tangent -> contact.evaluate` 的中央差分
接触微扰，不是能量筛选失败或已收敛联合根。最后完整迭代日志为 iteration=6、
完整求值计数15、自由力残差0.7464189018213343 N、门槛1.0166916460658443e-5 N。
由于候选状态未保存，本次不能独立重算其材料/接触能量分解。

## 应修复的证据与发布边界

**异常候选遗漏。** `solve` 只在收敛、预算退出或回溯失败时写完整 journal；
`newton_direction` 的组件异常会直接抛到 worker。v1 macro 实际只有
`step_001/iterations.json`，没有最后完整候选 NPZ/JSON。
下一版本应在切线异常分支保存当前完整 `r` 后再抛出，并记录扰动轴、符号及
对应接触失败信息。保存的是拒绝/未提交证据，不能改变接受时间。

**新增状态尚未读回逐项核对。** 复用的 `verify_snapshot` 核对全部粒子属性，
但状态只核对 time、step、classification、crossing_ids、contact_history、free；
没有核对新 `interval_energy`、`bore_observation`、energy 或来源哈希。
下一次生产式续跑前，应在DG发布适配器读回并验证这些新增状态，不应仅凭
文件可打开推定累计门槛记录和穿孔证据完整。

**发布后异常的事实恢复。** `accepted.npz` 原子发布之后到 worker
`history.append` 之间有短窗口。若此处异常/硬停，磁盘可能有已验收并发布的
checkpoint，而内存 history 与 manifest 显示零接受。该窗口不会把坏状态验收，
但会遗漏已发布事实。异常报告应从已发布且经哈希/状态验证的文件恢复接受索引；
不要仅依赖内存历史，也不要把不可恢复的 UNPUBLISHED 文件计入。

上述修复属于输出证据与异常边界。它们不降低原始能量/物理门槛，也不改变
v1 已记录的失败结论。v1 无接受checkpoint，因此本次窗口问题没有造成实际
接受状态歧义。
