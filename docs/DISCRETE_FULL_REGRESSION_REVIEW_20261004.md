# 当前完整 unittest 回归的失败审计

本轮完整发现并执行的 Python unittest **没有全通过**：1237 项，11 FAIL、1 ERROR、
3 expected failures；测试体 83.364 s，外部 84.789763 s，退出码 1。
原始证据为 `logs/discrete_contact_process_20261004/full_regression/stdout.log`
及同目录 `receipt.json`。这是组件与证据检查，不是 1237 个 AFSD 物理验证。
`setUpClass` 出错还会阻止该类方法执行，不能由计数推称所有方法均已完成。

本审计逐项读取原日志、旧日志、冻结文件及当前源哈希，保留所有原失败。
另经授权执行了 3 个有界组件回放子进程，核对 Q7 线程环境与旧 refined 失败；
没有修改算法、测试门槛或物理输出，没有推进 AFSD 时间。

## 逐项分类

以下表格的数值均取本次完整日志；“历史”表示找到对应原始失败证据，不等于风险消失。

| 当前非通过项 | 当前原因/量值 | 分类、历史证据与当前影响 |
| --- | --- | --- |
| `test_coupling_benchmarks.LegacyRed.test_original_friction_must_resolve_true_reversal` | 旧摩擦反转终速 −2，而独立参考 −2/3 | 历史刻意保留的旧算法红灯。源文件明确 `--legacy-red intentionally FAILS`；`logs/contact454_20260930/final_suite/legacy_reversal_RED.log` 已失败。新 DG 没有调用该旧反转函数，但不能删除反例。 |
| `test_mechanics_review.OpenMethodRisks.test_rigid_rotation_sph_boundary_must_not_create_strain` | 未修正梯度的虚假应变相对量 0.4208896 > 1e−11 | 历史开放方法红灯；`logs/contact454_20260930/final_suite/original_scientific.log` 已失败。该测试直接使用未修正 `operators`。新过程随后重建修正矩阵，故此项不是新路径使用修正梯度失败的直接证据。 |
| `test_mechanics_review.OpenMethodRisks.test_unequal_density_table1_operators_must_be_power_conjugate` | 相对功缺口 0.0976246 > 1e−12 | 同一旧日志已失败；测试未传 `correction` 的旧空间配对。新过程使用修正后的功共轭配对并独立检查虚功，不能以新检查通过改判旧失败。 |
| `test_mechanics_review.UnresolvedThermodynamics.test_prestressed_compression_needs_consistent_free_energy` | 相对自由能非共轭量 0.9999985 > 1e−10 | 历史、仍真实存在的材料风险；同一 2026-09-30 日志已失败。常 G Jaumann 更新与随体积变化的储能定义不一致，减小 dt 不能消除该瞬时不相容。新 DG 仍复用此材料更新，法向接触改进没有解决它。 |
| `test_spatial_probes.BoundaryRed.test_pointwise_constant_stress_boundary_patch` | 点态边界残差 0.3261572 > 1e−10 | 历史真实空间边界风险；`logs/boundary_step1_20260930/original_boundary_RED.log` 已失败。该测试已经使用修正算子，不能归因于“仅旧未修正梯度”；弱形式功一致不保证自由表面逐点牵引准确。对当前过程的边界精度仍有关联。 |
| `test_upstream_solid_limits.UpstreamSolidLimits.test_pure_xz_shear_must_return_to_three_dimensional_yield` | 上游纯 xz 剪切限幅后仍超屈服，相对误差约 1 | 上游采用守卫，源码明确直接调用固定 PySPH 方程并保留 FAIL；`logs/engineering_e1_20261001/upstream_red.log` 已失败。当前项目自己的三维 J2 更新不等于此上游限幅器。 |
| `test_upstream_solid_limits.UpstreamSolidLimits.test_stationary_pressure_pair_must_have_zero_power` | 静止正压粒子对上游能量项为 2，参考为 0 | 同一上游旧日志已失败；禁止将该上游式直接当成已验证的机械热源。当前模型使用独立塑性/摩擦功及导热账。 |
| `test_event_response.ActualEventResponse.test_two_free_axis_cross_coupling_independent_ode` | B[1,1] 为 0.0755199944032，独立 ODE 参考 0.0755606707296；差 4.06763264e−5 > 2e−5 | 既有连续参考精度失败。`logs/process_event_20261003/event_response_initial.log` 已失败，不是有意制造的通过项；旧事件响应的可用性不能替代 ODE 精度。新 DG 不直接调用此事件响应求导。 |
| `test_event_response_refined.RefinedResponse.test_real_pid1_old_unavailable_refined_available_and_pure` | `REFINEMENT_BUDGET_EXHAUSTED`，两次细分仍不可用 | 既有算法可用性失败。`logs/process_event_20261003/event_response_refined_initial.log`、`preflight_refined/event_response_refined.log` 均在相同断言失败；2026-10-03 preflight 总结也明确 FAIL。源/测试与 `docs/PROCESS_EVENT_FREEZE.json` 一致。本次单线程定向回放仍失败，不能说线程修复或 DG 已解决它。 |
| `test_process_corner_q4.CornerQ4Tests.test_12_j2_staggered_shared_state_partial_release_and_free_residual_shape` | 卸载后 ep 与卸载前 ep 最大差 4.91845974e−12 > 1e−12，8/27 点超限 | 既有严格数值回归失败。2026-10-01 `q4_tests_initial.log`、`q4_tests_final.log` 和后续 `process_friction_q5_20261001/q4_compatibility.log` 均已有同一 ep 断言、同一最大量级失败。Q4 源、测试、底层塑性源与 Q4 冻结哈希一致。不可仅因数值小而改判 PASS。 |
| `test_process_plastic_friction_q7_recovery.BoundedPlasticFrictionRecovery.test_02_original_twelve_correction_failure_replays_exactly` | 旧 12 次修正失败轨迹逐位比较不等；首项条件数 27.794800671844232 对 27.79480067184424 | 本轮原环境的数值再现失败；不能先验归为历史红灯。2026-10-02 多份日志该项曾 PASS。定向控制回放现已确定：相同用例/命令/源/输入，未设置三线程变量时 FAIL，按旧三项 = 1 时 PASS。旧完整日志中的 FAIL 仍保留，后者只证明旧环境下这一组件的逐位重现。 |
| `test_engineering_integration_e3.IntegrationEvidence.setUpClass`（ERROR） | `changed source: src\\afsd_backend\\pysph\\contact_cone.py` | 旧证据绑定当前源码的哈希门槛失效，不是求解器本轮运行崩溃。E3 档案要求旧 hash `4993fd7d…`，当前为 `d9d877e4…`；后者已经存在 2026-10-01 Q4/Q7 冻结。本次逐项比对 E3 source inventory，仅此一项不一致。档案未改动，但该 reader 不能据此为当前代码签发原证据检查通过。 |

11 个 FAIL 因而可分为：7 个有原始历史证据的旧方法/上游红灯、3 个有历史证据的
算法/精度失败（ODE、refined、Q4），以及 1 个在本轮原环境出现、随后经受控回放定位
到线程运行设置的逐位重现失败（Q7）。这些类别不是豁免或自动接受清单。

## 三个重点用例何时已变化

- **event_response_refined**：最晚在 2026-10-03 初始实现测试已经失败；初始和 preflight
  原始日志均是 `REFINEMENT_BUDGET_EXHAUSTED`。当前源/测试哈希与当日冻结相同。
  本轮在旧单线程环境仍失败，未找到“此前通过、本轮 DG 导致失败”的证据。
- **process_corner_q4**：最晚在 2026-10-01 的 initial/final 回归已经失败。
  当前最大 ep 变化与旧 final 相同量级且打印值相同，相关冻结字节相同。
- **q7_recovery**：2026-10-02 initial、final、后续 feasible predictor/onset 兼容日志均有 PASS。
  当前全量回归首次在所查日志中发现 FAIL；不能靠“旧模块”这个名称归为原有失败。
  旧 16 项物理源/配置、诊断恢复函数和测试字节均未变。旧命令记录明确三线程变量均为 1，
  本次完整 runner 没有显式固定这些变量。后述控制回放确认环境设置足以切换该项结果。

## 有界线程环境控制回放

新证据保存在 `logs/discrete_contact_process_20261004/full_regression_replay/`，
`receipt.json` 包含完整命令、线程环境、子进程退出码、外部耗时及源/输入哈希。
每个子进程外部上限 60 s。Python 3.11.16、NumPy 2.4.6；没有安装依赖。

| 子进程与日志 | 环境 | 实际结果 | 外部耗时 |
| --- | --- | --- | ---: |
| `q7_current_environment.log`，仅 Q7 exact 用例 | OPENBLAS/MKL/OMP 三项均 unset | 1 test，1 FAIL，退出 1；首差与完整回归相同 | 0.649496 s |
| `q7_single_thread_control.log`，完全相同 Q7 exact 命令 | `OPENBLAS_NUM_THREADS=MKL_NUM_THREADS=OMP_NUM_THREADS=1` | 1 test PASS，退出 0 | 0.604672 s |
| `q7_and_refined_single_thread.log`，Q7 exact + refined PID1 | 三项均 1 | Q7 PASS、refined FAIL；2 tests、1 FAIL，退出 1 | 2.011124 s |

前两者的相同命令对照支持本次 Q7 exact 漂移由线程运行设置引入；没有定位到底层哪一个
BLAS/LAPACK 调用，也不声称任意平台都能逐位重现。该控制没有放宽断言、覆盖旧失败或将
原本拒绝的 Q7 候选接受。refined 仍失败表明其预算问题不能用此环境解释。
显式固定单线程是数值再现约束，不是物理质量、时间或模量缩放。

## 对新 DG 联合路线的影响

本次完整日志中 `test_discrete_contact` 17 项和 `test_discrete_process` 3 项均显示 PASS；
这里没有发现这些新组件的直接失败。它们覆盖本身声明的范围，不能证明完整过程精度，
也不能覆盖旧 ODE、Jaumann 自由能、修正 SPH 点态边界以及 Q4 卸载问题。
仅调用图未显示 Q7/refined 直接依赖新 DG；这不等于全局环境与测试顺序永远无影响。

原完整回归的最终状态仍是 **FAIL**。后续若在固定环境重新执行完整回归，应另存新日志；
不得用单项 PASS 替换本次完整结果，也不得把历史红灯改成 expectedFailure 来制造全绿。
真实短程试跑仍需独立模型合同、区间能量、接触可接受性、真实孔口事件及成本门槛；
任何局部通过都不能授权“已完成 AFSD 沉积过程”的结论。

## v3 冻结后固定单线程完整回归

本节是独立的新执行记录，不替换上述最初 1237 项结果。v3 物理实现冻结后，通过
监督子进程执行 `scripts/verify_discrete_project.py --worker`，在启动 Python 前显式设置
`OPENBLAS_NUM_THREADS=MKL_NUM_THREADS=OMP_NUM_THREADS=1`，外部硬上限 540 s。

结果：**实际 1246 项、10 FAIL、1 ERROR、3 expected failures，整体仍 FAIL**。
测试体 78.871 s，外部 79.824998 s，子进程退出码 1；没有 AFSD 过程推进。
Q7 exact 用例现在 PASS，其余最初 10 个 FAIL 及 E3 `setUpClass` ERROR 均再次出现，
没有发现新增失败名称。不能把减少 1 个失败描述成解决了一项物理模型问题。

开始前 loader 列出 1251 个方法；E3 `setUpClass` 出错使该类 5 个方法未执行，故
实际执行 1246 个，与 unittest 计数一致。没有通过项或未执行项的模糊合并。
此时 `test_secant_equilibrium.py` 尚未出现在发现清单，未纳入本次结论。

独立新证据目录为 `logs/discrete_contact_process_20261004/full_regression_v3_single_thread/`：

- `stdout.log`：完整原始输出，保留所有失败及 traceback。
- `receipt.json`：命令、环境、退出码、耗时、结果和源核验摘要。
- `source_test_snapshot.json`：启动前 `src/`、`tests/`、`scripts/` 全部现有 `.py` 文件哈希。
  结束核对，既有文件变动 **0**，运行期间新增 `.py` 文件 **0**。
- `planned_tests.json`、`executed_test_ids.json`、`not_executed_test_ids.json`：计划、实际
  执行及受 class setup 错误阻止的方法身份；原始 test 行另存供核对。

这次完整执行确认固定环境下当前回归的实际状态；未撤销历史物理风险、旧证据哈希
限制或 refined 可用性失败，也未认证 5 ms 试跑之后的孔口输运与沉积。
