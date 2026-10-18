# 统一可见图求值、转接和求根：独立只读审查

本次只读审查 `visible_quasistatic.py`、`visibility_equilibrium.py`、`check_visible_trial.py` v2 和 `VISIBLE_SUPPORT_PROCESS_CONTRACT_20261004.md`，追读原 QS/DG 接受、能量、切线和 Broyden 辅助函数。没有执行真实物理、求根、试算或额外组件基准；没有修改被冻结实现。本结论是代码审查，不代替运行结果和过程验收。

结论：所审查代码中，统一可见图确实进入完整力热求值及冻结搜索切线，没有发现“只删除输出边”“只过滤内力”“继续用旧 B/旧热矩阵”或跨图进行 Broyden 秩一更新的问题。显式转接及所用原接受链也保留能量和被动观察历史。下列外围约束和未覆盖项须与这一结论同时保留。

## 已确认的接口

- `visible_quasistatic.py:83–91`：实际端点 PySPH NNPS 后立即过滤，随后密度固定点每轮 `_geometry` 都使用同一 i/j/g；H、B、rho、L 和材料率从该图产生。原图另存 raw_i/raw_j/raw_g，没有进入后续力热组装。
- `visible_quasistatic.py:116–131`：隐式导热使用同一子集的原径向 g，内力使用同一 i/j/g 及该图的 B。材料应力和塑性热均来自该 L；没有通过修正机械梯度替换热学径向梯度。
- `visible_quasistatic.py:157–162` 与原 `quasistatic_process.py:160`：切线先核其 i/j/g 与 result.visibility_graph 完全一致，再由原 gradient_map 组装 D。体积、材料状态及接触也来自该完整候选。原搜索切线仍明确冻结 D、V、几何和热反馈导数，不是全耦合精确 Jacobian；这项原有限制没有被包装成新证明。
- `explicit_transition` deep copy 原 state，新增明确模型转接元数据并保留 time、step、能量前缀、接触历史、crossing 和 observer。粒子数组没有被转接函数改写；evaluate 只修改 extract 得到的副本。类初始化仅在调用旧构造器时使用临时旧分类，随后 self.state 恢复新分类。
- 新 `accepted_copy` 经原 DG 接受函数先对相同位移重新完整求值，再用 `require_interval_energy` 在原 records 后追加有符号 raw 与绝对工作尺度，逐步及全部前缀检查；原 QS 接受还检查力、粒子接触 sweep、罚壳和单步原始能量，并写实际材料/温度字段。新 state 由 trial.state 派生，故不会退回旧模型分类或清空 431 步历史。没有扣除算法耗散或将余项注入热源。
- `visibility_equilibrium.py:141–149`：实际接受的搜索点拓扑 fingerprint 变化时不调用 `_broyden`，设置强制重建原因。若尚需下一方向，下一轮先新建完整冻结切线；若两个刷新预算已尽则拒绝。跨图候选恰已通过完整力收敛时可以直接返回，因为没有再使用旧矩阵生成下一方向。拒绝的回溯候选不改变当前搜索状态。24 更新、100 完整求值、2 刷新以及真实力而非近似力的判断仍在。

fingerprint 是规范化持久 PID 的无向边拓扑，不包含 g 的浮点变化。对于冻结模型内的 Broyden，这允许同拓扑下几何、B、温度、材料随完整位移变化，并只在离散支持拓扑改变时禁止 secant；这与登记范围一致，不要求每次 g 更新都重建。

## 具体缺口与使用边界

1. **禁提交是新包装入口保证，旧通用接口并不会自动识别 diagnostic_only。** `visible_quasistatic.py:50–53` 在调用 core 前明确拒绝 diagnostic-only，当前 harness 也测试了该入口；但直接调用 `core.accepted_copy` 不包含这一标志检查。相同地，完整 retained-edge 路径证书并不在此 trial 或 accepted_copy 中运行，必须由新 driver 在发布前强制核验。root 已明确将 driver 全程限定为 `visible.accepted_copy`，且发布前再断言新 classification、diagnostic_only=False 及资格。未审查完成的 driver 不能仅凭此模块报告视作可发布。

2. **“所有字段等价”大于现 harness 的实际断言范围。** `check_visible_trial.py:23–31` 加强了 i/j/g、B、volume、rate、internal_force、材料反馈温度、热导/配对热及热总账，但调用的 `batch.candidate_equivalence` 对 contact 只逐项比较 ndarray，忽略字符串、标量和嵌套 sweep 字典；机械能量比较也只有 plastic_work、global_raw 和 absolute_work_scale 三个汇总项。没有逐字段比较所有 material tuple、particle_energy_J 或所有能量分项。源码未显示这些项发生变化，但当前 PASS 不应写成每个嵌套字段均已断言相等。可通过独立补充审计或缩窄声明处理，不能把聚合相等当作一般逐粒子能量相等证明。

3. **harness 失败状态与直接 helper 来源记录较弱。** v2 main 无捕获后写 FAIL 的外层处理，异常后 manifest 会保留 RUNNING_DIAGNOSTIC_ONLY；新 driver 必须只接受精确 PASS 终态，不能以目录存在或父等价中间字段准入。运行前冻结集合记录 core 和列出的新文件，但没有把直接调用的 `run_batched_process.py` / `run_discrete_post_start.py` 加入本次结束时复核集合；原父 implementation_source 在入口检查一次包含相应历史来源，并不等价于本次全过程直接 helper 均被同一冻结集合覆盖。建议外层准入固定和复核这两个 helper 及资格产物 SHA，不必修改旧证据。

4. **契约要求的条件数记录未完整覆盖改变图的真实诊断。** 旧失败候选的新图诊断保存 B、余项和能量，但 `old_failure_map_diagnostic` 没有显式 moment/B 条件数摘要；原 `_geometry` 仍执行 1e6 条件数拒绝，故不是绕过退化门槛，只是诊断记录不全。可从已保存 B 只读计算，而无需重新求值。

5. **搜索失败 journal 不含完整 graph 证据。** visibility controller 复用旧 journal，仅保存候选位置、材料、力、接触和汇总能量；iterations 中有 topology fingerprint，但 journal NPZ 没有 raw/kept i/j/g/B 或 blocked 原索引。完整候选的 graph 归档需要外层 driver 负责；尤其求根中途失败的最后候选，不能宣称旧 journal 本身已保存全部新图证据。可只读按位置和已冻结参数重建几何，但这不是“当时实际图已完整落盘”的同义词。

## v1/v2 和物理解释

原 v1 harness 使用公开 restore 读取 NOT_RESTARTABLE 参考而被拒绝，是原防护有效。v2 改为仅诊断的显式 load 并断言 restart_allowed=False，随后 explicit_transition(..., diagnostic_only=True)，并验证新包装接受入口拒绝。这个例外没有将旧候选变为合格父，也没有增加接受时间。v1 失败目录须保留，不能用 v2 PASS 覆盖旧状态。

删除边改变实际空间算子。固定端点可见图把终点已阻断的边排除整个 BE 增量；retained-edge 连续路径证书只保证保留边沿登记路径可见，不能证明被删边在区间早期的物理作用已得到时间积分。拓扑切换可以使受力变化不光滑，完整残差收敛、原始功热账和同终时局部步长比较仍是必要检查。它们也不能代替空间收敛：独立制造二次场在 250 粒子分辨率的导数 RMS 仍约为解析连续域参考尺度的 45.2%。

审查版本 SHA256：

| 文件 | SHA256 |
| --- | --- |
| scripts/visible_quasistatic.py | bae4f5b12e751da53ce426987b1d49e098529260469dea934c5db5bf3c374540 |
| scripts/visibility_equilibrium.py | 1371979ffba97a1dd41c99fbfccd8f73108ab3ac9fcdc546b3b835fe937d5949 |
| scripts/check_visible_trial.py（v2） | 46a73eb31579a7b9fae12d5b40af028fe5b8c576544619e68dbe9fda7e9e39a5 |
| docs/VISIBLE_SUPPORT_PROCESS_CONTRACT_20261004.md | 85ab1c995cb19ff04fd7aff6123e810b63d546352816369607672230f28c87ae |
