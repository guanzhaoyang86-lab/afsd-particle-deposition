# 新独立沉积模型的受限分段诊断推进

用户已明确授权继续真实求解沉积。本入口只使用
`afsd_deposition_setup_20261007` 的新独立矩形域/较长棒料，不接旧342粒子结果，
不迁移旧3ms诊断父态，不改冻结材料/接触/热/空间算子，不给自由材料预设运动。

`run_afsd_deposition_stage_20261007.py --stage cost` 创建独立初态，实际测残差、
切线、线性求解、标量/批量等价及峰值内存。未求力根的成本映射不能证明稳定步长。
先要求真实组件/算子收据，包括setup、控制器、沉积证据观察测试。

`--stage stage --cost ... --end-time ...` 从该cost初态开始；`--parent ...`
只能指定同一新模型、同一冻结实现已通过全部宏步粗细比较的阶段目录。
检查父阶段manifest/supervision、源与全部artifact哈希、每个入选fine状态与observer、
screen哈希，并沿parent链核对独立起点和连续父态。状态仍为
`diagnostic_only=true/accepted=false/restart_allowed=false`；阶段收据仅授权同模型
的screened diagnostic continuation，不构成正式接受、实验验证或全流程收敛声明。

每一宏步以完全相同的材料/机器/能量/穿孔父态，求一个dt粗步和两个dt/2细步。
各候选执行原力、热、质量/身份、能量完整前缀、标量/批量等价、接触壳与整路径门槛；
等时比较保留原全局、全部粒子局部与穿孔ID门槛。只有整对通过才选择两个fine态，
失败态不进入screened_history。数值/物理候选或时间误差失败，可从原父态减半重新求解，
受min-dt/max-halvings限制，绝不放宽阈值。明确的几何/接触壳/支持边路径失败
保存事件诊断后停止，不能自动删除边或以不断缩dt冒充事件积分。

可设置initial-dt/min-dt/max-dt；未指定initial-dt时续用父阶段建议值。
通过后只有全局及每粒子局部误差占原限值最多25%时才可按growth-factor试探增步
（1–2，默认1.5），下一个宏步仍重新求完整粗细比较，失败仍回原父态减半。
这只是局部时间误差控制，不能证明累积全流程误差已收敛。

每阶段有明确global end-time及max-macrosteps，且最多推进到当前机器阶段边界；
不得跨界使用一个候选或一个阶段。进入走刀前须有实际screened新料肩下质量，
且这些eligible crossed肩下ID本身相对于独立初态具有正累计塑性应变；不能用全体
粒子的epmax替代。该门槛仍不证明穿孔后才开始塑性或已经形成冶金结合。
不足时保存失败，不按命令结束就声称填满。`--stop-on-first-shoulder`可在第一个
通过比较的宏步观察到新料肩下质量时停止。离肩或后方几何留存不等于冶金结合。

每子进程预测≤450s、内部540s、外部600s、CPU单线程。显式max-candidate-attempts
额度包含所有真实候选调用，失败和重试同样扣额度；每个完整宏步至少3次。
预测为30+1.5×总候选额度×实测候选成本，而非忽略重试或将最坏减半数无界相乘。
运行中继续检查剩余候选额度和时间预算。
达到有限步数、阶段边界或通过宏步后的剩余预算界限，可保存通过的诊断前缀供下段
继续；预算用尽明确标记`PASS_SCREENED_PREFIX_BUDGET_BOUND_NOT_TARGET`及
stage_target_reached=false，只授权之前确实通过的前缀，失败尝试仍FAIL且不进入
screened_history，并记录下一建议dt。旧源和所有失败输出均保留。

manifest的`screened_history`仅含本段通过宏步的fine状态，每宏步恰两行，带
checkpoint/sha256/time_s、observer_path/observer_sha256、macro_index/half_index、
screen_path/screen_sha256。parent_stage_manifest及其sha256重构全链。
`screened_continuation_checkpoint/sha256`永远指末个screened fine态。
正式接受步数与时间恒为零，实际诊断物理时间独立报告。
