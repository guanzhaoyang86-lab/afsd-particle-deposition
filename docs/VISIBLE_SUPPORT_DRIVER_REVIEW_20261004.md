# 可见支持短窗驱动：只读事务审查

审查 `scripts/run_visible_support_process.py`，包括资格 v3 路径、RecordedTrial 和失败 Newton 最后完整图保存后的版本，SHA256 `4ce4601ab6581202c5187932dda649d0f8e9ba17d07ccd7f5c80f00165b5a94c`。没有执行 driver、预检或物理求值，没有修改驱动或冻结源码。最终路径资格 adapter 尚待落盘，不在这里冒充已具备运行资格。

唯一发现的具体发布状态缺口：`publish` 先原子写出 `accepted.npz`，调用 snapshot verification，返回后 worker 才计算 SHA、append history 并写 accepted_history。若在原子重命名与历史记账之间发生 I/O 异常或硬超时，磁盘上可能已有 `restart_allowed=True` 的合格 checkpoint，而 history 尚为空或少一步。当前普通异常终态仅按 history 是否为空选择 `FAILED_NO_NEW_ACCEPTED_STATE`，硬超时 receipt 也只保存内存 history/rows；会漏报这个已发布前缀。

这不会绕过两 half 的物理/几何/时间筛选，因为真正发布仍在共同筛选后，但会造成终态和可重启证据不一致。建议异常及硬超时清单枚举实际 `step_*/accepted.npz`，存在未记账文件时使用“已发布 checkpoint 待独立读回核验”状态，禁止断言无新接受态；由只读重放检查 screen、父 SHA 和原子 checkpoint 恢复已合格单 half。发现已即时告知 root，此报告不声称该项已修复。

其余重点核查未发现阻断性错误：

- 39.375 ms、global 431 的父 checkpoint 和原 432 帧读回摘要均固定 SHA；原参数、基础源、实现源和完整 prefix 有入口及后续源检查。cost/coarse/fine 都从该真父独立分叉，模型转接不新增物理步骤。
- cost 固定 1 步、coarse 固定 8 步；只写显式不可重启候选，不进入 fine 发布分支。fine 16 half 先按物理顺序生成两候选，只有共同全局/局部/strict-cross-ID 筛选通过后才开始发布；失败 pair 的两个候选不会混入已接受历史。
- 每 half 经过新 `visible.accepted_copy` 的 fresh 求值及原能量前缀，再通过标量接触对照、连续粒子壳和固定保留图路径证书。发布再次断言新 classification、diagnostic_only=False、已取得父等价资格、路径 CERTIFIED 及 pair screen 通过。
- 局部时间筛选累计肩接触、cross IDs 和本窗口被阻断边端点，并扩展到 coarse 与双方 half 的 raw/visible 图一跳并集。两分支是整个窗口同一真父起步，不是将 coarse 状态接到 fine；每对共同终时检查时刻相同。
- fine 真正发布的第二 half continuation_parent_sha256 指向第一 half 的实际 accepted 文件 SHA；下一 pair 继续以上一真实 accepted SHA 为父。未提交候选的数值父 hash 与发布后的元数据父 hash有明确不同用途，没有将未接受 coarse checkpoint 当作 fine 父。
- 全程保留最初 431 条 interval_energy records，原 DG 接受还验证每个后续 state 的全部记录与步数；没有每窗口清账或更改能量阈值。
- 450 s 预测门槛、525 s 内部搜索截止、540 s 工作进程硬截止和外部 600 s 均在；下一根成本判断含不断更新的最慢实测完整候选成本。Broyden 与同父同初值 Newton 回退各使用已登记的 24 更新/100 求值范围，最多一次回退；未见递归减步或隐式扩窗。

这份审查确认的是受控驱动顺序，不认证求解结果。资格文件/源码若随后变化，需要以运行冻结集合和最终只读证据为准；本报告没有授权跳过尚缺的路径 adapter 或组件终态。
