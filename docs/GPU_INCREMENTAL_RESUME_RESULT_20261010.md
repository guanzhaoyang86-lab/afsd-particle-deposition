# 结束授权超时修复与续跑（2026-10-10）

## 14:49最新实际进展

恢复任务已完成首次历史核验（533.070秒），旧第56段重演三个终态粒子数组全部逐位相同。现有6段完整PASS至56.3076171875ms/380细步，名义dt已恢复474.609375µs；第7段继续运行，监督56536/子50836存活，原截止16:44:18不变。最新完整manifest和检查点SHA独立核对一致，粗细比较及来源检查通过。

首次新料穿孔出现在55.3583984375ms。最新4个初始eligible新料ID为1151、1203、1229、1281，合计1.35e-6kg（1.35mg）；肩下新料仍0。穿出新料半径约0.508–0.520mm，尚未达到肩下环区；这证明发生了新料穿孔，不能证明形成了沉积层。后续观察横向铺展和持续留存。最近未穿孔ID已换为1333，距出口64.759µm；该数值换了观测对象，不是原已穿孔前沿倒退。

以下为启动时记录。

用户要求“继续启动”。旧任务第57段的三个候选资格和粗细比较均通过，但worker在结束授权处等待45秒超时，控制器随后发送回复时出现BrokenPipe；该段不作为完整父态。最后完整第56段为55.00244140625ms/370细步，新料距出口2.897133µm、穿孔/肩下均0。

## 修改与验证

新增`afsd_incremental_worker_resume_20261010.py`将每次授权等待上限设为180秒，并按原本段525秒截止截断，发送/等待/接收共用同一截止，接收后再次确认未超期。无真实有效授权仍停止。原物理求解、540秒hard timer、600秒外层限制和450秒成本预测不变。

新入口`run_gpu_incremental_resume_20261010.py`仅更换明确版本标识、合同、组件收据及worker导入，整个控制器AST经过等价核对。旧源及失败证据未改；仍用冻结的物理source_sha256。新行政来源单独绑定，历史锁、增量核验、步长恢复和肩下观察保持原逻辑。

`scripts/check_gpu_incremental_resume_20261010.py`实际执行16项测试全部通过，返回码0。覆盖模拟60秒回复、180秒上限、525秒剩余预算、错误身份/阶段/锁状态、断开/过期回复拒绝、实际Windows spawn/Pipe开始与结束双授权、整个物理阶段仅增加两处deadline参数、完整第55→56段父态及SHA绑定。旧组件13来源与日志字节SHA重新核对不变。

组件证据：`outputs/gpu_acceleration_20261007/incremental_resume_components_20261010_v1.json`和同名`.log`；相关新源、测试和合同均已封存。

## 已执行后台启动

2026-10-10北京时间14:31:32.787820启动，监督PID56536，实际截止16:44:18.787820；使用原计划剩余7966秒，未新增四小时。初始实时状态`INITIAL_FULL_HISTORY_VERIFICATION`，进程存活且CPU时间增长。bootstrap为`PASS_DETACHED_JOB_POLICY`、`kill_on_job_close=false`。首次完整历史核验后，从旧第55段重演旧第56段，要求三个终态粒子数组逐位一致，然后从第56段继续；不会直接使用未完整的第57段。

实时状态：`outputs/gpu_acceleration_20261007/incremental_deposition_resume_20261010_v1/progress.json`；新会话最初仍需一次全历史核验，后续才复用本会话持锁的历史内容，不把核验时间计作模拟进展。首个完整新结果尚待生成。

实际启动命令（已启动，不要重复执行）：

```powershell
& '.envs/pysph/python.exe' -B scripts/run_gpu_incremental_resume_20261010.py --detach --output outputs/gpu_acceleration_20261007/incremental_deposition_resume_20261010_v1 --parent outputs/gpu_acceleration_20261007/incremental_deposition_20261010_v1/attempt_0055 --replay-stage outputs/gpu_acceleration_20261007/incremental_deposition_20261010_v1/attempt_0056 --wall-seconds 7966 --maximum-dt 0.000474609375 --target-time 0.22 --shoulder-observe-ms 5
```

最后完整候选最长耗时10.004757秒，沿用预估30+4.5×最长候选=75.021407秒，低于450秒。1381粒子，dx=0.5mm，恢复dt=118.652344µs，仍为工程近似诊断；全部原物理检查和首肩下后至少5ms/5段观察限制保留。实际新穿孔与沉积须看新完整结果，不预先宣称成功。
