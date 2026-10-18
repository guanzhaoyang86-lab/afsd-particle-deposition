# 用户授权的 GPU 中断恢复

用户在获知旧任务实际中断后明确要求“继续”。本次从`campaign_10h_20261008_fast_v1/attempt_0003/fine_2`完整PASS状态恢复，物理时间25.6953125ms、156细步，SHA256为`5265b7ed9245a97198a3f482951d02ddc8b718e1bd6dc5f65d5f23541d705831`。中断的attempt_0004没有完整结果，原文件保留，新任务从父态重算该区间。

## 已执行检查

新行政恢复入口18项测试PASS，命令exit0、0.567s；监督循环与冻结fast监督器AST完全相同。检查包含原截止、准确父态、旧进程退出、来源绑定、完整外层收据，以及拒绝已有超时/失败/停止记录。未重做未改变的物理基准，也没有放宽求解检查。

检查收据：`outputs/gpu_acceleration_20261007/resume_supervisor_components_20261008_v1.json`。冻结新增文件为`scripts/resume_gpu_campaign_20261008.py`、`scripts/start_gpu_resume_detached_20261008.py`、`tests/review/test_gpu_resume_20261008.py`及恢复合同。

## 后台启动与实际状态

原工具进程的Windows Job flags为0x2800，包含允许breakaway和KILL_ON_JOB_CLOSE。两组15秒独立探针完成，breakaway子进程查询到flags0。此结果识别出原启动方式的一项风险，不能证明前次任务退出就是该原因。

新入口使用系统允许的CREATE_BREAKAWAY_FROM_JOB | DETACHED_PROCESS | CREATE_NEW_PROCESS_GROUP，隐藏窗口、重定向标准流，无提权或系统配置变更。子进程先记录Job策略；若仍有KILL_ON_JOB_CLOSE则不执行求解。实际bootstrap为PASS，PID56532、flags0、kill_on_job_close=false；in_job仍true，可能还有外层Job。

新任务目录`outputs/gpu_acceleration_20261007/campaign_10h_20261008_resume_v1`，监督于2026-10-08 20:54:49开始，原绝对截止2026-10-09 02:19:45不变。20:55读回RUNNING_GPU_CHILD，心跳更新，首段dt0.000474609375s。最近完整段耗时196.422s，最近最大候选23.694s，原下一段预测136.625s；1381粒子，原450/540/600/650s预算和失败停止规则不变。

同目录名前缀的`.launch.json`和`.bootstrap.json`记录命令、PID、源哈希和启动策略；`.stdout.log`记录后台输出。运行状态以目录内`progress.json`结合实际进程是否存活判断，不能只读RUNNING字样。退出写`final_summary.json`。

本次起点穿孔和肩下新料仍为0。所有新结果仍为diagnostic_only，accepted=false、restart_allowed=false，正式接受时间增量0。

## 首段恢复结果

20:58首段完整PASS并独立读回，物理时间26.169921875ms、158细步；监督自动开始第2段。首段总195.322s，其中预检76.249s、worker93.534s。粗步在原Broyden失败后按既有Newton备用路径于15次迭代、24次求值收敛，越过上次中断记录所在的位置；两细步及原配对检查随后通过。

新料最近中心距孔口180.8845um，穿孔及肩下仍0。最新检查点SHA256为`e397b8bbaaac76675f3de00918594e20cfb2358d67497bd4c0127b729c3457fc`。完整读回记录和绑定文件哈希见`outputs/gpu_acceleration_20261007/resume_readback_20261008_v1.json`。这证明本次恢复成功推进了一段，不能据此确定此前退出的系统原因。
