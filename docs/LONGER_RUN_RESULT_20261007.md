# 延长运行与新增孔口穿越（2026-10-07）

**最新：后台已于北京时间 19:04:09 自动停止；最后接受状态为
62.1875 ms / 1143 步，8 个合格粒子穿孔 / 2.7 mg，未达到 100.1875 ms。**
停止原因及最新观察见文末；下文运行中记录均为此前快照。

用户要求延长实际计算，让更多粒子穿过并观察沉积。沿固定 PySPH 的
342 粒子热态工程近似分支继续，物理参数、自由材料受力、质量与能量账
均未改变。100.9375 ms 是本轮最初提出的延长目标，**尚未达到**。

## 原步长续跑及停止

从 50.9375 ms / 687 步接续 `tight_path_process_v2_20261004/window_002`：
粗支 16×62.5 us 全部通过，退出 0，监督耗时 51.2082 s。
细支新增 8×31.25 us 至 51.1875 ms / 695 步后，step009 几何路径检查失败，
退出 2，监督耗时 37.6699 s。独立读回退出 0 / 24.8318 s，明确
`window_completed=false`，仅确认八步有效前缀；旧失败候选未被晋升。

失败为 `MAXIMUM_ANALYTIC_BATCH_QUERIES`：256 查询上限用尽，深度 8/32。
7797 条支持边中，7796 条已有完整路径证书；唯一剩余边 PID195–209
尾部 1.953125% 时间段未完成认证。未发现严格工具交叉，不代表未检查
区间已经安全。只读取证见 `outputs/longer_run_diagnostic_20261007/v1/summary.json`。

## 保持物理与判据的有界步长细化

新协议 `TIGHT_PATH_REFINEMENT_CONTRACT_20261007.md` 在运行前登记，
从上述已读回的真实 695 步父态重新求解；粗 15.625 us、细 7.8125 us，
每窗最多增加 0.25 ms、最多四窗，完整读回通过才允许下一窗。
256/32 路径上限、力/热/全前缀能量、标量等价、局部/全局时间比较均未放宽。

新增入口为 `scripts/run_tight_path_refinement_20261007.py` 及独立
`scripts/readback_tight_path_refinement_20261007.py`。静态适配只改变父态、
路径、来源接线及时间安排。原冻结源码/合同未修改。独立 reader 的 12 个
证据审计算法函数 AST 与原版相同；真实父态接线复核通过，26,116 次来源
条目检查、21.469 s、零物理/NNPS/几何重算。此检查及历史组件来源 SHA
核验不等于重新执行历史整包回归。

首窗粗、细及独立读回均退出 0；实际接受 32 步至 **51.4375 ms / 727 步**。
PID322 新增合格孔口穿越，累计 **6 点 / 2.025 mg**，初始基线为
5 点 / 1.6875 mg。该结果来自真实接受粒子状态，不是增加动画帧或插值。
恢复四窗最终均完整通过，共 128 个新细步；连同原窗新增八步，本轮已确认
新增 136 步 / 1.25 ms，至 **52.1875 ms / 823 步**。全历史观察退出 0，
824 个真实帧，仍为六点 / 2.025 mg；这六点的肩下环区质量与实际肩接触质量
均为零。观察和已目视核对的结果图位于
`outputs/deposition_observer_20261007/after_refinement_004_v1/`。

## 长序列监督入口

`scripts/continue_tight_path_longer.py` 顺序调用原控制器与独立读回，
完整父窗、成本预测、原每阶段时限和失败停止均保留，另有总墙时限制、
进程树清理、`STOP_REQUESTED` 及原子进度文件。已做 AST/父态/停止路径
检查及 Windows 父/子进程清理检查。实际以 `--start-window 3 --end-window 51`
启动前检查时，正确拒绝旧 window002 的不完整前缀，退出 2，物理阶段为零；
收据 `outputs/longer_run_20261007/v1/progress.json`。该入口绑定原系列，
没有在后台运行到 100 ms，也不能跳过当前原系列失败。

## 复跑入口和结果边界

项目环境 `.envs/pysph/python.exe -B`；`PYTHONPATH` 在子进程中清空，
`OMP_NUM_THREADS`、`OPENBLAS_NUM_THREADS`、`MKL_NUM_THREADS` 均为 1。
所有原有结果目录拒绝覆盖。新恢复入口的参数为 `--window N --stage coarse`
与 `--window N --stage fine`；独立 reader 使用 `--window <window-dir>`
与 `--output <new-readback-dir>`。实际命令和退出码在每窗的 supervision JSON。

全历史观察入口 `readback_deposition_progress_refinement_20261007.py` 仅增加
新读回状态与来源绑定；保留从零时刻起的全部实际帧和 ID，不截断原塑性/接触历史。
绘图使用既有 `D:/py/anaocnda/python.exe` 和原绘图脚本，没有安装依赖。

当前基板宽 4 mm、小于肩径 5.4 mm，料棒为有限库存、没有生产级持续补料。
更多孔口穿越不等于肩下铺展或完整沉积；918 点分支和 ramp 诊断时间没有
与本分支拼接，原 FKS/Jaumann/空间及实验局限保持。

## 后台长序列已启动，目标尚未完成

四个恢复窗完整通过后，另登记 `EXTENDED_TIGHT_PATH_CONTRACT_20261007.md`，
从真实 823 步父态启动 `outputs/extended_tight_path_20261007`。恢复既有
62.5/31.25 us 粗细安排并不预认其稳定；每窗仍须完整求解、时间比较和独立
读回。最多 48 个 1 ms 窗口，目标 **100.1875 ms**，最长五小时；任何失败停止。

新物理入口 `run_extended_tight_path_20261007.py`、独立 reader
`readback_extended_tight_path_20261007.py`。reader 的 12 个核心审计函数
AST 未变，实际父态接线检查 PASS / 29,777 次来源条目 / 28.3107 s。
收据 `outputs/extended_tight_path_20261007/readback_wiring_preflight/summary.json`。
后台监督 `continue_extended_tight_path.py` 的 11 项纯控制流/守卫检查通过；
它们不属于物理验证。

北京时间 2026-10-07 18:21:22 已启动后台监督 PID **52064**，确认首窗粗支
正在实际推进。启动记录 `logs/extended_supervision_20261007/launch.json`。
此段记录是启动快照，**实时状态以**
`outputs/extended_tight_path_supervision_20261007/v1/progress.json` **为准**；
其 `latest_audited` 才是最后独立复核端点，不能把运行候选当完成结果。
停止文件为同目录 `STOP_REQUESTED`，在阶段边界生效。原候选/失败证据保持。

交付前新长序列首窗粗、细及独立读回均退出 0，完整接受到
**53.1875 ms / 855 步**；较本轮初始新增 **168 步 / 2.25 ms**。
首窗读回 SHA256 为
`78d2655b90e365a9c165fc563dd98d201eadc1bd77e28d8e632fd392cd63f081`，
末 checkpoint SHA256 为
`a5e6da84a08fa94d5688728dff8786aac9c75016106681b14cb0326d1d12a165`。
独立全历史观察退出 0，856 帧，**7 点 / 2.3625 mg**（本轮新增 PID322、214）；
已穿孔材料的肩下环区和实际肩接触质量仍为零。最新静态观察及已目视检查图：
`outputs/deposition_observer_20261007/after_extended_001_v1/`。
交付前监督仍 RUNNING、第二窗 fine 正在进行；本报告不是后续后台结果的自动更新页面。

## 19:04 停止后的诊断与最终接受状态

长序列 window001–010 均完整接受并独立读回，末态为
`outputs/extended_tight_path_20261007/window_010/fine/step_032/accepted.npz`，
SHA256 `9942cc906f54ffed95ccf4ccd880edda25d2dd0baa1a6e13c9af5943fe5db968`。
相对本轮最初 50.9375 ms / 687 步，累计增加 11.25 ms / 456 步。

window011 粗支完成 10 个未接受参考步，随后 step011（62.8125–62.875 ms）
被整步支持图几何检查拒绝。`visible_path_audit.json` 状态为
`SOLID_OCCLUSION_FOUND_ON_RETAINED_EDGE_PATH`，failure_kind 为
`STRICT_SOLID_INTERSECTION_ON_FIXED_RETAINED_PATH`：PID195–209 连线在
fraction=0 与工具实体相交长度 7.34094944e-8 m。只使用 2/256 次查询，
没有未解区间或数值域失败；本次不是此前 51.1875 ms 的查询预算不足。
这说明候选采用的支持连线不能覆盖完整步内几何；不等于粒子中心穿进工具。
保存图对照进一步确认：step010 的终点图把同一 PID195–209 边列为 blocked，
step011 终点图则将它保留。当前候选用终点支持图检查完整路径，因该边在
步初仍被遮挡而拒绝；这是遮挡转可见的图切换问题。独立连续接触壳检查
仍为 CERTIFIED_NO_HIDDEN_SHELL_CONTACT。不能保证只减小步长即可解决。

监督原样停止，fine=NOT_RUN，没有晋升任何 window011 候选。
window011 独立读回退出0，但 window_completed=false、新增接受0，
末态仍为 window010 的 62.1875 ms。其零增量是零接受步的结果，不能拿来
断言材料停流。原 coarse 错误、候选和检查上限均保留。
后台总墙时 2567.25 s（约42分47秒），状态 PHYSICAL_STAGE_FAILED，退出2；
19:04:09 停止，原 PID52064 已退出，先前长任务完成时刻的估计失效。

只读全历史观察已另保存至
`outputs/deposition_observer_20261007/after_extended_010_v1/summary.json`，
退出0、1144真实帧、34.0814 s、无物理重算。8 个已穿孔粒子全部仍在方孔
投影内，距离投影边界约0.439–0.815 mm；肩下环区/实际肩接触质量均为0。
它们穿孔后有真实新增塑性应变，范围0.0475–0.1292；不能说完全没有材料响应。
还有67个初始合格粒子未穿孔，当前库存并未耗尽。

目前直接优先项是工具运动时支持邻接关系的遮挡/重新可见事件处理与受控
时间推进；延长总时长本身不能越过此次拒绝。完整沉积还需处理4 mm基板
小于5.4 mm肩径以及持续供料问题；这两项是适用范围限制，尚不能据此
单独认定本次没有铺展的因果根源。局部热态、材料/空间近似的旧限制仍保留。
