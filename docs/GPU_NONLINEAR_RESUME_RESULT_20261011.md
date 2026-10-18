# Nonlinear continuation on 2026-10-11

2026-10-15约17:56最新：短棒真实事件停跑已修复并恢复CPU续跑，详见
`CPU_SHORT_BAR_EVENT_REPAIR_RESULT_20261015.md`。三组同父容差和跨组比较全部通过，
首续跑段已73.763671875ms/393步，6粒新料穿孔2.025mg、肩下0。
实时目录`outputs/short_bar_event_repair_20261015_v1/run`，监督47984，原18:15:41截止。
下文停止/重启失败记录均为修复前历史，未改判或删除。

## 2026-10-15 五小时短棒续跑（事件失败后停止）

**17:11实际重启结果：** 用户再次明确“重启”后，新行政目录
`outputs/short_bar_restart_20261015_v1/`从完整stage035以原474.609375微秒步长重新执行。
15项行政测试、49段完整父链核查通过；17:08:33隐藏后台启动PID32040，实际Job无
kill-on-close。17:11:18首段stage036再次`GeometryPathStop`，监督exit 2、来源匹配，
最终全父链核验PASS后退出。没有新增完整段，检查点仍72.814453125毫秒/388步，
5粒穿孔1.6875毫克、肩下0，原完整NPZ实际SHA不变。新旧失败的path/router JSON及
`SOLVED_UNCOMMITTED_candidate.npz`的SHA分别完全一致，证实原样重启复现相同事件。
未使用失败候选、未放宽任何检查、未自动反复重启，18:15:41原截止不变。
新`restart_result.json`记录启动/结束与验证结果；下一步仍需显式适配短棒事件处理。

**约17:00最新复核：** 实际于14:35:50结束，状态`STOPPED_CHILD_FAILURE`、exit 2，
没有Python求解进程。16个完整阶段至stage035的72.814453125毫秒/388步，5粒新料穿孔
（1100、1116、1148、1164、1196），总质量1.6875毫克，肩下新料仍0。
末checkpoint SHA为`fc36a5ee9a80a643d0f0436005d1fc5c4c1c54ae801ca5db2bb64fffae57b4dd`；
根实际核对checkpoint、manifest、原监督与observer一致，原停止终审记录为PASS。

失败stage036的首coarse触发`GeometryPathStop`，支持边1059–1085发生真实
`BLOCKED_TO_VISIBLE`，起始严格工具交叠0.238572589微米；路径检查仅用2/256次查询，
shell连续检查已通过，故不是查询预算或超时。诊断弦fraction约0.275026596–0.275026597
仅定位未接受候选上的事件，不是已经推进的物理时刻；失败段没有合格细步，不能接续。
当前短棒CPU入口不含耦合事件求解，旧长棒GPU入口绑定原parent/receipt与旧授权时限，
不能直接替换。用户本次要求停则启动，但本轮未盲重启：须先把已有事件处理显式适配到
短棒完整父态并验证，再通过原物理和粗细检查后续跑。原18:15:41截止未延长，旧源码、
失败证据和完整检查点不变。只读核查证据为
`outputs/short_bar_five_hours_20261015_v1/resume_assessment_20261015_1700.json`。
下文运行描述为较早快照。

首次穿孔时间已独立只读核对：前一完整fine采样55.6875毫秒仍0，stage026/macro0001/
attempt00/fine1在55.9248046875毫秒首次记录1粒新料（ID1196，0.3375毫克），其observer
实际SHA匹配`bfed64710212e4775a09352db08590403a9aed57abfe08cba7c04dece906ac2a`。
56.162109375毫秒累计0.675毫克，56.63671875毫秒达到1.35毫克，并保持至59.009765625毫秒。
stage020–027全部68个保存fine采样的肩下新料质量均0。与旧长棒记录的55.3583984375毫秒
首次穿孔相比，短棒保存的首次观测晚约0.5664毫秒，未显示更早穿孔。此为保存状态的观测比较，
CPU/GPU与时间步路径不同，不是连续事件时间误差界，也不能据此认定棒长造成显著差异或沉积改善。

**约14:00状态更新：** 后台仍运行，8个完整阶段至stage027的59.009765625毫秒/328步，
stage028子PID42140继续，监督PID16220存活。短棒自身已4粒新料穿孔（1116、1148、1164、
1196），总质量1.35毫克，当前肩下新料0；不等同于肩下沉积或持续铺展已经建立。
根实读stage027监督exit 0、来源一致、manifest与NPZ哈希、observer均通过；checkpoint SHA
`219880b409e2ad4e3d9373a7cf03049759c8d05be6b6fa5f76e1435ce93e2bce`。
本次只读核查，没有重启或修改运行中的源码；18:15:41截止保持。

**13:29实际进展：** 新stage020已完整完成4个宏步，原粗细检查全部通过且无失败，
到48.396484375毫秒/268步；监督exit 0、来源一致，NPZ实际SHA
`54c6b8c26e90098c48eb5de1e6c80c95774d2133975103cde7e5d303a6f997e2`。
根独立读取末检查点及observer，新料穿孔/肩下均0，记录于新目录`startup_verification.json`。
程序已自动启动stage021子PID7668，监督PID16220继续运行；截止仍为18:15:41。
以下13:24内容为启动快照。

北京时间13:23:47实际启动独立后台监督PID16220；13:24:19心跳RUNNING_CHILD，
新首段stage020子PID42160，此快照尚无新完整阶段。`run.bootstrap.json`确认实际
kill_on_job_close=false及来源SHA匹配，`run.launch.json`保存启动命令与绑定来源。
截止仍为18:15:41，启动准备没有重置五小时。

新15项调度测试实际PASS，33个完整祖先的一次预检通过；独立核验末完整父态与新
行政代码无阻断，证据分别为新目录`prelaunch_checks.json`、`resume_audit.json`。
首段dt=0.5毫秒、12候选/最多4宏步、原预测437.918秒；原450秒预估上限与
540/600秒内部/外部限制、650秒外层监督和结束余量全部保留。原CPU物理来源、
配置与门槛未改，新campaign/test/launch在启动后冻结。

暂停方式为在本轮`run/PAUSE_REQUEST.json`写入请求，存在即锁存。正在计算的段
仍受原超时及总截止约束，完成后先核验并保存完整状态，再STOPPED_USER_PAUSE，
不再启动下一段；失败段仍不可提升为父态。该控制流已用调度测试验证，不代表新物理结果。

用户在前次暂停请求之后明确恢复五小时。本轮授权北京时间13:15:41至18:15:41，
准备和核验计入。检查时Python进程数0；旧两小时进度最后更新于14日22:05:57，
仍写RUNNING_CHILD，但没有退出和最终审计收据，因此不将其视为仍运行或正常完成。
上次暂停处理回合曾中断，现有证据不能确定后台进程的具体停止原因。

有效父态为旧`short_bar_two_hours_20261014_v1/run/short_4p0mm/stage_019`，
46.396484375毫秒/260步，监督exit 0且来源一致；checkpoint SHA为
`0c50016612dadb4a71cc1654c8916695797eafce40a7a11e08b9b0785ce398ab`。
stage020留下47.396484375毫秒/264步运行中记录，但缺完整supervision且未授权接续，
故不提升为父态。原12个完整阶段增加22.125毫秒，监督累计2731.591秒；旧文件保留。

新目录`outputs/short_bar_five_hours_20261015_v1/`复用已验证的CPU行政调度与原求解器，
原物理配置、来源和门槛保持。60与120.439208984375毫秒为观察节点，原220毫秒
stationary阶段末尾为诊断上限。原配置的comparison_physical_end_limit_s约120.439毫秒
标签保留；新protocol明确超过该点是本轮五小时授权下的短棒继续诊断，并非宣称已完成
更长的等时长短棒比较。未延长原stationary持续时间，也不进入走刀阶段。

## 2026-10-14 两小时短棒续跑（已启动）

**21:24实际运行更新：** 已完成stage008、009两个阶段，自动进入stage010；最新完整状态
26.021484375毫秒/172步，新料穿孔及肩下粒子均0。首个0.5毫秒配对未通过原时间误差检查，
冻结runner自动同父缩到0.25毫秒后通过并继续，未造成整轮停止。末stage009监督exit 0，
来源一致、原配对通过，NPZ实际SHA为`48603a0e5ac932a421557614748cb0794de5eca70d67d62fd8cd5f25cf7cb2cc`。
根独立读回证据保存在新目录`startup_verification.json`；监督PID51704存活，当前子PID52908。
截止时间仍为23:06:13，不因完成初始核验而重新计时。下列21:18内容是启动时快照。

用户明确要求续跑两小时。本轮从短棒24.271484375毫秒/158细步完整状态接续，
授权时间为北京时间2026-10-14 21:06:13至23:06:13，准备与核验计入，不另行延长。
原短棒CPU配置、成本记录与冻结求解器继续使用；原120.439毫秒长棒检查点保持。
60毫秒作为精确观察节点，剩余预算内继续，物理诊断上限120.439208984375毫秒。
新后台控制位于`outputs/short_bar_two_hours_20261014_v1/`，运行状态为`run/progress.json`。
北京时间21:17:49实际启动独立后台监督PID51704，21:18:04心跳为RUNNING_CHILD，
首段stage008子PID52556，此快照尚无新完整阶段。`run.bootstrap.json`实测Job的
kill_on_job_close=false且campaign来源SHA匹配；`run.launch.json`保存实际启动命令及来源。
后台为当前短棒原CPU路线，未切换计算实现；原120毫秒长棒不参与此次续算父链。

`campaign.py`和`test_campaign.py`的12项调度检查实际通过，完整初始父链预检exit 0；
证据见`prelaunch_checks.json`。独立静态审查核对绝对截止、单实例、完整prefix接续、
失败拒绝与60毫秒之后继续逻辑，无阻断。首段原预算quota15、最多5宏步，预测442.143秒；
恢复此前实际通过的0.5毫秒步长，末271.484375微秒仅为旧目标截尾，依据写入protocol。
每段仍由冻结CPU脚本执行原540秒内部/600秒外部限制，外层650秒监督及截止余量保留。
只有完整exit 0、来源/产物/物理配对检查通过的状态才能接续；未知失败不会跳过或改判。
原模型、配置和物理检查不变，后台控制源码启动后冻结；未声称此次已形成沉积。

## 2026-10-13 启动结束后的对照结果（已完成）

- **实际推进已完成。** 短棒从11.5625毫秒接续到20毫秒，再到24.271484375毫秒/158细步；跨过启动段并观察4.271毫秒恒定送料。7个新阶段、27个完整宏步全部通过，81次候选计算无新增失败，已按目标停止并保存。
- **短棒孔口推进量有所增加。** 在20→24.271毫秒这段，相同压头下移9.05555微米，孔口初始同坐标25粒组平均下移：长棒7.63596、短棒8.40218微米，短棒约多10.03%。从零累计到终点则分别为25.89331和27.31314微米，相差5.48%；两种百分比统计区间不同。
- **仍没有形成出料或肩下沉积。** 两组各读取该区间全部19个保存细步态，新料首次穿孔质量均0，当前及曾观测肩下新料质量均0。终点最近未穿孔间隙长棒188.915、短棒185.849微米；最近粒子身份可能改变，不能将此间隙当作单粒子的单调运动或用它线性预测穿孔时间。
- **当前判断：** 缩短棒料影响了局部推进，但尚未验证沉积收益。恒定指令仅观察4.271毫秒，不代表材料流动已经稳态，也不能排除后期棒长作用或判定原120毫秒问题根因。

| 观察量 | 6.5毫米长棒 | 4.0毫米短棒 |
| --- | ---: | ---: |
| 20毫秒孔口固定组累计下移（微米） | 18.25735 | 18.91096 |
| 24.271毫秒孔口固定组累计下移（微米） | 25.89331 | 27.31314 |
| 恒定送料区间孔口固定组平均下移速度（毫米/秒） | 1.78766 | 1.96704 |
| 24.271毫秒压头轴向反力（牛） | 8.98983 | 10.75147 |
| 24.271毫秒新料穿孔 / 肩下粒子数 | 0 / 0 | 0 / 0 |

上述位移、速度和压头位移比均不是质量流率或沉积效率。肩下记录是保存细步的几何采样，
不能排除采样间隙内未保存的短暂事件，也不代表冶金结合。长棒复用原GPU结果，短棒沿用
原CPU求解器且自适应步长历史不同；同父CPU/GPU重演只验证对应一步实现一致性，不能
用作后续累计误差上界。两组原粗细步检查通过，仍不构成空间收敛证明。

结果目录：`outputs/length_comparison_post_startup_20261013_v1/`。
`comparison_plan.json`保存预定目标与预算；`short_run_summary.json`保存7阶段完整链、
退出及源校验记录；`comparison_20ms.json`和`comparison_24p271484375ms.json`均实际
执行等时比较并exit 0。`constant_command_interval.json`实际读回两组各19细步，exit 0，
保存观察结果、文件哈希及采样限制；分析仅读取状态，没有新增物理解。
`reference_adapter_checks_v3.json`记录原GPU完整源/历史/产物核验、9项拒绝检查及原replay
自测通过；适配器与区间分析均经独立审查。

长棒参考为`campaign_10h_20261008_v1/attempt_0017`和`attempt_0026`，未作短棒父态。
短棒20毫秒保存于新`short_4p0mm/stage_005`；最终保存于新
`short_4p0mm/stage_007/macro_0004/attempt_00/fine_2/NOT_RESTARTABLE_candidate.npz`，
SHA256 `a1843b8115e1d5ee182085b98c008b6d3a925ba462c4d760c035892b13d6ba45`。
沿用原450秒预测/540秒内部/600秒外部限制，物理源、配置及检查门槛未改；本轮worker
1162.068秒、监督1225.928秒，任务至汇总1891.396秒，低于预定3600秒。
23:54:48 +08求解已退出，Python进程数0。后续须从完整合格状态经过同样检查接续，
不是将诊断候选文件冒充正式生产重启文件。

原主算例120.439208984毫秒/4粒新料穿孔/肩下0保持完整，检查点SHA
`987b3919c11a02a36f0c27facb80822ab3a6a1ea1f63a46907315b29ae225df0`已再次实读一致。
以上为最新状态；以下均为较早阶段记录。

## 2026-10-13 重点复核报告

- **计算记录完整。** 已重新读取最终比较，并核对14段连续保存链和304项文件校验值，结果一致；旧120.439毫秒检查点未改变。本轮仅复核，未增加模拟时间，当前没有求解进程。
- **缩短棒料尚不是已验证的修复。** 终点孔口固定粒子组平均位移仅增加0.235%，这是位移差，不是沉积量增长。棒长确实影响局部响应（孔内整体下移约多5.27%），但尚未带来新料穿孔或肩下沉积。
- **主要未解问题仍是持续出料与肩下铺展。** 旧长程为4粒新料穿孔、肩下0；新对照只到11.5625毫秒，仍在20毫秒启动段内，不能据此判定后期棒长作用或原长程问题的根因。
- **下一步重点：** 补齐送料速度稳定后的等时对照，使用新料穿孔质量和肩下累计质量判断改动是否有效，再决定是否调整主模型。

复核证据位于 `outputs/length_comparison_20261013_v1/`：
`review_comparison_11p5625ms_20261013.json`、`review_integrity_20261013.py` 和
`review_integrity_20261013.json`，两项实际读回执行均exit 0；后者记录52个完整通过的宏步、104个细步状态与171次候选计算的连续性。独立审查另外确认：跨组用初始坐标匹配，未误用不同几何下的相同粒子ID；最近未穿孔前沿的长短差仅0.1785微米。末宏步粗细差不能视为累计时间误差界或空间收敛证明。没有修改物理源、配置、阈值或旧结果。

## Executed independent bar-length comparison on 2026-10-13

The request to make the next effective advance has been executed as an independent
6.5 mm / 4.0 mm bar-length comparison through 11.5625 ms, under
`outputs/length_comparison_20261013_v1/`. This supersedes the earlier preparation-only
NOT_RUN status below. It is not a continuation, truncation or reset of the original
120.439208984 ms checkpoint: that original campaign remains stopped and intact.
Both new controls start at zero from the same prescribed all-hot initial state,
with dx = 0.5 mm, common 20 ms smooth startup and identical stationary feed/rotation.
The only changed resolved physical parameter is bar_length. Both new configurations
have a 1 mm future track, while the archived reference has 3 mm; this unused future
difference does not enter the stationary interval being compared. Particle counts
are 1381 and 1256, at the same resolution; this is not particle refinement.

The frozen remote CPU runner and its existing component receipt are reused.
Both cost stages passed: long/short measured total cost-stage times 4.502/3.593 s,
candidate forecasts 11.482/9.685 s, and peak working sets 1.114/0.935 GB (decimal).
These are initial map/tangent cost estimates, not proof of nonlinear step cost or
stability. Each actual stage retains forecast <=450 s, internal 540 s / external
600 s bounds, one CPU thread and all original nonlinear, global/local coarse/fine,
crossing-ID, geometry and source checks. Measured subsequent candidate costs,
including failed trials, determine the next stage quota; no physical threshold changed.

New long and short stage_001 both reached 3 ms with six complete macros (18
candidate evaluations), exit 0. The new/archived long initial states and 3 ms
endpoints have all 40 particle properties bitwise equal, as well as the same
physical parameters and saved machine/contact/energy state. All 134 common
source hashes match. This permits reusing archived CPU long states as explicit
read-only equal-time references, avoiding duplicate long simulations. They are
never used as parents for the short case. The archived CPU chain extends to
11.5625 ms; the later archived 20 ms GPU state is not part of the current CPU
comparison and is not covered by this bitwise claim.

Short stage_002 retained failed temporal-screen trials at 500/250/125 us and
recovered with a complete passing 62.5 us pair, reaching 3.0625 ms. Stages 003/004
then reached 3.25/3.6875 ms. Stage_005 reached 4.125 ms with all three macro pairs
passing, stopped at DECLARED_MACROSTEP_BOUND, and did not reach its requested
4.8125 ms endpoint. Its complete screened prefix is explicitly used for the next
bounded stage toward 4.8125 ms. Stage_006 reached 4.8125 ms and stage_007 reached
5.65625 ms, both with three passing macros and exit 0. Stage_008 then reached
6.5 ms. The explicitly bounded follow-up targeted the archived CPU endpoint at
11.5625 ms, with intervening equal-time comparisons and an additional 1800 s
wall ceiling; every child retains its original 450/600 s forecast/runtime limits.
Stage_009 preserved an original temporal-screen failure, halved from 281.25 to
140.625 us, and passed five macros to 7.203125 ms before the candidate quota
stopped the stage. Its qualified prefix was continued in a new directory.
Stage_010 completed 7.765625 ms; stage_011 completed 9.34765625 ms; stage_012
completed 10.61328125 ms. Stage_013 retained another temporal-screen failure at
474.609375 us, halved to a passing 237.3046875 us pair, and saved a complete prefix
at 10.8505859375 ms. Stage_014 completed two passing 355.95703125 us macros to
the final 11.5625 ms / step 104, exit 0. The bounded experiment is now finished;
there is no continuing physical worker. No failed coarse/fine state was used as
a parent. Normal temporal refinement inside the frozen runner was allowed;
unknown/source/geometry failures were not bypassed.

Final short checkpoint:
`outputs/length_comparison_20261013_v1/short_4p0mm/stage_014/macro_0002/attempt_00/fine_2/NOT_RESTARTABLE_candidate.npz`,
SHA256 `3678594a7eee7087714135d2b0b2ab7b574eec29cb9d61e3b2fc74a7add3727c`.
The short chain completed 14 stages, 52 passing macros / 104 fine steps and
171 candidate evaluations. Five failed temporal pairs are retained, with all
subsequent selected pairs passing. Ten stages reached their own requested targets;
four preserved qualified budget-limited prefixes that were subsequently continued.
Every stage has supervisor exit 0 and source_inputs_still_match=true. Stage worker
time totals 1902.4459061 s and supervisor time 1973.2948679 s; including the cost
run these are 1906.0387452 / 1977.3113103 s. These sums exclude separate preflights,
readback, plotting and conversation time. `short_run_summary.json` records the
per-stage manifest/supervision hashes and metrics, SHA256
`f08d99a77052d6862cdcaa69894372d9c4c16a76fa4543f58b7e6cd7c4901448`.

The final state remains an engineering diagnostic with accepted=false, restart_allowed=false;
the existing explicitly screened diagnostic continuation mechanism is unchanged.

`compare.py` verifies the original stage/parent receipts, source and checkpoint
hashes, and compares only equal physical times. Cohorts use initial xyz: 192 common
coordinates have different particle IDs across the independently initialized bars,
so equal IDs across cases are not equivalent material locations. Identity is fixed
within each case. The initial z = 1.75 mm outlet layer results are:

| Time (ms) | Common cap feed (um) | Long downward motion (um) | Short downward motion (um) |
|---:|---:|---:|---:|
| 3.0000 | 0.13237 | 0.04977 | 0.08139 |
| 3.2500 | 0.16716 | 0.08719 | 0.11989 |
| 3.6875 | 0.24125 | 0.19853 | 0.20932 |
| 4.8125 | 0.51966 | 0.45788 | 0.46583 |
| 5.65625 | 0.82347 | 0.72838 | 0.73672 |
| 6.5000 | 1.21899 | 1.07628 | 1.08582 |
| 7.765625 | 2.00016 | 1.76866 | 1.78102 |
| 9.34765625 | 3.31732 | 2.93389 | 2.93979 |
| 10.61328125 | 4.65496 | 4.11201 | 4.11349 |
| 11.5625 | 5.82455 | 5.12626 | 5.13828 |

The short-bar advantage shrinks from 0.03271 um at 3.25 ms to 0.01079 um at
3.6875 ms and remains below 1% of the long-bar motion at 6.5 ms. Thus the initial
64% displacement advantage at 3 ms is not evidence of sustained feed improvement.
Last-macro coarse/fine differences are retained as
local numerical indicators, not cumulative error bounds, confidence intervals or
spatial convergence. Both cases still have zero supplied crossing and zero supplied
shoulder occupancy in these early states. No deposition claim or final choice of
bar length is justified yet. At the final matched endpoint the displacement
advantage is only 0.01203 um, about 0.235%, compared with the initial 63.5% at
3 ms. It is not monotonic. Final cap reaction magnitudes are 8.61861 N (long)
and 9.44966 N (short), so short-bar reaction is about 9.64% higher. The upper
free spans have different initial lengths (3.0 / 0.5 mm between layer centers):
their absolute shortening must not be compared as the same strain measurement.
The final nearest supplied-center gaps remain 237.497 / 237.319 um, and neither
control has supplied crossing or supplied shoulder occupancy.

These observations do not support selecting a shorter bar as an established
deposition fix. They also do not rule out a later length effect: 11.5625 ms is
still inside the 20 ms startup, with only 5.82455 um imposed feed. The original
120 ms backflow/compression symptom, stationary sustained flow, spatial sensitivity
and deposition remain unresolved. Any later decision needs matched later-time
evidence; the current result is not extrapolated to the original long trajectory.

Executed readback commands use `.envs/pysph/python.exe -B
outputs/length_comparison_20261013_v1/compare.py --long-stage <long stage>
--short-stage <short stage> --output <new JSON>`; archived references additionally
require `--archived-long-reference`. Comparisons at 3, 3.25, 3.6875, 4.8125,
5.65625, 6.5, 7.765625, 9.34765625, 10.61328125 and 11.5625 ms exited 0.
Missing explicit archive opt-in, unequal times and a changed temperature parameter
were separately confirmed to be rejected. `archive_identity.py` executed exit 0
and saved the explicit long-reference proof in `archived_long_identity.json`
(SHA256 `48fca3e5c6768939dab0a82fe0eecf23171d3edd5c34caced2a29b1da8568575`).
`plot_comparison.py` creates new, non-overwriting PNG/PDF directories with a
source-hashed summary. The preliminary 4.8125 ms PNG and actual PDF page rendering
were visually checked; final figures include the final matched endpoint in
`plots_final_11p5625ms_v2/length_comparison.png` and `.pdf`, with source-hashed
`summary.json`. The completed short chain is summarized in `short_run_summary.json`.
The v2 layout places the upper-span annotation clear of the final data point;
`plot_comparison_v2.py` preserves the original plotting script, and
`visual_qa.json` records actual PNG and Poppler PDF-render inspection. Both final
formats passed visual inspection. The long/short endpoint step counts are 96/104:
times match, but adaptive time-step paths differ, so no cumulative-error equality
is claimed. A final direct hash of the original 120.439208984 ms checkpoint still
equals `987b3919c11a02a36f0c27facb80822ab3a6a1ea1f63a46907315b29ae225df0`.
The final process check reported zero Python processes. No solver source was
edited for this work.

## Layer motion and contact-force diagnosis on 2026-10-13

The user requested continuation of the diagnosis. The readback now covers the
complete restart parent and all 83 completed macro endpoints: 84 states spanning
89.589599609–120.439208984 ms. It does not advance physical time. New artifacts
are in `outputs/flow_diagnostic_20261013_v1/`: `analyze.py`, `diagnostic.json`,
`contact_readback.py`, `contact_readback.json`, `plot.py`, and
`flow_diagnostic.png/pdf`. Every checkpoint hash was checked, persistent identity
and mass were checked, and saved force candidates matched the endpoint positions
and deviatoric stresses. The contact readback verified force decomposition,
constraint reactions, active-contact geometry and the original residual limit.

The key distinction is that pressure/load **does reach the outlet**, but its
deformation response has not established sustained supplied-material outflow.
The earlier paired-position observation alone could not establish this distinction.

| Fixed initial layer z0 | Mean axial displacement over the interval |
|---|---:|
| 6.25 mm, actuator | -65.401172 micrometres |
| 4.25 mm | -14.433212 micrometres |
| 2.75 mm, initially above the tool | +19.503510 micrometres |
| 2.25 mm, initially in the bore | +46.691145 micrometres |
| 1.75 mm, initially near the outlet | +13.315211 micrometres |

Positive z is upward. Each row tracks the same 25 IDs, rather than particles
selected anew by current height. The average separation between the initial
2.75 and 6.25 mm layers decreased by 84.904682 micrometres. Upper layers also
increased their RMS radius measured about their own centroids; this measure
removes rigid rotation and centroid translation. It is not a cross-sectional
area or a mass-flux measurement.

All 175 originally upper free-bar particles have positive plastic-strain
increment in the latest fine step and endpoint temperatures 723.3805–723.7637 K
(450.23–450.61 C). Their mean equivalent-plastic-strain increment over the
whole observation window is 0.03072955. This confirms plastic deformation of
the upper free bar in the current hot-start model. Its plastic-work increment
is 0.001250134 J, only 5.6149% of the total 0.022264527 J. Therefore the geometric
shortening must not be described as most input energy being lost in the upper
bar; the rotating tool also supplies work and most plastic work occurs elsewhere.

At the final endpoint, the cap applies -10.916716 N axially, the tool applies
-3.665946 N to the material, and fixed substrate supports supply +14.582652 N.
Summing pairwise internal forces over material above horizontal cuts gives
upward resultants of 10.916708 N at z=3.0 mm, 10.168582 N at 2.5 mm,
10.278570 N at 2.0 mm and 14.582653 N at 1.5 mm. The corresponding opposite
loads act on the material below each cut. These are discrete partition forces,
not calibrated continuum traction integrals. The load path is not broken.

In all 84 selected final fine-step intervals, the net contact force on material
currently in the bore points downwards (-5.9063 to -3.4495 N). This does not
support a claim that the bore walls as a whole axially oppose feeding in these
samples. Individual contacts can have other directions. The four crossed IDs
have zero direct tool-contact force in every sampled interval; their inward
motion cannot be attributed to direct tool friction on those four particles.

Shoulder contact is inactive throughout these samples, consistently with the
saved finite-face distances and 250 micrometre contact radius. The all-particle
minimum centre-to-shoulder distance changes from 298.066 to 292.045 micrometres,
leaving 42.045 micrometres before contact at the final sample. This minimum is
nonmonotonic and changes ID: the final nearest ID 1189 is initial lower-bar
material, not one of the supplied crossed IDs. Its endpoint change cannot be
extrapolated to a contact time. The four crossed particles are much farther
from the shoulder face (minimum 834.070 micrometres at the end).

The geometric under-shoulder region contains 4–12 original lower-bar centres
across these endpoints (8 initially, 4 finally), while supplied-material
occupancy remains zero. Original lower-bar material was already present below
the outlet at initialization and must not be relabelled new deposition.
Neither the contraction of the four crossed particles nor their zero occupancy
proves that all old material stopped filling the surrounding region.

This narrows the useful follow-up to the deformation response of the hot free
bar and the actual spatial discretization near the outlet. It does not justify
changing friction, prescribing lateral particle motion, or changing temperature
inside the existing checkpoint history. The current configuration explicitly
assumes a fully hot deformable upper bar; guided/cold-feed alternatives need
an independently labelled setup and evidence. Source physics and existing
checkpoints were left unchanged.

Executed commands (all final invocations exit 0):

```text
.envs/pysph/python.exe -B outputs/flow_diagnostic_20261013_v1/analyze.py
.envs/pysph/python.exe -B outputs/flow_diagnostic_20261013_v1/contact_readback.py
python -B outputs/flow_diagnostic_20261013_v1/plot.py
```

The PySPH environment lacks matplotlib, so plotting uses the existing system
Python. No package installation was needed. The PNG was visually inspected.
Contact conclusions apply to the sampled final fine-step intervals, not every
coarse trial, event substep or unsaved instant.

### Actual saved-graph response probe

`gradient_probe.py/json` applies small diagnostic fields to the archived
endpoint graphs, volumes and correction matrices. It performs no nonlinear
solve, constitutive update, NNPS reconstruction or checkpoint modification.
The two-state probe ran in 0.394 s, exit 0; all source artifact hashes were
unchanged before/after. The independent affine-gradient error is at most
1.49e-15, its elastic-reference energy matches the analytic formula, and
reconstructed actual last-step symmetric gradients match saved D to 2.77e-15.

For the same within-region RMS of neighbor displacement differences divided by
separation, an initial-layer alternating axial field produces the following
elastic-reference energy relative to an affine axial field at the final state:
upper free bar 5.58%, initial bore region 3.48%, initial lower bar 6.74%.
After removing a best affine trend from the actual final step's axial
displacement, the selected alternating pattern accounts for 0.56% of upper-bar
nonaffine squared displacement and 20.72% of lower-bar nonaffine squared
displacement (lower-bar start value 3.88%).

These are descriptive operator/projection measurements, not actual plastic
energy fractions, a numerical FAIL threshold, or proof of an hourglass root
cause. Manufactured fields include fixed degrees of freedom, so this is not
an admissible constrained mechanical eigenmode test. It supports targeted
spatial sensitivity work in the lower region but does not explain most upper
deformation by this one alternating pattern. No stabilization was added.

### Prepared independent free-bar-length controls

`prepare_length_controls.py` was executed, exit 0. It creates two separate
zero-time configurations in `outputs/flow_diagnostic_20261013_v1/length_controls/`:
6.5 mm bar / 1381 particles and 4.0 mm bar / 1256 particles. The only differing
resolved physics parameter between them is bar_length, which also changes
the location of the same 0.5 mm top material actuator. Both retain a deformable
bar, 450 C initialization, the same tool, material, dx, cap depth and identical
rotation/feed programs. Both initial geometry/identity/inventory checks pass;
this is initialization evidence only, **not a physical comparison result**.

The common unused future track is 1 mm in both diagnostic plans so both finite
inventories pass planning checks. Through the stationary phase at 0.220 s,
command integrals match the original setup; the declared comparison end limit
is 0.120439208984 s. The test is explicitly independent: the existing deformed
state is not cropped, cooled, re-IDed or restarted as a shorter bar. The
production/original 3 mm plan is unchanged.

`length_controls/preparation.json` records the exact configurations and audit
results. Nonlinear cost measurement and physical comparison are **NOT_RUN**.
The next execution should first measure cost and pass the original force,
geometry, energy and time-screen checks for each setup, then compare layer
motion, reactions, upper shortening and supplied crossing/shoulder metrics at
matched physical time and cap displacement. A changed response would identify
sensitivity to free-bar length; it would not establish that the shorter bar
represents the real machine. No long-run authorization was extended here.

## Feed-to-outlet readback on 2026-10-13

The user asked for the first next step. A bounded read-only comparison used
the complete restart parent (old event-resume attempt 39, 89.589599609 ms)
and final nonlinear-resume attempt 86 (120.439208984 ms). Both checkpoint
hashes matched their manifests; persistent particle IDs, masses, eligible IDs
and crossed IDs agreed. No physical step was run and no frozen source was
changed. Executed command (exit 0):
`.envs/pysph/python.exe -B outputs/flow_readback_20261013_v1/readback.py`.
The script and `readback.json` in that directory retain the cohort definitions,
source hashes and measured values.

Both states lie in stationary_feed_attempt (20–220 ms): axial feed is
2.12 mm/s downward, rotation is 300 rpm, traverse is zero. The 25 actuator
particles have prescribed z motion and free x/y motion. During the 30.849609375 ms
interval every actuator particle moved down 65.401171875 micrometres, matching
the command integral to 2.13e-17 m. Final actuator height is 6.015868877 mm,
still above the 2.5 mm tool top. This interval does not support an explanation
that the actuator stopped feeding or reached the outlet.

The fixed cohort of 24 eligible uncrossed particles initially at z=1.5–2.0 mm
had mean upward displacement 16.529569 micrometres, with individual displacements
from -118.837151 to +86.418632 micrometres. The four previously crossed IDs
1151, 1203, 1229 and 1281 all ended closer to the stationary tool axis:
radial changes were -78.075252, -108.511395, -49.437623 and -32.059530 micrometres,
respectively. Their mean axial displacement was -0.322499 micrometres, with
substantial individual upward and downward movement. Crossed count remains four;
shoulder occupancy remains zero in the campaign observer.

These are two-state displacement observations of fixed cohorts, not net mass
flux measurements or proof of monotonic motion. They narrow the next diagnostic
to deformation/force transfer through the free bar and outlet contact. They do
not establish a defective contact law, jamming, or inadequate particle count as
the root cause. Next inspect layerwise axial displacement, radial size, outlet
velocity/stress and contact reactions before selecting a model change or a
resolution comparison. No new long run was launched.

## Completed outcome, inspected on 2026-10-13

No Python simulation process is running at inspection. The latest session
finished normally at 2026-10-11 23:40:32 +08:00 with STOP_WALL_RESERVE, exit 0,
83 complete passing stages out of 86 attempts. The decision had 151.111210 s
remaining versus the measured minimum starting reserve of 190.775297 s;
after the final audit, 74.402128 s remained before the authorized 23:41:47
cutoff. This is a budget-reserve stop, not an unresolved child failure.

Latest complete stage is attempt_0086, physical time 0.12043920898437391 s,
step 765. Actual checkpoint SHA256 agrees with the final summary:
`987b3919c11a02a36f0c27facb80822ab3a6a1ea1f63a46907315b29ae225df0`.
The continuation advanced 30.849609375 ms from its complete restart parent
at 89.589599609 ms, over 10887.895334 s of session wall time including startup
and final verification. Attempts 11 and 29 were nonlinear nonconvergence
rejections, while attempt 77 failed the temporal pair screen. All three were
retained, halved from the same complete parent and subsequently recovered.

Final supplied crossing count remains four (1.35 mg). Current and ever-observed
shoulder occupancy are both zero. Nearest uncrossed PID is now 1242, gap
16.364646019 micrometres; it is not the prior PID 1099, so changes in the nearest
gap must not be interpreted as one particle's monotonic progress. Maximum radius
among the four crossed centres is 0.496653755 mm, with no shoulder deposition
demonstrated.

The final history audit reports 38117 unique, expected and locked files, one
final byte audit and failed=false. The original authorization has expired.
This inspection did not restart or extend the computation. General automatic
checkpoint discovery/restart and Linux portability discussed previously remain
unimplemented; saved qualified evidence remains available for a later explicitly
authorized continuation.

## Hardware and server-readiness assessment at 21:48 +08:00

The session remains active, supervisor 9848 and child 45556 alive. Latest complete
attempt 25 is 0.0962341308593744 s, step 647, with 24 passed stages and the same
four crossed particles and zero shoulder occupancy. This demonstrates preserved
advancement from the 0.08958959960937454 s restart parent, rather than a restart
from the initial physical state.

Attempt 25 worker time is 78.6321105 s. GPU inverse/solve calls sum to
5.4982714 s, including host/device transfers and synchronization in the actual
instrumentation. If that entire measured portion were halved while all other
work stayed unchanged, this worker would save only 2.7491357 s. This is a scoped
estimate, not a benchmark of another GPU. Matrix assembly, nonlinear control,
geometry and validation retain CPU paths; selecting hardware requires a common
checkpoint benchmark of the whole stage. This follows the bottleneck-first
approach in the [NVIDIA CUDA performance guide](https://docs.nvidia.com/cuda/cuda-c-best-practices-guide/index.html).

Existing screened complete checkpoints preserve prior progress, but the current
entry is not a general unattended recovery service: its source parent and
deadline are explicit fixed values, and only qualified failure classes retry.
It does not automatically discover and restore any latest complete stage after
an arbitrary crash. The CUDA DLL binding, live history-lock implementation,
detached launch and process termination paths are Windows-specific. An unchanged
copy cannot be claimed ready for a Linux server. A future server preparation
needs general checkpoint discovery, interruption/restart tests, and platform
adaptation with equivalent physics checks. Scheduler requeue alone is not
application recovery: [Slurm documents that requeue restarts the batch script
from its beginning](https://slurm.schedmd.com/sbatch.html#OPT_requeue).

This assessment did not deploy, restart, alter the frozen implementation or
extend the current 23:41:47 authorization.

## Running observation at 21:43 +08:00

Supervisor 9848 and current child 10828 are alive. At the 21:43:11 snapshot,
21 complete stages passed and attempt 23 is running. Latest complete attempt 22
is at 0.09481030273437444 s, step 641, checkpoint SHA256
`bd3933209a24706aa6f6499a936ae4ddc35ef85da8c20a989bdb9d43c9b394b5`.
Independent actual-NPZ, manifest/log/screen/checkpoint hashes, original global,
local and crossing-ID gates, supervision exit and unchanged 636-source checks
passed. The physical advance from restart is 5.220703125 ms and 42 fine steps.

The new recovery has now worked in actual computation: attempt 11 failed at
474.609375 microseconds, recording nonlinear line-search nonconvergence; the
same complete parent was retried at 237.3046875 microseconds and attempt 12
passed. After cooldown and small-error results, attempt 22 returned to
474.609375 microseconds and passed. The failure remains retained; no tolerance
was relaxed. There is only one rejected attempt among the first 22 attempts.

The same four supplied IDs remain crossed, total mass 1.35 mg. Shoulder
occupancy and cumulative observed shoulder residence remain zero. The nearest
uncrossed particle is PID 1099 in both the restart parent and this snapshot:
its gap increased from 46.337873584 micrometres to 54.402920892 micrometres.
The maximum radius among crossed centres is 0.551078520 mm, inside the square
bore projection whose half-width is 1.252 mm. These observations establish
continued numerical advancement but do not yet establish sustained new outflow
or lateral spread beneath the shoulder. A gap extrapolation to zero is not
supported by the current motion.

Estimated startup-to-first-stage overhead was 1611.372 seconds (26.86 minutes),
computed from process start to first receipt timestamp minus that stage's
reported wall duration; it includes restoration and bookkeeping, so it is not
the preflight timer alone (649.226 seconds). Through attempt 22, physical
progress averaged 8.540 ms per wall hour in the continuation loop and 4.930 ms
per wall hour including startup. Recent receipt completion intervals were
98.65--103.88 seconds. Attempt 22 reports 77.870 seconds in the worker and
5.310 seconds summed across recorded GPU linear-algebra groups; this is scoped
instrumentation, not a full CPU/GPU profile. It shows why GPU availability alone
does not imply proportional speedup of the complete campaign.

Current geometry remains dx=0.5 mm, bore width=2.504 mm, 1381 particles. No
spatial-refinement comparison has been performed in this continuation, so
particle-count sensitivity remains unresolved. This status check only read
current evidence and updated documentation; computation and the 23:41:47
deadline remain unchanged.

## Verified launch at 20:39 +08:00

The detached supervisor started at 20:39:04.702537, PID 9848, in
`outputs/gpu_acceleration_20261007/nonlinear_resume_20261011_v1`.
At 20:39:19 it is alive in INITIAL_FULL_HISTORY_VERIFICATION with no GPU child
or new completed physical stage yet. Bootstrap is PASS_DETACHED_JOB_POLICY,
kill_on_job_close=false, limit_flags=0. The absolute cutoff is 23:41:47;
the actual new process initially has about three hours remaining.

Unified component validation passed all 37 tests (22 administrative restoration
and binding tests, 15 nonlinear-failure recognition/refinement tests).
Both groups independently passed review. The actual saved-state fixture
restores 120 observation frames, 59 complete cost rows, six retained failed
cost records and the exact prior STOPPED policy, then returns READY at the
same parent with 237.3046875-microsecond dt. There is no new physical PASS in
that fixture. A preparatory check verified 4,011 historical artifacts, and
the formal component checker reverified the inherited closure and unchanged
636-entry physical implementation.

The seven new administrative source files and component receipt are frozen.
Receipt: `outputs/gpu_acceleration_20261007/nonlinear_resume_components_20261011_v1.json`;
SHA256 `610aa876a680df1d626b77c8e6e9e870d54858110fbdf924a47f57f96c1ccb8a`.
Controller source SHA256 is
`50c8568beba4dfb23971edeb3946859cf0dba986f328a633bf7b0f59ba67755d`.

Executed launch command (exit 0):

```powershell
.envs/pysph/python.exe -B scripts/run_gpu_nonlinear_resume_20261011.py --output outputs/gpu_acceleration_20261007/nonlinear_resume_20261011_v1 --wall-seconds 28800 --max-stages 100000 --detach
```

The wall-seconds option is capped by the remaining original absolute deadline;
it does not authorize an additional eight hours. The process will perform its
startup history verification before attempting the real GPU step. Whether that
real smaller step converges remains unverified at this launch observation.

## Preparation and retained evidence

The user explicitly requested restarting the current computation after the
nonlinear failure. The original authorization cutoff remains 23:41:47 +08:00;
the continuation does not reset the eight-hour budget.

The complete parent is
`outputs/gpu_acceleration_20261007/event_resume_8h_20261011_v1/attempt_0039`,
at 0.08958959960937454 s, step 599, with 1381 particles. Checkpoint SHA256 is
`aedf8bfed9c67f38ec12bf38e1e9ec77e37cbba39914f63064dfedbfb537c474`.
Original global/local/crossing-ID checks and direct checkpoint hash verification
passed. Four supplied particles have crossed, mass 1.35 mg; shoulder occupancy
is zero.

The rejected next coarse proposal exhausted the original Broyden Armijo refresh
budget and the original Newton line search. Its final force residual was
6.386702887629103e-4 N versus the unchanged 1.0502203271584808e-5 N limit.
The failure remains recorded and no rejected state becomes a continuation parent.
New administrative recovery is being prepared to retry from the full parent at
237.3046875 microseconds, retaining all original candidate and pair checks and
the unchanged 636-entry physical implementation.

The previous audit's one extra locked file is accounted for by the final
`attempt_0040_failed_cost.json`: it was hashed and locked before STOP, with no
later stage to register its expected hash. The final byte audit still reread
all 30284 locked files. The new restoration explicitly registers that cost
record; it must not generally permit unexplained audit-count discrepancies.

At creation of this note, unified tests, component freeze, and background launch
are pending. No new physical advancement has yet been established.
