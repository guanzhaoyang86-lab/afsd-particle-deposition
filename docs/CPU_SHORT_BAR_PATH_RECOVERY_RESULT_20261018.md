# 事件后路径证明修正与一小时续跑

约11:52 CDT核查最终停止：监督44628实际已11:29:18.101退出（exit2），
本轮46分51秒，完成64pair至87.0527343749995ms/586步，较起点推进7.341614ms。
末SHA `c7b99e13974f40419aca9c28948303e9a08ab929d769379c5bfd0822f7844c70`，
root用原新verify_pair对220源、收据及粗细状态实际重算PASS，仍7穿孔2.3625mg/肩下0。

真正首错在piece_000065_attempt_000001/fine_2/proposal：Broyden的Armijo搜索在
3次Jacobian构造/2次refresh、37次评估后失败，Newton fallback也发生
complete residual line search/evaluation budget failed，封装为RefinableCandidateError。
分支监督46.906秒，timed_out=false，源与manifest/protocol SHA均一致。
求解未完成所以没有proposal/visibility_event_routing.json；controller retry_evidence
误将它按多开启几何事件读取该文件，导致二次FileNotFoundError覆盖首错。
这不是本轮旧几何包络问题复发、也不是累计减半次数限制或硬截止。

需要补充有完整证据约束的未收敛整pair回滚减步重试和缺文件的正确分类，
不能使用失败fine2或放松原物理门槛。证据v5/status_20261018_stopped.json。
本次仅核查记档，未启动或改冻结源。原截止11:42:26 CDT已过，细粒子仍只讨论目标，
未生成/启动新的0.25mm物理解。

11:23 CDT状态复核：用户怀疑停跑，实际监督44628仍运行，进度从58段fine1更新到fine2，
没有run.exit.json或错误。已完成57个完整段，末86.22216796874951ms/572步，
正常宏步长118.65234375us已恢复。root再次独立verify_pair/220源/原粗细状态重算PASS，
末SHA `d428d7b9deac2b37dfb75b4319e66790498de20e074bbef443a492e83fe12fe7`。
仍7穿孔2.3625mg/肩下0；证据同v5目录status_20261018_1123.json。
本轮仅读状态，未重启或变更粒距，截止11:42:26 CDT及开段650秒余量不变。

10:46 CDT实际成功更新：新首piece_000001_attempt_000001完整粗细pair已通过，
79.72595214843717ms/460步，13次物理解，178.424874秒。
checkpoint SHA `223dfaa007343697d97dad063de4cadfbdd3f6f5225c1eb74c1ff991c11caf1a`。
root用新verify_pair独立实读220源、全部收据及原粗细状态重算PASS，证据
`outputs/short_bar_event_continue_20261018_v5/independent_first_pair.json`。
实际fine1事件后路径输入SHA与旧v4失败路径逐项相同，新35queries完整证明，
故本次已真实跨过旧阻断，而非只有几何回放或启动成功。后台随后第二pair完成，
progress保存至79.75561523437466ms/462步，第三段coarse子15632已启动；监督44628继续。
仍7粒穿孔2.3625mg、肩下0，未宣称形成沉积。截止11:42:26 CDT不变。
启动实际220来源与旧207完整继承证据见同目录startup_verification.json。
以下为启动阶段记录。

2026-10-18用户要求继续跑，沿用上一轮一小时预算。新v5于UTC15:42:26.885103
（America/Chicago 10:42:26 CDT）启动，固定截止UTC16:42:26.885103（11:42:26 CDT）。
命令`.envs/pysph/python.exe -B outputs/short_bar_event_continue_20261018_v5/launch.py --wall-seconds 3600`
父启动exit0，后台监督44628，首coarse子2676已进入实际求解。启动收据实际
PASS_DETACHED_JOB_POLICY_AND_SOURCE，kill_on_job_close=false。本次检查时首个完整pair
尚未通过，不能把几何复核或启动成功当作跨事件后的物理进展。

起点为旧v4首完整段79.71112060546842ms/457步，SHA
`e93edb7ca310b551f7a76c2f3636ce5ef8b3172bc64cf7506622001bbf094793`。
原verify_pair重读全部207源、原artifact及粗细实际状态PASS，7粒穿孔2.3625mg/肩下0。
旧失败piece2和全部旧源不变。恢复三条待完成绝对时间边界，nominal118.65234375us，
下一小段14.83154296876us；只有完整粗细pair通过才消费边界、更新检查点。

原阻断：1091–1173开启后的保留边连续路径，原通用Rodrigues矩阵范数造成
6.21345e-13m误差包络，大于保存的起点约5.7604e-13m投影间隙。新模块保留原严格
相交/退化/数值域/曲率/旧包络检查，仅在原棱柱证书失败时补充角度相关范数界限，
显式保留绝对sin/cos及角度运算误差，eta仍8192epsilon。推导与局限见
`VISIBLE_SUPPORT_ANCHORED_CONTRACT_20261018.md`，仍为工程浮点包络而非形式化区间证明。

真实旧失败30499边同输入SHA复核：旧256查询/深32未证，新35查询/深5全部证明，
新增分离证书最小余量1.920946e-14m；0物理解、0NNPS、旧artifact未改。
实证`outputs/short_bar_event_continue_20261018_v5/saved_geometry_replay.json`及full.json，
可复现脚本同目录`replay_geometry.py`。该旧弦未作为新的物理状态升格。

root最终199项组件测试全部PASS，unittest3.341s、整段3.713s、exit0，证据
`outputs/short_bar_event_continue_20261018_v5/component_checks.json`。包括125旧回归、
20几何（含15旧fixture重定向）、4CPU/receipt边界AST核验、33恢复行政与17启动器。
独立源码/推导审查无阻断。另有普通双精度数值交叉检查，不称为高精度证明。
新CPU driver只变路径provider与显式协议；新reader只变显式driverimport；原CPU
力热/本构/非线性/定位/blocked-prefix和粗细阈值未修改。

新配置SHA `06f49cf0dd35d5fd0f04ea1d230f6a766cd10ce6746ae6043d0d9f441dbf1c29`，
controller SHA `8dc31d829f7fb257f1aaeab87af11b20e94193560d9dbdca5e7e3d2bcc7662a9`。
每子450预测、540/600限制及开段650秒余量不变；可能早于硬截止收尾。
当前完整成本102样本预测366.615s，未采纳失败piece2的未完成细分支。
进度在新run/progress.json，暂停用新run/PAUSE_REQUEST.json；勿改冻结源或重复启动。
