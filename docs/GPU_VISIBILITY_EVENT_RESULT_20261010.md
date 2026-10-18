# 遮挡激活事件修复与验证记录

2026-10-10，用户要求解决十小时任务旧attempt_0044的真实支持边激活失败。原始最后完整checkpoint为旧43：71.37646484375 ms、480步、SHA `3e9caf7d52dc71b02f7ea92e6d5d93b25ab804162992de1d499e7ca61638c1dd`；4个新料穿孔、肩下0。旧44失败及所有来源冻结保留。

## 实现与边界

合同为`docs/GPU_VISIBILITY_EVENT_CONTRACT_20261010.md`。新变体用同父真实求解定位blocked/visible两侧，以统一支持图延迟指定incoming边激活；每个所用图进入原密度、材料梯度、力、切线、热传导及全部原物理检查。位置容差5e-8 m，时间容差min(1e-7 s,宏段dt/1024)。另对选中高端声明路径证明incoming持续blocked至最后容差窗口，随后真实子步恢复该边，保存逐个真实细子步与不可变历史。

这是显式工程近似，不宣称零延迟精确混合系统解或任意邻域拓扑事件都能处理。所有状态仍diagnostic_only、accepted=false、restart_allowed=false，正式接受增量0。保留450秒预测限、525秒worker截止、540秒硬限及600秒外部限。首次事件阶段为新模型的资格和成本实测，旧普通候选耗时不能证明新事件成本；实际超限仍停止。

## 已执行组件验证

`outputs/gpu_acceleration_20261007/visibility_event_components_20261010_v1.json`：116 tests，exit 0，PASS_INCREMENTAL_COMPONENTS。包含真实PySPH制造场的力/热隔离与恢复、守恒、真实失败支持边回归、解析运动定位、连续blocked前缀正反例、不可变子步与事件收据负向、原控制器限额和源码代际检查。未把制造场当作AFSD沉积验证。

另以旧43三个已保存分支的真实起点NNPS和终点邻域并集只读核验：三个分支均无观察到的支持边切换/未解析事件，可用于无事件重演。

## 真实GPU验证

启动目录：`outputs/gpu_acceleration_20261007/event_validation_20261010_v1`。
监督PID64124，2026-10-10 17:40:52启动；DETACHED/BREAKAWAY，实测kill_on_job_close=false。原绝对截止2026-10-11 01:37:12.059551+08:00不延长。

18:05结果：该验证任务已停止，监督已退出，最终全历史20,070文件审计failed=false。真实结果如下：

1. 旧42重演完整旧43：PASS，3个分支各40个粒子属性的shape/dtype/数组值均独立核对逐位相同；源、父SHA、历史、步数与宏alias逐字节绑定通过，worker108.1766秒。新完整父态为attempt_0001，71.37646484375毫秒/480步。
2. 旧44事件区间：attempt_0002 FAIL，worker333.7184秒。粗段locator耗尽10次真实求解，仍未得到合格宽度，fine未执行。full masked候选的弦提示168.727794微秒，而真实小步至174.487794微秒仍blocked。原严格交叠量从168.682794微秒时80.079901纳米降至174.487794微秒时74.783928纳米；弦提示不能当真实物理事件时间。
3. 同父容差1/2和1/4敏感性检验：未执行，不能基于失败粗段继续。
4. 长跑：尚未恢复，所有失败候选不作为父態；最新实际穿孔4、肩下0。

## 真实交叠量割线修订（已执行，前缀证明预算失败）

新增`run_gpu_visibility_event_secant_20261010.py`与单独locator，通过已完成的blocked端点严格交叠长度下降率给出下一次真实求解时长。只改时间建议及同父masked探针的非线性初值，不插值状态，不放宽10次定位/时间/位置/原物理几何门槛。普通分支和原gate保留AST；10项新driver、23项locator、53项行政测试已PASS。全stage计数上限明确登记为25，覆盖粗细两路各10次定位与3次普通提议/2次后侧求解，墙钟原限仍绑定。具体见`GPU_VISIBILITY_EVENT_SECANT_20261010.md`。

新行政入口一次完整历史核验后，同父执行1/.5/.25并从选定.25完整末态续跑。新的完整pair经验成本估算包含全部probe/定位/检查/worker/监督耗时，保留原50%裕量、450秒预测限及525/540/600硬边界；该经验预测不保证未知未来成本。

正式新收据`event_sensitivity_components_20261010_v1.json`，86 tests PASS，全部新来源冻结。`event_secant_sensitivity_20261010_v1`于18:09:59启动，监督PID51092，bootstrap PASS/kill_on_job_close=false；18:10处于INITIAL_FULL_HISTORY_VERIFICATION，进程活跃。截止仍01:37:12.059551。不重复启动。仍需真实tol1/.5/.25结果，不能称已跨越真实事件或形成肩下沉积。


18:25真实结果：监督51092已退出，status=STOPPED_SENSITIVITY_CHILD_FAILURE，全历史最终审计完成。sensitivity_one worker208.14秒、监督211.59秒；8次真实探针得到low249.367323865微秒/high249.457323865微秒，括区90纳秒、实际状态最大差4.009508纳米，通过定位。high自身router事件延迟约94.58纳秒，也满足100纳秒。但incoming_blocked_prefix在256总查询时UNRESOLVED：104 strict、101 box、51 contracted，仅50个原box成功证书，覆盖至0.9995533718124937，距目标0.9995991298293011尚差11.414672纳秒对应路径；无可见样本、退化或数值域失败。没有after_event、完整coarse/fine/pair，不能续接该状态。

## 原包络盒查询效率修订（待真实运行）

新增afsd_event_blocked_box_20261010只省去无成功证书的contracted备用查询，沿用原strict/fixed-parameter box证明、浮点/退化/深度32/预算256。同真实saved high、同input_sha256以224次获证，113 strict/111 box、56连续证书完整覆盖；越过事件至1.0拒绝，预算223/深度17拒绝。21几何tests PASS。新run_gpu_visibility_event_box_20261010与secant版所有函数AST完全一致，仅proofprovider/注册路径更换，3driver测试PASS。新行政三尺度同池入口待统一组件冻结和实际运行。当前最后完整物理时间与穿孔/肩下计数均未变。


18:38成本复核与新行政登记（尚未启动）：旧450秒预测门槛配合30+1.5×完整pair，要求实际完整pair≤280秒；而真实coarse未完成前已耗时208.14秒，历史无事件fine1+fine2也约55–59秒。box只修查询预算，原证明耗时约0.04秒，不能解决这一墙钟资源缺口。新v3明确将行政预算登记为900秒预测上限、960秒worker软截止、990秒硬定时和1080秒外监督，原十小时绝对截止不延长，所有物理/几何/精度门槛原样保留。另发现原candidate_forecast_s把含全部probe的宏分支成本用于一次Newton重试预留，量纲不符；新driver在每次完整求解保存后写只读physical_solve_cost.json，绑定父指纹、checkpoint/SHA和实测单次耗时，供行政入口按单次求解记账。删去该日志语句后所有物理函数AST与冻结secant版完全一致，更新后的3项driver测试PASS。完整pair成本仍包含全部probe和监督，不丢弃昂贵历史以压低估计。新行政测试和独立复核完成后才统一冻结。


18:43统一组件收据event_box_components_20261010_v1.json：95 tests PASS，独立审查71行政+3driver复跑PASS、341项旧闭包与v2完整失败SHA一致。新来源全部冻结。event_box_campaign_20261010_v1后台启动，监督PID61728，先全历史核验再执行同父三尺度敏感性与续跑。实际结果尚待，不把组件PASS当真实事件完成。


18:58:30 v3真实结果：STOPPED_SENSITIVITY_CHILD_FAILURE，监督61728退出，worker245.184秒/监督248.847秒。8探针定位、90纳秒/4.009508纳米、低/高端资格与224查询blocked前缀全通过；两端40粒子属性与v2逐位相同。真实after_event/proposal恢复PID1334–1371支持边，活动边33720条，post单次求解19.884秒。随后原可见路径证明耗尽256查询，found/degenerate/numerical_domain均空，shellCERTIFIED、两router无事件；整段仍FAIL，无fine/event_resolution/pair_screen，不能报告完整跨事件成功。10份single-solve日志全与实际同目录checkpoint/SHA绑定，最大51.448707秒，旧证据完整保留。

19:06诊断：恢复边在工具z=0平面距孔边约22.9皮米，整段末约71.7纳米；原全向体坐标速度包络5.33156微米/全fraction主要为切向运动，迫使证明极细分箱。同保存弦前半至1/2048均256耗尽，前100纳秒248次通过；这种插值仅几何诊断，绝非状态接续。准备v4第二阶曲率四点凸包/SAT补充证明，四个保守wall矩形分别严格分离、浮点/曲率包络向外，原strict/一阶proof/256预算保留。新driver仅改可见路径proofprovider，所有物理解和25/10计数原样，行政900/960/990/1080及原绝对截止不变。未启动v4，须真实失败路径与制造几何测试、独立审查后冻结再运行。


19:11 v4新proof真实post33720边同inputSHA以53/256查询获证，19strict+17inflated+17SAT、depth7，全部边[0,1]连续覆盖，旧44真相交仍2strict拒绝。统一event_prism_components_20261010_v1.json：95 tests PASS（15proof/3driver/77admin）；独立数学审查及542项v3闭包fresh-byte复核PASS。全部v4源冻结，event_prism_campaign_20261010_v1启动监督PID58764。先全历史核验再同父三容差实际验证与选定续跑，绝对截止不延长。当前仍待实际完整pair和后续敏感性结果，未形成新的可续接物理末态。


19:30首次完整真实事件验证PASS：sensitivity_one原粗细全局/局部/crossing IDs全通过，18次真实求解，coarse/fine1/fine2实际子步2/1/2、事件1/0/1，worker439.167秒。完整细末态71.85107421874992毫秒/483步，SHA8fae514013b72eb4e28a965fb4a06a456561e6fb82f79e148b6ca423b72a9fd6。全局位置差7.277271纳米、温度差4.594883微K、应力相对差0.0518069%；局部最大应力差0.6872998472293538%（PID1321）。coarse真实post所有粒子属性与旧v3保存post逐位相同，来源/父态/能量机器历史增量与恢复边独立核对正确，新path53queries通过。fine2为5探针、85.264纳秒/3.798纳米，blocked prefix124次通过。监督已进入sensitivity_half，同父比较仍待half/quarter与全部跨尺度screen。此时不能把三个实验累计为三段物理进度，也不能说常规长跑已恢复。


## 19:33 第二档实际停止与下一版收敛修复

prism campaign 已停止，one完整PASS保留；half worker207.614s，coarse 8探针真实括区45ns、位置差2.0047679665nm，但hi自身路径最早可见后的延迟60.086032082ns超过50ns。lo249.39009304025963µs、hi249.4350930402596µs，router start .9997591115534306/end .9997591124847531。原后置门槛正确拒绝；未执行后续prefix/after/fine/quarter，失败不作父态。

新delay版本正在让定位器同时收敛原括区/位置/自身延迟三个条件，窄括区delay未过继续真实同父二分，原10定位/25总solve及全部物理gate不变。冻结旧源和证据，新版本需测试和实际三个容差完整验证，当前没有运行中的GPU长跑。


## 19:42 delay v5 组件封存与实际启动

统一121tests PASS（35定位器、4驱动桥接与前缀资格、82行政）。独立审查复跑39物理侧tests PASS。旧源/收据全部原SHA不变，新event_delay_components_20261010_v1.json已封存。回调提供保守自身delay及原blockedprefix证明成功资格，locator四条件同时达标才返回；原独立后置delay/prefix和强收据验证仍保留。原256纯查询预算且无可见/退化/domain错误时只允许剩余10定位内再求真实解，不延长单次证明预算或接受插值状态。

新event_delay_campaign_20261010_v1监督PID60096已启动，sourceSHA7b48355a753f3373a2dfdd24fa6cddf40e6ccf67d123844cb4cdabcdf75ae278；先一次完整历史核验，再同父三个尺度新GPU实验、全部跨尺度比较，选quarter并继续原授权长跑。900/960/990/1080秒边界、25真实solve、原绝对截止2026-10-11 01:37:12.059551+08不变。此刻尚无新实际GPU结果，不能宣称已通过新half/quarter或恢复普通续跑。


## 20:00 delay v5 第一档完整 PASS

初始全历史628.589593s PASS（20069文件、2375913哈希请求、666JSON、389checkpoint），同池继续。one完整18solve、worker447.155887s、监督约452.07s，候选分支234.587989/24.349869/178.083729s，细末态71.85107421875ms/483步，SHA9f540046a9b56bd02aeb982b86b862b931f9886b4338853d4a82f8cf6742a3db。独立636source+447artifact SHA、strong receipt、原全局局部屏均PASS，3终态各40个粒子数组与v4 one逐位相同。coarse 8探针90ns括区、4.009508nm、delay94.578675ns≤100、prefix224；fine2 5探针85.264414ns、3.798123nm、delay40.264423ns、prefix124。两个post均恢复1334–1371且全33720边prism53queries获证；前置/后置prefix证据除walltime外完全相同。half已20:00:32启动PID6556，剩余两档及比较尚待完成，普通续跑尚未开始。


## 20:08 delay v5 半容差完整 PASS

第二档完成20次真实solve，worker469.737s、监督475.190s，最大single52.593s。独立636源/487artifact及监督SHA、strong receipt、原NPZ全局局部重算PASS。粗分支第9个真实probe得到22.5ns括区、1.0024nm位置差、25.184417842ns自身delay≤50，prefix188，post全33720边59queries获证并恢复1334–1371；细2第6probe45ns括区、2.0042nm、delay20.325618611ns≤50，prefix140，post53queries获证。前置/后置prefix除walltime外完全一致。旧第8高端60.086ns delay仍保存为未选探针，未续接或改判。细末态71.85107421874992ms/483，SHA8b7851bddefc4e92d3c7ccb641289fe25973e19b87fa41669c7647aa8a5ecb0c。quarter于20:08:42自动启动PID62220，尚待最后一档/跨档比较/普通续跑。


## 20:16 三档敏感性全部通过并恢复普通续跑

quarter19真实solve，worker455.825132s，636source+467artifact SHA/strong receipt/原global与local从NPZ重算全部PASS。coarse8probe、19.675603ns括区、0.8765nm位置差、own delay8.425719505ns≤25、prefix196、post65queries；fine2 6probe、22.5ns、1.0021nm、delay7.820612636ns≤25、prefix156、post59queries。两路后侧实际33720边图恢复1334–1371。

独立从三档真实NPZ及事件证据重算compare_group，与sensitivity_comparison.json完全相等：三组配对、各三角色共9路原endpoint screens及6项事件比较全部PASS，事件括区均重叠，gap=0。三档共同父态SHA ae6221bd1a785d242ce3a4cde52b5c5c1b3c4685ee3afeddaef821edcce2dcce，物理间隔只计一次。选quarter末态SHA a6bc6d81c5cba08f5e825c76a5d4400219b47d0ad0a40b788aea749cd78d559e，71.85107421874992ms/483；selected_physical_interval_count=1、adaptive passes=1。新料穿孔4、1.35mg、肩下0，未穿孔PID1333距出口32.263994µm。

20:16:34已在同一监督60096、同一历史活锁池内自动开始ordinary attempt_0001（child48708），确切父态为quarter；未重复startup或将另两档累加为时间推进。原截止01:37:12保持。首个普通段完整验收仍待完成；此阶段仍为工程近似诊断输出，非生产接受或冶金结合证明。


## 20:19 首个普通续跑完整通过——本次修复验收完成

ordinary attempt_0001确切父态为选定quarter，细步483→485，物理时间71.85107421874992→72.32568359374991ms，dt0.474609375ms只加一次。3真实solve、分支1/1/1、0事件；worker136.288038s/监督139.458101s，分支74.625419/33.417742/18.764105s，single最大71.494s。全636源/97artifact及监督manifest/log SHA匹配；strong receipt通过，NPZ重算原global/local与保存完全一致且PASS。三路可见路径17/11/11queries获证，union router均无事件。末checkpoint SHA142803c330fb7015fd199bde43c0a64fdf07004d5442abc37e33f01e76ecfc50。

保存态observer一致：4粒新料穿孔（1151/1203/1229/1281），总1.35mg，肩下0、留存0；最近未穿孔PID1333距出口32.074322502µm。新事件分段积分、连续几何证明、延迟/前缀共同收敛和三容差比较已实际跨过原失败区间；普通推进也已恢复。20:19后台自动启动attempt_0002，监督60096/子34088；原绝对截止2026-10-11 01:37:12.059551+08保持。所有绑定版本/组件收据冻结，旧失败保留；无需再启动另一任务。输出仍diagnostic_only、accepted=false、restart_allowed=false、正式接受增量0；未将这次求解器修复等同于已形成肩下沉积或冶金结合。


## 2026-10-10 22:59 状态核查：完整PASS之后因成本预测停止

20:48:31以STOP_COST_REASSESSMENT正常结束，exit0。7个普通段全部原pair_screen PASS、supervision exit0且manifest/log SHA匹配，无物理FAIL，未启动attempt8。第7段15次物理解/2支持切换；分支总576.017074s、worker585.551811s、监督590.505583s，完整pair取最大590.505583s。下一段经验预测30+1.5*590.505583=915.758374s，大于900s门槛15.758374s（1.75093%），因此停止。此前最大475.190269s，对应预测742.785403s。停止时仍余17320.898秒授权时间，原截止01:37:12尚未到。

独立载入最新完整checkpoint，实际1381粒、75.17333984374985ms/499步，SHA b1f950f4f5b947ed8d13f2b256f33b7ee7fa9f50b8dbbc5ba2f9b8541d8e9c67匹配。最新global stress_relative0.00509938，原global/local均PASS。4新料穿孔1.35mg、肩下0，最近未穿孔PID1307距出口31.755169438µm（最前ID已变化，不应当作同一粒子速度）。最终历史审计final_audits1、failedfalse、24355expected/locked一致；closed=false是close前快照，后续进程已正常退出，当前无python运行。本轮仅回答状态并更正交接记录，未改冻结源或重新启动。后续续跑须处理成本调度策略，保留全部物理门槛与原授权截止。

## 2026-10-10 23:55 用户确认90分钟，行政续跑入口已启动

用户明确确认“一个半小时（90分钟）”。授权窗口固定为2026-10-10 23:44:34至2026-10-11 01:14:34+08，准备与初始核验计入预算；此次较早截止优先于旧01:37:12。新入口只改变行政成本边界：forecast1200、physical soft1260、worker hard1290、controller outer1380秒。保留30+1.5倍全部历史完整pair最大耗时的预测、实际single-solve记账和剩余窗口准入。原915.758374秒预测可在新边界内继续；未删除高成本记录或放宽物理检查。

统一32项测试及独立源码/实际历史审查通过，组件收据`outputs/gpu_acceleration_20261007/event_resume_90m_components_20261010_v1.json`已封存，新5个行政文件冻结。物理源映射仍为原636项，generation6、.25事件容差、步长策略与沉积观测历史保持一致。真实旧policy和20帧observer重建与保存值一致，414项末段证据及继承2054项绑定核对通过；旧三档敏感性只读核验和重算比较，不重新执行GPU实验。合同见`GPU_EVENT_RESUME_90M_20261010.md`。

23:54:59新`event_resume_90m_20261010_v1`监督29104实际启动，bootstrap为PASS_DETACHED_JOB_POLICY、kill_on_job_close=false，独立读回进程存在且CPU累增。当前状态INITIAL_FULL_HISTORY_VERIFICATION，尚无新GPU子任务/新PASS。待续父态为旧event_delay_campaign_20261010_v1/attempt_0007，checkpoint SHA b1f950f4f5b947ed8d13f2b256f33b7ee7fa9f50b8dbbc5ba2f9b8541d8e9c67，物理75.17333984374985ms、499步；4新料穿孔1.35mg、肩下0。进度明确initial_parent_verified_in_this_session=false，初始核验结束后才会启动下一段。此时仅确认后台续跑程序已启动，不能声称新增物理推进或沉积。

## 2026-10-11 00:34 续跑9段通过，仍未新增穿孔或肩下沉积

实际监督29104/子64904均存活，progress为RUNNING_GPU_CHILD，9段完成、attempt10进行中，截止仍01:14:34。初始全历史preflight855.151009s通过，session protocol于00:12:28写入；截至00:33:46，9段从75.17333984374985ms推进至79.44482421874977ms（增4.271484375ms），499→517细步，dt始终0.474609375ms，无缩步/失败。近期每段连同核验约2分钟，较长的一段在求解结束后的核验耗时更多；不能把物理时间毫秒和实际运行分钟混用。

最新attempt9实际NPZ共1381粒子，SHA f1b1deb581cdb557581e7928cfa7622496660b880dcf7d8183d8c60dee353b20，与manifest/progress一致。独立有界核对636物理source映射与旧父7/本轮父8/实际NPZ implementation完全相等；global/local/crossing-ID检查通过，global stress_relative0.00497601782665，局部最大应力相对差约0.0397900424677，仍在原0.1门槛内。监督exit0、未中断、源输入匹配，manifest/log及13项关键artifact SHA一致；本段3次物理解、0支持事件，worker91.163s，监督95.091s。此次未重复全历史核验、未额外执行物理解。

saved observer与实际保存观察一致，38帧；新料穿孔仍4个（1151/1203/1229/1281）、1.35mg，肩下0、肩下留存为空。最近未穿孔仍PID1307，出口距离由本轮起点31.755169438µm增至33.295834697µm；这段距离变化不能支持单调趋零或给出沉积到来时间。已穿孔料最大半径0.542864mm，相对首次保存穿孔态最大径向增量约16.280µm；尚不足以构成肩下沉积证据。00:34 GPU瞬时62℃、显存4689MiB；当前保持原后台任务，未重启或延长授权。
