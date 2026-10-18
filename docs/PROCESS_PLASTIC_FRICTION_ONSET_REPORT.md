# Q7边界起滑入口：实际修复与部分路径推进（2026-10-02）

本轮实现命名变体Q7_BOUNDARY_ELIGIBILITY_V1，原先.05s处的起滑阻塞在实际
27点三维小块上已越过。接受态由6步/.05s增加至16步/.1125s，共7个内部事件；
已完成预压、正向切向加载并进入反向加载。但完整.2s程序仍 **BLOCKED**：
请求.125s时触发negative normal outside fixed-contact cone scope。
这不是AFSD沉积时间，不能称完整论文或连续体定量验证。

## 最小改动及依据

新脚本probe_process_plastic_friction_onset.py只改变事件入口：

1. 已接受接触在原力容限内靠近摩擦圆边界，且容量大于该容限时，
   获得自然映射求解资格。该资格不是指定滑动；方程仍可返回零滑移。
2. 同一整体耦合探测若明确向内卸载（margin>原ft），撤回资格，
   让随后反向到另一侧边界仍作为内部事件定位。
3. 零容量初始态不预释放；真正内部接触保持原strict二分。
4. 材料、整体法切切线、功账、14次分支预算和原子commit全部复用，
   接受条件不变；不裁剪力/负耗散，不改质量/模量/时间，不加入热源。

旧Q7/Q7R代码、配置、测试和失败输出均未修改。没有调整生产默认。
新onset_provenance记录新增入口/配置和恢复父来源；原16来源也保留。
新独立reader逐项核验精确配置差异、候选/接受trace与原逐粒子功账。

## 先验参考与实际测试

参考先登记于PROCESS_PLASTIC_FRICTION_ONSET_REFERENCES.md；容限、成本与限制
先登记于CONTRACT。对应测试/日志在tests/review及
logs/process_plastic_friction_onset_20261002/。

| 项目 | 实际结果 | 日志 |
| --- | --- | --- |
| 独立接触参考 | 10 PASS | independent_final.log |
| 原入口起滑反例 | 1 RED保留，expectedFailure，不计PASS | independent_old_red.log及最终合并 |
| 实际小块入口/零载荷/标签检查 | 3 PASS | actual_initial.log |
| 来源/候选trace/防篡改检查 | 14 PASS，合成数据不是物理解 | readback_tests_initial.log |
| 新测试合并 | 27 PASS＋1原RED，7.529s | onset_tests_final.log |
| 原Q7＋Q7R兼容 | 39 PASS，18.825s | q7_compatibility.log |
| 原Q6兼容 | 31 PASS，1.189s | q6_compatibility.log |
| 唯一4增量完整尝试 | BLOCKED，45.687s，非timeout | n3_i4.log |
| 实际输出独立读回 | PASS_PARTIAL_READBACK，exit2，.399s | n3_i4_readback.log |
| 8增量/空间/热/旋转工具/AFSD/生产 | NOT_RUN | n3_i8/NOT_RUN.json |

合并测试28项并不代表28PASS；focused与合并运行不重复累计。
其他历史Q4ep/Jaumann闭循环/空间FAIL保持原判，本轮未重跑全部历史套件。

关键独立参考：

- 单位切向弹簧初态在边界，正向壁移.2m：Q=.2J、材料不被强制拖动；
  反向小壁移-.2m：先stick，Q=0，壁功=弹性能变化=-.18J。
- 真正内部事件t=.5s；两接触整体耦合第二接触事件t=2/3，不是独立对角近似的1/2；
  总Q=2/5J，壁功44/75J。
- 共享自由二维节点的一滑一黏基准：Q=(.36,0)J、壁功.072J、材料功-.288J。
- 全反向先拒绝未经分段的相反力候选，再实际定位t=2/3，Q=1J；
  不以“Q非负”取代轨迹/分段检查。

容限带内正负初始偏差和横向转向另作离散端点参考，不冒充精确连续功。
例如-4e-5N初态偏差的单步中点Q与连续解相差约-4e-6J；
这些人为软弹簧偏差不是实际小块误差估计。任意曲线路径精度仍未验证。

## 实际小块结果与预算

不变：27粒子、9固定法向接触、54自由DOF，1cm立方体、常Y=90MPa制造J2，
E=70GPa、nu=0、rho=2700kg/m^3、mu=.1；不是AA6061/FKS标定。
每相4增量，计划16请求区间；实际10个请求（最后失败）、265个返回根的诊断候选、
16个接受区间、17帧。候选不是额外物理步。
新接受的内部事件约.06694718/.06964713/.07118869/.07205148/.11000767s，
原预压2事件仍保留。字段为实际三维应力/位移/塑性与接触账，不是动画驱动。

最后accepted .1125s：

- P=.251944605192J，Q=.014912570219J，max ep=.004181979456；
- Wwall=.317587396091J，Ws=.302674825878J，Es=.055076775459J；
- Ws-Es-P=-.004346554773J告警仍保留，未补热；不授予材料自由能认证；
- 最小法向反力212.144812N；接触根最大6.828055e-5N，小于原9e-5N；
- 历史raw摩擦增量最小约-6.435e-20J，作为舍入尾原样保存；
- 质量.0027kg，condM=3.385745，选定数组235764bytes，不是峰值RSS。

运行前估计20–50s，实际进程45.687s，外timeout60s；没有超600s或自动重试。
这是部分准静态程序的成本，不提供显式CFL稳定dt或论文尺度成本结论。
按部分结果不能声称完整程序保证60s内完成；下一次任何运行须重新估预算。
共同.2s尚未达到，8增量不运行、不用不同失败终时比较收敛。

## 新失败证据与下一阶段

attempt_0016_0.125.json显示：
.11250125s边界探测成功，随后TARGET .125s报负法向，iterations为空；
原REJECTED_manifest诚实标candidate_available=false，没有伪造完整候选。
已完成一次事前登记的纯诊断，外层.416s、exit0；恰好一次实际branch/trial，
accepted记录逐字段重建成功且全部状态指纹前后一致。未改生产restart防护。
捕获ID19/feature1初猜N=-6.97211055N，前一accepted为+212.144811664N；
九个规定法向gap仍全0，候选自由机械残差约.185315（不是cone root）。
实际调用上下文为原ConeBlock.branch第182行、iteration0，尚无Newton修正或线搜索。
这确认首个停止点是初猜超出正法向定义域，不是已求得平衡解开离；
加Newton次数对尚未开始的求根无效，但正法向目标解是否存在仍未知。
原raw trial、法向/切向、ID及约束账另存REJECTED_DIAGNOSTIC_NOT_RESTART，
不能当成最后有效物理态，也没有把负力裁零。详细结果见
PROCESS_PLASTIC_FRICTION_ONSET_FAILURE_DIAGNOSTIC_RESULT.md及
outputs/process_plastic_friction_onset_20261002/negative_normal_diagnostic/。

下一修复应先区分“固定法向求根的初始猜测不可行”与“平衡解需要正常脱离”：
前者只处理候选初始化/步内可行性，后者才需要有独立基准的单边法向主动集。
参考原probe_process_friction_cone.py的branch/state及cone_response，
保留本轮失败/新旧源码快照，未获证据前不放宽负法向条件或直接裁零。
下一独立基准建议（设计，**NOT_RUN**）：制造线弹性K=[[2,1,-1],[1,3,-1],[-1,-1,2]]
N/m，位移(x,t,q)单位m，q为自由DOF，N=-(K*u)_x。x=-1、t=0时q=-.5，
N=1.5N。目标t=2.5时保留旧q的初猜N=-1N，但自由平衡给q=.75、N=.25N；
目标t=4的自由平衡则q=1.5、N=-.5N，不能伪称存在受压固定接触解。
这个独立消元参考可同时捕获“假脱离初猜”和“固定接触确实不可行”，
不是对实际Q7目标解的证明。下一测试建议
tests/review/test_process_plastic_friction_feasible_predictor.py；
补丁仅限companion的branch初始化/有界局部延续，先隔离回归再考虑完整重跑。
本轮不增加第二种接触修复，不恢复真实AFSD生产；真实旧AFSD仍2us/191步无沉积。
重跑命令见logs/process_plastic_friction_onset_20261002/COMMANDS.md。
