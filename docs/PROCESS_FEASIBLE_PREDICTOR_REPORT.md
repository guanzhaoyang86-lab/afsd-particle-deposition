# 可行初猜：已验证最小修复，尚未提交新物理步（2026-10-02）

本轮只完成一项：旧起滑入口后的负法向初猜恢复。
新增companion solve_with_feasible_predictor；没有修改旧求解器/参数/失败记录。
工程近似PySPH基线commit仍1cd3dea9cfaac7d9f1bda2a2491ce4b42a31dbe5。

## 结果先行

独立解析回归13 PASS＋1旧guard RED；一次实际冻结目标对照中，
旧入口在ID19 N=-6.97211055N处停止，新入口1次预测修正后进入正法向域，
再4次原cone修正返回分支根，总5次，未超原14预算。
实际对照进程.6751822s，内部.2193021s，外timeout60s，exit2。
但ID19–26八个强制黏着接触已超摩擦圆，必须沿现有事件流程分段；
状态为ROOT_RETURNED_EVENT_SPLIT_REQUIRED_NOT_ACCEPTED，不是完整接触解。
accepted保持.1125s/16步，未增加物理进度；所有新候选标REJECTED且不可重启。

## 最小补丁与实际调用

原：Onset.advance -> advance_onset -> Recovery.branch -> ConeBlock.branch
-> trial（当前J2/成对corrected内力）-> state -> cone_response。
guard在第一个Newton残差之前看到负法向；此前只有接触点套用了目标位移。

新companion函数：原初猜 -> 原raw trial -> 若N<0，暂时冻结全部接触分量，
对非接触自由DOF用原tangent/reduced_correction做有界机械修正；
首个通过残差下降线搜索且全部N>=0的候选进入原state/cone_rows求根。
不是把任意N正的候选直接接受。预测和cone共用14修正，线搜索8不变；
支承/壁位移、释放掩码、材料历史、真实目标dt全部不变。
预测器不是物理子步，不累计塑性历史、不加入预测能量差/残差热。
实际commit/advance在诊断owner上均禁止，没有用修改accepted位移喂旧branch。

文件/接口：
- scripts/probe_process_feasible_predictor.py：纯分支函数，owner/t/released/
  maximum_iterations/observe；未接入生产默认。
- scripts/diagnose_process_feasible_predictor.py：精确旧/新各一次无提交比较，
  复用旧诊断恢复、旧reader/assess；保存每个原始trial而不更改旧文件。
- tests/review/test_process_feasible_predictor.py：独立线性消元回归。
- CONTRACT/REFERENCES：执行前预算/判据和独立推导；原始合同不改。

## 独立参考与实际测试

K=[[2,1,-1],[1,3,-1],[-1,-1,2]]N/m，x=-1m，
q为自由变量、s为规定壁切向位移。独立消元q=(s-1)/2,N=(3-s)/2。
s=2.5m：旧猜N=-1N；根q=.75m,N=.25N；Es=117/16J，
外功105/16J，摩擦功0。s=4m：根N=-.5N，必须拒绝。
mu=100仅为保持该弹簧例黏着的制造输入，不是AA/AFSD标定。
测试参考为独立常数/公式，不调用被测函数生成答案。
合同采用abs1e-11 m/N及1e-10J；实际实现比较用1e-11亦更严于功合同，
小良态线性参考，不是新的物理验收门槛。REFERENCES的1e-12为审查建议，
不是试后调门槛。N<0 guard仍无放宽。

| 检查 | 实际结果 | 日志（logs/process_feasible_predictor_20261002/） |
| --- | --- | --- |
| 独立初始10项＋旧guard反例 | 10 PASS＋1 RED，.366s | independent_initial.log |
| 最终独立含输入/法向边界检查 | 13 PASS＋1 RED，.367s；包含初始10项，不重复计数 | independent_final.log |
| 旧Onset兼容 | 27 PASS＋1原RED，7.273s | onset_compatibility.log |
| 旧Q7/Q7R兼容 | 39 PASS，23.925s | q7_compatibility.log |
| 旧Q6兼容 | 31 PASS，1.174s | q6_compatibility.log |
| 唯一实际目标对照 | 分支根返回，但须事件分割，exit2 | actual_comparison.log |
| 独立无求解器数据复核 | PASS_INDEPENDENT_READBACK_NOT_PHYSICAL_ACCEPTANCE | independent_readback.log |
| 完整.2s/8增量/空间/热/工具/AFSD/生产 | NOT_RUN | outputs/process_feasible_predictor_20261002/NOT_RUN.json |

还检查了预测＋cone共享预算、released接触在预测阶段不动但在cone阶段恢复自由、
正初猜路径与旧分支逐数组相等、无状态污染、零预算/奇异切线/非有限量明确失败。
N=0但自由残差不零的反例不能接受；N>0但摩擦圆超限也不能接受。
原本构/闭循环/空间失败没有重判或删除，其他历史套件本轮NOT_RUN。

## 实际数据：为什么可以说修复有效，但不能说目标步完成

冻结source/config/input校验通过；旧/新分别精确重建同一个accepted，
原材料、ID、壁程、release mask均相同。旧1次raw、新6次raw；没有额外目标搜索。
新raw最小N依次（N）：-6.97211055, 50.70480928, 13.55946842,
14.57922535, 14.59915777, 14.59925832。
预测阶段内自由残差.18531461 -> .10066083；后续cone恢复全部耦合求平衡。
最终分支root=2.9240111e-10（阈值1e-8），机械平衡3.7972574e-11。
仅已released的ID27接触根2.6093e-6N，在原9e-5N内。

然而强制黏着ID19–26的摩擦余量均为负，最严重ID23超限138.833853N，
不是舍入尺度；分支根不等于全接触模式正确。必须先定位各自起滑时刻，
不能把八个接触直接整段释放后提交来省去事件历史。
raw候选Q总量约1.6616717e-4J，最小分项-1.38727e-19J原样保留；
候选塑性功.00906984J，均未计入accepted。不能以Q非负证明事件可接受。

原accepted的Ws-Es-P=-.004346554773J仍为材料诊断警告。
输入输出源文件hash与accepted粒子/参考/数组/元数据完整指纹在前后全部相同。
全部候选/raw另存REJECTED_DIAGNOSTIC_NOT_RESTART，生产恢复防护未改。
独立后处理只导入NumPy/标准库，按ID定位后重算法向、摩擦圆余量、signed Q、
规定壁程、分支自然映射残差和4类accepted数组指纹，结果与诊断一致；
不调用被测branch/contact/material/assess。详见INDEPENDENT_AUDIT及原始复核日志。

## 证据分级与下一项

CONFIRMED：本例旧失败是可恢复的初猜定义域阻塞；不需要增加Newton预算。
CONFIRMED：新分支根存在且法向均正，但八个黏着模式超摩擦圆，必须事件处理。
RULED_OUT（仅本次停止点）：“初猜负N已经证明该分支需物理开离”。
BLOCKED：全接触可接受的.125端态、完整.2s路径，仍无证据。
NOT_RUN：接回现有事件/commit后的路径推进、正常开离再接触、AFSD生产。

下一项不是重写接触或本构：把已测纯分支函数接回现有Onset事件入口，
先验证原子提交与起滑定位，保留原负法向/事件/功账保护。
若有新正常脱离证据再单独评估；不能从某个预平衡失败自动推断完整cone无解。
完整路径执行前须重新估成本；本次.675s两分支不代表整段0.2s的成本。
真实AFSD仍2us/191步且未形成持续沉积，不报告完成度提升为70%。
命令见本轮COMMANDS.md，独立数据复核另见INDEPENDENT_AUDIT。
