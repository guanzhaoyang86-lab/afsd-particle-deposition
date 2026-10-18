# 一次两次出生时间/预装有限长棒对照合同 V1

这是 `TWO_INSERTION_REFERENCE_PROPOSAL_20261005.md` 的有界实现。唯一允许的新物理工作为
一次968点完整map/tangent/rank成本候选、分期出生fine4步和预装长棒fine4步；总8根。
旧粗链保持原文件/分类/FAIL历史，仅读`two_insertion_diagnostic_20261005/v1`及独立读回。
全部新状态沿用独立热QS诊断模型，附study_scope=`TWO_BIRTH_TIME_AND_PRELOADED_REFERENCE_V1`，
accepted=false、restart_allowed=false；没有入口过程升级、没有继续窗口。

物理参数来自独立8mm initial_v2，feed/traverse/omega/G/rho/cp和标签作用不变，无任何缩放。
H=1.953125 us，终时4H=7.8125 us。分期出生t0和2H各25点，z3.25/3.75mm，T723.15K、
rho0、S/ep/初接触史0、初速(0,0,−feed)；midtime第二出生前的库存943，终态968。
固定world-z>=2.5mm且component1控制区仅z驱动；不移动控制面以制造释放。
预装参考从t0即拥有IDs919..968及同样初值，50点属于初始库存，source_events为空，
初始热库存明确比原918多6.530625J（Tref293.15K）；这不伪装成出生事件。

新固定库存source-work接口 `fixed_inventory_source_work.evaluate` 必须核parent/prepared/actual
持久ID集相同，且prepared旧数组未变；只在此依据下认定质量/动量/热/弹性对流输入0。
压力流0依据显式外部零牵引，反力功从实际力·规定轴位移重算，不能置零或重复加内部截面功。
完整material/contact/thermal库存差和raw代数必须与原实际map能量账一致；负反力功与raw缺陷
不裁剪，QS动能不补热。新接口不能接收birth事件或修改冻结birth evaluator。
6项独立制造算术检查覆盖热/弹性恒等、卸载负功、raw缺陷保留、排序、真出生拒绝、篡改拒绝；
这些是组件资格，不是实际热流或塑性验证。

各根沿用完整visible map与切线，fresh、scalar、实际图/修正矩阵/弹性参考rank、连续shell/path、
原接触扫掠/事件白名单、逐粒热账和全prefix机械能量均保留。弹性参考rank中固定归一化dt
只用于组装线性参考应变导数并精确抵消，不改变实际物理H；它不证明非线性切线稳定性。
每步V2外部verifier核8gate文件、全部audit SHA、实际PA/history/time/sourcework绑定，然后seal。
每步真实NPZ读回与reverse sort，下一步按ID恢复历史；最终birth observer不可退化成legacy。

时间对照：原粗第1/2步 vs 新分期fine第2/4步，同出生时刻、同PID库存，使用原global和local
场门槛（位置5um、温度0.5K、S相对0.1/分母floor1MPa、ep相对0.2/floor1e-5）。
local固定原150料棒ID，加粗细实际图一跳，并保存绝对S差/未加floor范数。
另要求fine全prefix raw绝对累计<=粗链累计+512eps×两支累计功幅度，禁止借放宽原allowance通过。
预装参考在2H/4H与分期fine比较原918 IDs，local原料棒及两支一跳与原918交集；使用相同场门槛。
所有额外点的差异另报告，不修改其真实早生/晚生运动史。分别报告完整反力功、源热/初始热库存、
raw/Wact及floor影响；控制反力和出口邻域统计为观察量，不捏造未发生的塑性或穿越。

CLI阶段固定 `cost_968`、`inlet_fine`、`preloaded_fine`、`compare`。前3阶段独立worker540/600s，
只有实测保守预测<=450s且至少2更新预算才允许运行根；两个fine阶段各4根，无fallback。
cost前参考943测量按DOF三次量级预测；cost后使用968实测map、tangent/linear、rank，以
`30+1.5*4*((min(100,1+8*M+16)+2)*e+max(M,3)*t+rank+10)`选择M<=24。
8根是总计划数，失败即停，不做替代路线或额外求根。compare为saved-only，无NNPS/物理求值。

即使全部通过，状态也只能是
`PASS_ELASTIC_ISOTHERMAL_TWO_BIRTH_TIME_AND_PRELOADED_REFERENCE_NOT_PROCESS_ACCEPTANCE`；
要按实际字段记录是否仍弹性等温，释放未出现写NOT_EXERCISED。它不注册持续入口、不推进接受
过程时间、不代表完整沉积；主线有限库存肩下铺展不依赖持续补料开发。
