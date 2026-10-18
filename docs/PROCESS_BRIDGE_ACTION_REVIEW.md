# 真实存档方向对现有数值修正的离线检验

2026-10-04。ALGEBRA_ONLY，不是新物理步、完整Jacobian或大步认证。
范围先登记于PROCESS_BRIDGE_ACTION_CONTRACT；使用650.393633ns父态、50ns dt、
真实保存的C+/map+和C-/map-，没有新evaluate、kick、材料/热更新或epsilon扫描。

## 实际执行与工具修复

初次check_saved_action在任何矩阵构造前因checkpoint state含ndarray而JSON序列化
失败，外0.3924092s，models=[]。原源码、FAIL和日志保存在
logs/process_bridge_action_20261003。没有物理解算失败或部分矩阵结果被丢弃。

事前单独登记PROCESS_BRIDGE_ACTION_REPAIR1，只将状态/参数JSON比较交给脚本已有的
递归plain转换，保留全部字段；复制到新目录，不改初次文件。唯一修复执行：

    .envs/pysph/python.exe -B logs/process_bridge_action_repair1_20261004/run_action.py

exit0，外2.6388493s，小于60s；status=COMPLETE_ALGEBRA_CHECK_NO_PHYSICAL_ACCEPTANCE。
两矩阵各对真实正/负候选apply一次，父态、参数、ID和来源不变；规定轴、proposal身份、
有限性及原线性后向误差界1e-12检查PASS。完整数据见同目录report.json/action_arrays.npz。
没有新accepted_history/checkpoint；accepted/restart_allowed/can_publish=false。

## 量的定义与观察结果

R=C-map；dC=C+-C-，dR=R+-R-。现有apply返回P(map-C)，因此
E=dC+correction+-correction-=dC-P*dR。比较同一有限扰动，而不是对连续导数求逆。
下面速度inf/L2比值均为norm(E)/norm(dC)，不是物理模拟误差百分比。

| 全域量 | 原结构P | 含连续性回代P |
|---|---:|---:|
| 速度inf比 | 0.357995607 | 0.366854414 |
| 速度L2比 | 0.0705854532 | 0.0680514209 |
| 最大速度E [m/s] | 1.39634413e-7 | 1.43089746e-7 |
| 正塑性粒子速度E [m/s] | 6.28567338e-8 | 5.45485260e-8 |
| 多面接触粒子速度E [m/s] | 1.44534268e-8 | 6.08914195e-9 |
| 有停滑事件粒子速度E [m/s] | 2.12737458e-8 | 1.11803302e-8 |
| 最大密度E [kg/m3] | 2.81120265e-8 | 6.53039027e-9 |
| 最大温度E [K] | 3.26053851e-10 | 3.26053851e-10 |

两者最大速度失配均在no_contact组，不能把它直接判作该粒子的摩擦切线错误。
无直接接触仍可通过邻域受到接触/塑性材料影响，不等于与接触完全无关。
连续性回代改善密度和多面组，但全域inf略增；这不是一个全面更好的完整逆。
温度和密度输入扰动均为零，故只报告交叉响应，不用零分母编造消减比。
温度E最大落在正塑性粒子；其物理源和BE账已独立重建，但本数值P没有完整热—材料切线。

## 结论与下一范围

CONFIRMED：现有P对这一个真实342点方向不是完整逆；它能够降低整体L2扰动，
但余项不只在接触点。两个冻结矩阵构造/代数解未污染真实状态，未改变任何受力。
SUPPORTED：继续仅围绕单面摩擦调节不能覆盖本方向的主要剩余误差。
POSSIBLE：非零初应力的客观率、局部塑性切线、当前几何变化以及热耦合均可能是
冻结构P漏项，必须依据实际源码和存档分解区分，不能由组名直接断言根因。
BLOCKED：完整真实342事件B/C不存在，StageEvent对照NOT_AVAILABLE；不能移植27点结果。
NOT_RUN：850ns父/500ns旧失败候选重算、新根迭代、新接受时间、真实输运。

本实验没有放宽raw根门，亦没有规定P必须满足任意新误差门。50ns单方向结果
不证明500ns收敛或时间准确度。应先定位一个有证据的遗漏搜索项并给独立基准，
再登记一次有预算的原失败大步回归；不续50ns刷时间、不反复调epsilon。
