# Q7R 有界恢复合同（2026-10-02，运行前登记）

只贯通Q7已有J2/固定法向双切向摩擦基准，不添加物理、不替换生产默认。
父归档PROCESS_PLASTIC_FRICTION_Q7_FREEZE.json在修改前1429项完整性PASS。
既有src/scripts/configs/tests/outputs/logs全部保留；旧HANDOFF和两个索引文档可更新。

唯一算法预算改变：新增命名变体Q7R_BOUNDED_NEWTON14，将每个分支Newton
maximum_iterations从12改为14。依据是上一轮独立、未提交的固定12/14诊断：
第13次root=1.637588335e-10，但ID21/feature1强制stick超圆5.691232926 N。
原12失败不改。不得进一步递增预算试到通过。

## 原样保留的验收

分支相对残差1e-8；逐接触力9e-5 N；滑移1e-11 m；功账1e-6 J；
raw摩擦耗散低于-1e-10 J拒绝（舍入尾原样保存）；事件时间1e-12 s、
二分40、每请求区间最多12事件；线搜索8；condK<=1e12、condM<=1e6。
材料仍常Y=90MPa制造J2，不是AA6061/FKS标定；E/质量/时间/几何/边界不变。
Ws-Es-P超过1e-6 J仍告警，不补热、不宣称能量认证。

实际路径：advance -> probe/solve_cone_modes -> branch(target) -> feasible_cone
-> 必要时二分事件 -> PlasticFrictionBlock.commit重新计算全部材料/接触检查
-> accepted输出。root收敛不可绕过event或commit。任何新失败即保留候选并停止，
本轮不加第二种接触修复，不把诊断checkpoint改为生产restart。

## 先跑的回归与独立判据

1. 冻结.05s接受态上原12次重现REJECTED；试探不污染任何已接受状态。
2. 固定14次得到分支根；用独立norm(tau)-mu*N检查ID21超限>5 N，
   并核对记录5.691232926 N（仅归档重放，不是物理解析认证）。
3. 实际commit拒绝该候选，材料ep/P/R/Q、ID和所有接受态指纹不变。
4. 新旧构造器除明确预算及其元数据外一致；显式超过14或非法预算拒绝。
5. 复跑Q7独立解析/事务/标签/读回21测试；Q6旧兼容另列。

## 成本和运行授权范围

27粒子、9接触、54自由DOF；先4每相增量，共16请求区间＋事件。
上一轮预压6接受步约6s，Q6整段约16s；Q7R粗估20–50s，不保证收敛。
每个子进程外部timeout=60s，单线程。不是显式动力学，无CFL/dt稳定性认证。
若4增量完整到同一制造终时.2s、独立读回通过且预算可负担，才跑8增量；
否则8为NOT_RUN。空间n5、旋转工具、热、供料、动态AFSD、生产均NOT_RUN。
不得比较不同失败终时来声称收敛。保留最后accepted与REJECTED分离。

## 最小接口和输出

新增probe_process_plastic_friction_recovery.py中的RecoveryBlock只控制预算。
命令入口临时、显式向原main注入该构造器（finally恢复），复用完整输出路径，
不复制物理或改旧源码；只能单进程使用。旧provenance/report的source_sha256
保留父实现16项，新增recovery_provenance.json补充实际入口/变体配置hash，
并由report.artifact_sha256覆盖。resolved_config.recovery_variant显式记录覆盖。
新增独立reader先核验补充来源、配置差异再复用原数据/功账检查，不能只读旧reader
而忽略新入口来源。accepted journal/checkpoint仍不是生产restart。
