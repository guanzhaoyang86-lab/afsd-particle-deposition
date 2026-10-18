# 原细化等级下30ms之后的三个有界窗口（运行前）

本控制器只准备延续，不保证30ms细化试验会通过。第一父必须是
finite_shoulder_refinement/fine完整128步/64对通过的30ms最后接受态，连同先前
22ms的153步共281条累计能量记录。任何RUNNING、FAIL、中途前缀、首穿越提前
停止或非零cross均拒绝；不得跳接旧25.25ms/179步接受态，也不更改任何旧失败。
30ms父的具体SHA由完整接受历史与checkpoint在preflight中相互核验并写入manifest。

仅登记30→38、38→46、46→54ms三个8ms窗口，CLI --window 1|2|3（默认1）、
--stage coarse|fine。每窗coarse64×0.125ms、fine128×0.0625ms，保持此前唯一
细化等级，不登记更小步长或无限窗口。每次单独启动一支；第二/三窗只从对应
前窗完整PASS_WINDOW_NOT_TRANSPORT且cross IDs为空、质量0的fine最终态继续。
旧粗/细路径的物理时间不能拼接或累加。

新来源冻结含原refinement全部来源/旧失败gate，以及新driver/合同和30ms完整
manifest/history/父checkpoint。首父逐条按refinement原冻结来源核验；后窗按
本新冻结来源核验，不能要求旧checkpoint预先包含新源码hash。每条fine接受态
核直接父hash、分支起点hash、源/参数/时间/ID字段、完整递增能量prefix；逐对
核screen hash、真实coarse等终时参考hash与非零cross禁入。每候选连续壳、scalar
等价、有限几何、端点NNPS边可见性四份资格证据及文件hash全部核验。来源冻结
在三个窗口间固定，不把后续结果回填第一父来源。此前SHOULDER事件IDs跨窗口
累积保留，脱离接触不遗忘。

复用冻结finite candidate helper；显式事务代码仅替换注册窗口与64/128计数。
方程/材料/热/接触、Broyden24/100/2刷新、最多一次同父同初值Newton24/100回退、
fresh接受、单步/累计能量以及全部物理拒绝保持。连续有限patch证明仍为旧256
几何查询加最多256解析补证；BORE/已验证SHOULDER之外、并列、实体、隐藏壳、
吸引或负耗散均拒绝。实际端点NNPS边严格实体内连续相交>1nm仍拒绝；0.1nm
内缩只诊断，不删除边、不改变力或热，也不冒充两端之间的扫掠边证明。

coarse仅完整NOT_RESTARTABLE数值参考，可内存推进，禁止发布接受。fine每对
两half先保存全部候选，再与本窗共同分叉父的coarse等终时作全局5um/0.5K/
S10%/ep20%比较，以及累计SHOULDER IDs在coarse和两个half实际NNPS图一跳并集
的逐点同门槛；S底限1MPa、ep底限1e-5。只在全部通过后按直接父链依次原子发布
两个half。中断可能留下已整对筛选通过的第一half；未通过的对不能发布。

每screen保存双方strict crossing IDs/质量/时刻。IDs不同则未分辨、停止不提交。
首次真实eligible中心穿越只有双方IDs一致且全局/局部全部PASS才接受当前对，
立即终止全部系列为PASS_FIRST_CENTRE_PASSAGE_TIME_SCREEN_NOT_DEPOSITION，
不能再进入下窗或称沉积完成。无cross时须完整128步/64screen到注册终时才为
PASS_WINDOW_NOT_TRANSPORT。被动身份/AS-RS不参与力学热学。

每支启动成本按先前实测最大完整根（fine还包括同窗新coarse）加保守实测guard
成本，乘64或128×1.5+30s。预测>450s直接拒绝；不截短注册物理窗去绕过对照。
内部540s、外部600s及每根剩余预算守卫保持。最大到54ms、首次合格穿越或任一
失败即停止。此合同和preflight只提供运行资格，不能据此报告任何新增接受时间。
