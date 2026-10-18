# 已验证首肩之后22→30ms有界窗口（运行前）

共同父态为first_shoulder_enclosed/halves第二个已接受半步，22ms、累计153步，SHA为
5368f67e717b4ce42d3927b8eca4d03a502a7ab59c5899d3d9b0b776c71d2d95。
首肩全局及局部时间筛选须PASS。真实历史连续罚壳门槛由0→21.875ms的151段只读链审计
PASS，加其后两个实际半步各自CERTIFIED、同源与直接父链共同构成；不把151段声称为153段。
全部153条累计能量记录、材料/接触/温度/被动身份继承，不清零或重初始化。

## 两支及发布顺序

coarse固定32×0.25ms至30ms。每步完整求根、fresh接受门槛、全前缀能量、几何及标量
对照通过后，仅形成NOT_RESTARTABLE参考候选；允许候选内存状态继续生成下一参考，
不能作为可恢复接受历史。全部32步完成才为NUMERICALLY_QUALIFIED_COARSE_REFERENCE_NOT_ACCEPTED，
accepted_steps始终0、accepted_time保持22ms，候选数与候选时刻独立记录。

fine从同一22ms父态独立固定64×0.125ms。每两步到达对应coarse时刻，先保存两个
未提交候选，再对coarse作全局位置5um、温差.5K、应力10%、ep20%的原门槛筛选。
同时对当前及此前窗口内coarse/fine所有SHOULDER粒子ID并集及其在该coarse、half1、
half2真实NNPS图的一跳邻域做逐点同门槛；应力Frobenius分母地板1MPa、ep地板1e-5。
初始已验证肩接触ID也计入，不递归扩张邻域，不因一个粒子后来脱离接触而遗忘它。

只有一对全局及局部筛选全PASS，才依直接父hash链原子发布两个fine checkpoint；
对照失败时这两步都无accepted，已通过的旧前缀保留。每对保存原始screen及哈希，
每步真实事件ID/feature、完整候选、接触力/法切向/滑动/历史/释放能和能量账均保留。
仅完整64步、32个逐对screen通过且终时30ms时为PASS_WINDOW_NOT_TRANSPORT。
发布中断可能保留已整对通过筛选的第一半步；不发布未筛选的第一半步。

每个同终时screen同时保存两支严格bore_observation的crossing_ids、质量及事件时刻。
IDs不一致表示穿越时间尚未分辨，该对两步都不提交并停止；该标签只观察、不反馈方程。
coarse仍完成登记的30ms参考序列；fine若首次真正eligible中心穿越出现，只有两支IDs
一致且原全局/局部门槛均PASS，才提交当前两half并立即停止，以
PASS_FIRST_CENTRE_PASSAGE_TIME_SCREEN_NOT_DEPOSITION记录，不要求继续到30ms或宣称沉积。

## 不变的拒绝规则与成本

复用未改的Broyden24更新/100真残差/2刷新；失败仅允许一次同父同初值原Newton
24/100回退。各候选fresh完整验收、全前缀能量和原实体扫掠不变。
显式将原event.audit_candidate helper的审计函数替为已通过的continuous enclosure，
物理方程、材料/摩擦/热学、根求解器和接受门槛均不改变。每一步还执行原完整标量/
批量候选对照。只允许孔壁与已验证SHOULDER；TOP、OUTER、BASE、新隐藏壳接触、实体、
并列、吸引法向、负耗散、温域/密度/位移/穿透或任何能量失败均立即停止。

每支预测使用最大已测完整根成本×本支步数×1.5+30s，超过450s不启动；fine预测还
纳入coarse实际根成本。内部540s、外部600s硬停，每根检查剩余预算。
失败不自动减步、延长窗口或改变阈值。源、合同、父checkpoint、参数、旧链审计及
同源新半步证书在运行和提交前绑定核验；coarse参考链独立hash，不与fine时间拼接。
30ms或肩面流动也不自动等同于合格孔口穿越、沉积或实验验证。
