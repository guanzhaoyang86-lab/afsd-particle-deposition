# 有限patch补证后的三个有界肩下窗口（运行前）

独立目录finite_shoulder_windows/window_1、window_2、window_3分别仅登记22→30、
30→38、38→46ms。CLI为--window 1|2|3（默认1）、--stage coarse|fine。
每窗coarse32×0.25ms、fine64×0.125ms，单次独立预算不变；不自动无限窗口或减步。
首次父为22ms/153步SHA5368f67e717b4ce42d3927b8eca4d03a502a7ab59c5899d3d9b0b776c71d2d95。
后窗仅当前窗fine完整PASS_WINDOW_NOT_TRANSPORT且严格crossing_ids为空、质量0时可启动。
父为该fine最后接受checkpoint；累计prefix长度随真实接受步增加，不重设为153或清零。

旧shoulder_window_22_30ms的coarse在step019几何证明UNRESOLVED后失败，所有候选
均未接受。该旧运行及其源不修改、不改判，其18/19两个NOT_RESTARTABLE文件只用于
新几何证书与独立诊断的同输入资格，不能当作新窗口初态。

新contact_shell_finite_enclosure先保留旧256次几何查询结果，再用真实曲线的
坐标box保守包络证明旧未决有限patch区间，额外解析预算256，无新增有限patch查询。
新脚本/测试/正式receipt、同一旧失败路径的CERTIFIED只读证书、独立step019诊断
及其输入hash全部绑定。采样未见接触不是证书，FOUND/实体/负耗散等原拒绝不变。
数学范围见FINITE_PATCH_CHORD_PROOF_20261004.md；精确segment-to-box优化仅为未来可选，
本注册使用正式通过的axis-box证书，不在运行中换方法。

## 同一窗口内两条路径与提交顺序

复用旧已审查candidate/search/local-screen helper；新worker以静态源码明确替换
登记起止时刻与prefix长度，显式注入新finite enclosure审计。力学、材料、摩擦、
热学、单步位移/穿透/温域/密度、fresh接受和全部能量前缀门槛不变。
coarse全部为NOT_RESTARTABLE数值参考，可在内存继续，不产生accepted文件；只有
32步完整结束才为NUMERICALLY_QUALIFIED_COARSE_REFERENCE_NOT_ACCEPTED。

fine从本窗同一父态独立推进；每两half与coarse相同终时比较。这是同一窗口起点
的两条路径对照，第2对以后不声称相同直接父态。每对先完整计算两个候选，继承
本窗起点全部累计记录，再作全局位置5um/T.5K/S10%/ep20%筛选，以及累计所有SHOULDER
ID在该coarse和两个half真实NNPS图一跳邻域的逐点同门槛；S底限1MPa、ep底限1e-5。
历史SHOULDER IDs跨本系列窗口保留，不因当前脱离接触而遗忘，不递归扩张邻域。
只允许已验证SHOULDER及BORE；TOP/OUTER/BASE、新隐藏壳接触、实体、并列、吸引法向、
负耗散、热学或能量失败均停，不提交该对。

每个完整候选保存后，对其实际 endpoint position/i/j 原始行顺序和真实工具位姿执行
audit_solid_edge_visibility 纯几何守卫。边与工具严格实体内部相交的最长连续区间
超过1nm即保存候选和带PID的solid_edge_visibility.json并停，不发布；无效输入亦停。
0.1nm内缩实体仅作边界敏感性诊断，不能屏蔽严格实体分支。该守卫不筛除NNPS边，
不修改力或热学；只证当前端点直线边，不能声称证明两时刻之间的扫掠图或自由表面
核一致性。11个纯几何fixture、22ms实际7899边无遮挡读回、7次成本及输入/源hash
由solid_edge_visibility_receipt.json绑定，每次新候选仍须实际重算守卫。

每screen保存两支strict crossing IDs、质量、事件时刻；IDs不同表示穿越时间未分辨，
该对未提交并停。首次eligible中心穿越在两支IDs一致且全局/局部screen全部PASS后，
才提交当前两个half并立即终止整个系列，状态PASS_FIRST_CENTRE_PASSAGE_TIME_SCREEN_NOT_DEPOSITION。
coarse仍可完成已登记参考窗；fine首cross合法提前停止，不强求窗口终时、更不称沉积。

仅screen全部通过才依直接父链原子发布两个fine checkpoint；screen、同终时coarse
hash写入接受状态/记录。失败或中断保留旧合格前缀与所有candidate；发布中断可能
留下已经整对筛选通过的第一half。完整64步/32screen且到达本窗终时才为PASS_WINDOW_NOT_TRANSPORT。

## 后窗资格与硬边界

后窗逐条核前窗64个接受checkpoint的直接父链、来源/参数/时间、递增能量prefix、
32个全局/局部/strict crossing screen、对应coarse参考hash；同时核每步真实候选
与接受字段一致、连续壳/scalar等价/几何/端点边可见性资格文件的hash和PASS状态。沿登记窗口递归
核验至原22ms父，不能从未筛选的coarse路径或中途失败前缀跨窗启动。

每窗预测按（此前已测最大完整根成本+实测最大端点守卫成本）×32或64×1.5+30s，
fine还纳入当前coarse成本；若旧根已含守卫，仍额外加一次此成本作保守预测。
>450s不启动；内部540s/外部600s硬停，每根检查剩余预算。
Broyden24更新/100完整残差/2刷新，失败仅一次同父同初值原Newton24/100回退。
最多46ms、首次合格中心穿越或任一失败即停止，无新GPU、三工况生产或未披露缩放。
