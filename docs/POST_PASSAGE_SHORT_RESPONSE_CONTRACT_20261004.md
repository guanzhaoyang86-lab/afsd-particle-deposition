# 首次合格中心穿越后的单次4ms响应验证（运行前）

仅当refined_shoulder_continuation某登记窗口fine得到
PASS_FIRST_CENTRE_PASSAGE_TIME_SCREEN_NOT_DEPOSITION时可启动。CLI指定
--passage-window 1|2|3（默认1）与--stage coarse|fine；这只选择原系列的父出处，
本新目录post_passage_short_response只有一个窗口、不能追加窗口或再次减步。
RUNNING、FAIL、未完成pair、coarse穿越及普通完整无cross窗口均不是合法父。
原特殊PASS必须在首次非空strict crossing IDs的等终时pair通过全部screen后停止；
逐条核其接受链、真实候选、四种资格、全局/局部/双方crossIDs一致及全部累计ledger。
先按父的旧冻结来源核验，再冻结新driver/合同和该已完成且不再变化的父manifest/
histories/checkpoint；不绑定RUNNING结果或会更新的reader文件。全部历史FAIL保留。

新窗口从该真实父时刻t0到t0+4ms，32×0.125ms coarse、64×0.0625ms fine，
同父各自整窗推进。保留全部历史能量与被动observer状态，不清零crossIDs、速度、
应力、塑性、热、接触历史或累计记录，不跨接其他分支。

相同方程和全部有限工具/材料/热/根/接受门槛保持：Broyden24/100/2刷新加最多
一次同父同初值Newton24/100；fresh完整残差和所有能量账；连续壳256旧查询+
256解析补证；scalar-batch等价；BORE和已验证SHOULDER以外feature拒绝；实体、
并列、隐藏、吸引、负耗散失败拒绝。真实端点NNPS边严格实体内连续长度>1nm
拒绝，0.1nm内缩只是诊断，既不删边也不改变力学热学。

每两fine half先保存未提交候选，与coarse相同终时比较全局5um/0.5K/S10%/ep20%。
局部范围比原先更严格：累计SHOULDER事件IDs与累计strict crossing IDs的并集，
在coarse及两个half真实NNPS图上取一跳并集，逐点按原门槛比较；S底限1MPa、
ep底限1e-5。该集合只决定观察和拒绝，不改变力或热。双方累计crossIDs每对仍
必须一致，否则未提交停止；再次发生穿越不提前停止。本次只能走完4ms或遇任一
失败停止，不能自动细化、忽略事件时间不一致或增加窗口。

coarse仅数值合格NOT_RESTARTABLE参考，不发布接受。fine整对全局/局部/IDs
通过后才依直接父链原子发布两half，失败保留既有前缀和全部候选。只有64步/
32对完整到t0+4ms才为POST_PASSAGE_SHORT_RESPONSE_NOT_DEPOSITION。

每候选和screen额外输出父已cross粒子的实际ep及相对父增量、温度及增量、真实
工具局部位置/径向r/z及相对父变化、工具位姿和当前observer。非零塑性/径向位移/
升温不能自动宣称塑性铺展或沉积；其物理解释需结合轨迹、接触和能量证据。

按此前最大实测完整根（fine含本次coarse）加实测守卫成本，乘32或64×1.5+30s，
预测>450s不启动。每根剩余预算、内部540s、外部600s硬停保持。准备脚本或通过
入口核验不等于实际试跑通过；原fine尚未特殊PASS时本脚本必须拒绝运行。
