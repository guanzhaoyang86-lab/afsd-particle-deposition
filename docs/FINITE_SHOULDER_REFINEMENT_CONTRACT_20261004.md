# 同父22→30ms整窗一次时间细化（运行前登记）

旧finite_shoulder_windows/window_1的0.25ms参考路径完整到30ms，但0.125ms路径
在25.5ms第14对局部应力screen失败：PID275相对差0.100688352超过原0.1。
全局应力/位置通过不能替代局部通过；旧接受26步至25.25ms与第27/28两个未提交
候选保留，不改判。本试验从原22ms/153步父重新独立推进，不能从25.25ms重新
计累计误差。共同父SHA5368f67e717b4ce42d3927b8eca4d03a502a7ab59c5899d3d9b0b776c71d2d95。

唯一登记窗口22→30ms：coarse64×0.125ms，fine128×0.0625ms。新目录
batched/finite_shoulder_refinement；CLI仅--stage coarse|fine。本注册只允许此次
一级细化，不能自动再次减步、换阈值、启动后续窗口或扩展三工况生产。

原有限窗口driver及其合同、源和全部几何/组件资格继续冻结。新preflight复用
finite.preflight(1,'coarse')取得原22ms父，逐文件核旧完整coarse参考、旧26步接受
直接父链和13对screen、28个候选四种审计、真实第14对失败及未提交27/28候选。
旧manifest/history、失败screen/候选和新driver/合同加入新冻结来源；此前状态
的implementation_source仍核其旧来源，不能要求旧checkpoint预先携带新脚本hash。

复用已冻结候选helper：相同Broyden24/100/2刷新与最多一次同父同初值Newton24/100
回退，fresh完整残差、材料/热/接触接受和累计能量prefix不变。新事务代码显式复制
已审的有限窗逻辑，仅调整64对和步长。每候选连续有限patch壳审计（256旧查询+
最多256解析补证）和scalar/batch对照必须通过；BORE与SHOULDER以外feature、并列、
吸引、实体、负耗散、隐蔽壳或能量失败均拒绝。真实端点NNPS边严格实体内部连续
相交超过1nm仍拒绝；0.1nm内缩仅诊断，守卫不筛除边、不修改力，仅证端点直线边。

coarse64步均为完整NOT_RESTARTABLE参考，不发布接受。fine每两half先保存完整
候选，与同一22ms分叉父态的coarse等终时比较：全局位置5um、温度0.5K、S10%、
ep20%，累计该窗SHOULDER事件ID与coarse/两half实际NNPS图一跳邻域逐点同门槛，
S底限1MPa、ep底限1e-5。两路径的累计shoulder seeds和观测均真实重算，不能沿用
旧失败路径预设接触事件；示踪与AS/RS只观察，不参与力或热。

每对保存两支strict crossing IDs/质量/时刻；IDs不一致停且两候选不提交。两支首次
eligible中心穿越一致且全局/局部通过后才按直接父链发布两个half，并立即终止为
PASS_FIRST_CENTRE_PASSAGE_TIME_SCREEN_NOT_DEPOSITION。没有中心穿越时，只有
128步/64对完整至30ms才为PASS_WINDOW_NOT_TRANSPORT，仍不能声称沉积完成。
中断可能留下已整对筛选通过的第一half，任何未通过的对不发布。完整153条旧能量
prefix原样继承；新参考/接受时间不能相加，也不合并旧失败路径的后22ms记录。

每次stage按旧两支实测最大完整根（coarse约1.46997s）与新coarse实测成本取最大，
另保守加已测guard成本，再乘64或128×1.5+30s；超过450s不启动。内部540s、外部
600s硬停及每根剩余预算检查保持。只有此次整窗结果决定后续是否值得继续；本脚本
不登记更小步长或新的窗口。实现与preflight检查不等于本物理试验已执行或通过。
