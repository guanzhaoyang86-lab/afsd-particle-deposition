# 工具遮挡失效前的单次2ms响应验证（运行前登记）

原post_passage_short_response的4ms coarse在39.875ms第20候选真实NNPS边
PID219–251穿工具132.760817nm后FAIL；其19条登记参考与第20完整拒绝候选原样
保留，不把原4ms改判PASS，也不证明或绕过工具遮挡失效。

本次只缩小被检验的物理范围到37.375→39.375ms：原coarse前16×0.125ms资格
逐条通过的NOT_RESTARTABLE候选作只读参考；从原合格首次穿越父独立求32×
0.0625ms fine，绝不从coarse或失败状态续算。共同父SHA
9f78a757bb9897e97598b09bee5842f29f6d6d7a2d217dfb39add0fb4a096b67，399条全历史
ledger原样继承，终态应431步。目录pre_occlusion_response/fine，CLI没有coarse
阶段、可调时长或更小步长，只有--preflight-only/--worker监督接口。

入口核原specialPASS父完整链；旧4ms manifest必须仍为FAILED_NO_NEW_ACCEPTED_STATE、
accepted0、登记候选19、真实NOT_RESTARTABLE文件20；逐条核19个候选旧来源/参数/
时间/递增ledger/直接数值父hash和四资格文件；第20必须仍是实际39.875ms遮挡
拒绝，证据绑定其完整候选/evidence。原19个候选hash、第20证据、旧FAIL manifest/
history及新driver/合同进入新冻结。旧coarse按旧postpass来源核，新fine按新冻结
来源保存；源码适配只改变验证范围和元数据，不改变力、材料、接触、热或空间图。

沿用postpass所有更严格门槛：每对全局5um/0.5K/S10%/ep20%，局部为累计SHOULDER
与累计crossIDs在coarse和两half实际图的一跳并集、S底限1MPa/ep底限1e-5；双方
strict crossing IDs必须一致；不因新增cross提前停。连续壳/scalar-batch/有限几何/
端点边可见性、fresh完整根和全prefix能量全部保持。超1nm实体边、其他非许可
feature、实体/隐藏/并列/吸引/负耗散及任何screen失败仍拒绝，不后处理删边。

两half先完整保存并比较，整对通过后才依直接父链原子发布。保留父已穿越粒子
的ep/T增量、工具局部r/z/位姿观测；非零响应不自动代表肩外铺展或沉积。只有
32步/16对完整到39.375ms才为PASS_PRE_OCCLUSION_2MS_RESPONSE_NOT_DEPOSITION。
失败保留已合格前缀与全部候选；无论结果如何，本注册不能超过39.375ms、不能
再次细化、不能自动开启额外窗口，更不能据此将39.875ms失效改判已修复。

按已有最大完整根加守卫成本×32×1.5+30s，预测>450s不启动；内部540s、外部600s
和每根预算保持。该限制窗可以检验首次穿越之后、空间模型失效之前是否出现
真实材料响应；其意义仅限此窗口和现有模型门槛。
