# 停滑/反向事件开发结果（2026-10-03）

项目尚未完成，真实AFSD没有新增物理步或材料输运。本轮交付的是经过独立基准检验的
停滑组件、参考时序、事务修复与一次有完整失败证据的局部耦合推进，不是论文复现成功。
后端仍为PySPH 1cd3dea9cfaac7d9f1bda2a2491ce4b42a31dbe5。

## 本轮确实完成了什么

| 项目 | 实际状态 | 证据 |
|---|---|---|
| 原189项+新增44项回归 | 233 PASS，外14.095s | logs/process_stop_20261003/preflight_initial/report.json |
| 1D反向、可行粘着、完整矢量/单侧方向 | PASS；独立解析reference，不裁剪 | test_hht_zero_event.py / test_hht_event_reference.py |
| 严格零端点接触积分与私有候选隔离 | 8 PASS，含负节点热原值 | test_hht_stop_contact.py |
| checkpoint、诊断/积分绑定、重启物理连续性 | 4+10 PASS，保留初始FAIL与源码快照 | test_hht_stop_stage.py / test_hht_stop_transactions.py |
| 独立3D反向参考 | 105共同样本到6.5us，原9项误差界PASS | reference_endpoint2/report.json |
| 耦合3D局部基准 | 新接受5→6us；停滑根FAIL | outputs/process_stop_20261003/scout_initial/report.json |
| 输出完整性/独立功热读回 | PASS完整性；不等于精度PASS | PROCESS_STOP_OUTPUT_REVIEW.md |
| 耦合事件、一般3D粘着、多面执行器、正塑性接入 | 尚未验证/未实现；不继承组件PASS | PROCESS_STOP_CAUSAL_REVIEW.md |
| 真AFSD连续过程、三工况生产 | NOT_RUN；仍850.393632910904ns，crossings0 | 旧AFSD证据保留 |

新增代码位于src/afsd_backend/pysph/hht_zero_event.py、hht_stop_contact.py、
hht_stop_stage.py、hht_stop_flow.py；来源、接口、实际调用和边界限制见
PROCESS_STOP_IMPLEMENTATION.md。没有更换本构、空间算子、物性或PySPH设施。

## 独立参考：确实跨过了首次反向，但没有解决全部接触

初始将PID4/22作为同时事件的方案失败：均值零掩盖了各自约±5.25e-13m/s的非零值。
原源码/失败候选保留。第一次数值修复按实测时间逐个定位并切换力，不投影速度。
第二次修复只排除已定位端点的舍入别名：完整矢量仍正向增长，只是极小横向分量过零；
用多项式及单侧导数证明，原判零界不放宽，原始滑移/功热不修改。

最终首次两个事件的粗细时间差分别4.236e-16s、4.350e-16s，均在预注册1e-12s界内。
旧参考前缀逐位保留；事件前后全部物理/能量账状态逐位不变。
共同认证终时6.5us，增加9个反向后采样。其后粗网格在约6.5657616us、细网格在
6.5564031us遇到可行静摩擦：所需力约.494/.568N，容量约2.705N，当前一般3D粘着
尚不支持。两个失败终时不同，不能宣称该后续事件收敛，也不能授权8us。
详细来源、测试、两轮修复和初始FAIL见PROCESS_STOP_REFERENCE_RESULTS.md。

## 实际耦合调用与第一次新阻塞

仅运行一次scout，从旧发布索引绑定的5us父态出发，绝不读旧REJECTED6us作为父态。
预测最多17次联合求解、45.852s；内部55s/外部60s。实际外4.576s、内部4.112s，
峰值工作集225767424B。发布目录仅两个记录：同字段定义迁移5us、一个真实接受区间6us。
5→6us新解用了14次根评估、核心.9041s；raw Q=4.72611071448e-6J，
速度残差3.43363805e-10m/s。新的来源/参考证据允许这次发布；不反向改写旧拒绝记录。

6→6.5us触发自身dense路径的停滑预测。私有duration试探3/4各9次评估收敛，
但仍是不可发布/不可重启的延拓候选。第三个探根区间dt=2.7254762326e-7s在18次
评估后回溯失败：PID22横向系数超自然共线界约5.08%，同时全局速度残差
3.83375476e-6m/s尚未达到1e-9m/s门。残差峰值在其他角点的y自由度，
不是PID22自身；不能只凭函数名或小速度就宣布唯一原因。
所有候选、符号热、残差和回溯日志保留。没有接受6.27/6.30us状态或重启事件。

## 已知精度与成本风险没有被隐藏

6us同时间独立比较：位置RMS约1.952e-8m，速度RMS约.041851m/s，完整应力RMS
约.330632MPa，最大温差约3.30e-7K。PID4的滑移z仍−.0226522m/s，而参考已
变成+.000559327m/s：粗时间积分的分支延迟真实存在。这里只报告RAW差异，不称收敛。
输出读回/能量记账PASS不能抵消这个物理时序误差。

单个接受1us的核心求解成本约.904s；简单线性推到1ms已约904s，尚不含I/O，
也不是AFSD成本预测。事件和精度问题未解决，禁止据此启动长程或生产。
本轮S1/S2、时间/空间收敛、8us/1ms、事件后checkpoint重算均NOT_RUN。

## 下一步不是继续盲重跑

依据实际源码和失败，先独立核验完整接触冲量修正对候选加速度的响应，包含真实dense
路径积分、端点修正、法向gap及PRIVATE延拓的符号；只加endpoint Coulomb导数不够。
验证后才决定一次saved trial5私有根诊断。接触导数遗漏目前只是SUPPORTED候选，
不是已证唯一根因。顺序事件、静摩擦及长期时序精度还需单独验证。
细节见PROCESS_STOP_CAUSAL_REVIEW.md；旧FAIL不删、界不松、自由粒子不被驱动到期望轨迹。

可重跑命令在logs/process_stop_20261003/COMMANDS.md。PROCESS_STOP_FREEZE仅证明
所列文件完整性，不证明物理正确、真实材料输运或项目结束。当前层级仍是
ENGINEERING_APPROXIMATE方法开发与局部验证，尚未达到完整AFSD工程过程交付。
