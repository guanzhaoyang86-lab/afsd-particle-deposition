# E3-1：工程近似AFSD耦合集成探针（2026-10-01）

## 用户目标调整与本轮范围（运行前冻结）

用户明确：重点是尽可能重构真实AFSD过程，不要求先达到弹性小块应力定量精度。
本轮允许带已知模型风险的APPROXIMATE短程探索；E2应力27.53%仍为FAIL，
不修改其阈值/报告，不再单独作为本次实现探针的否决条件。
旧Jaumann变体积能量及空间边界问题不因此消失；不得声称应力/温度定量可信。
此授权不包括长时、GPU、三工况生产，亦不允许放宽数值安全停止条件。

唯一活动后端PySPH 1cd3dea9cfaac7d9f1bda2a2491ce4b42a31dbe5。
使用已有run_afsd_debug.py/SolidSolver及显式event_aware_v1，不改物理核心、
材料参数、空间算子、热学、几何或接触容量。本轮不是再次实现替代本构。

## 固定问题与真实性边界

原configs/debug/starved.json及其extends，所有数值/单位/来源保持不变：
dx=.5mm，CFL=.1，层高1.5mm=3dx；657粒子(432基板+225方棒)，均可变形。
方棒2.5×2.5×4.5mm、基板6×6×1.5mm、肩半径2.7mm，旋转方孔开放。
300rpm、走速2.12mm/s、进给1.06mm/s。仅棒顶执行器x/z规定运动，y与其余粒子自由。
底座/工具刚性；材料旋转不得通过标签或规定自由速度获得。
预热初温723.15K为假设输入，工具/底座为恒温热库；FKS无量纲b1及低率floor显式近似。
杆/板初始共享SPH邻域，即预结合接口；没有预测接合生成/分离，也没有无限入口。
真实代码包含导热、塑性功、摩擦热、接触热，均按现有独立计账；未加残差热或裁剪。
几何是缩小热启动DEBUG，工具外圆/顶面边缘接触不完整，进入未支持工具实体即失败。

## 前置检查与运行预算

先验证ENGINEERING_E2_FREEZE.json、重跑run_coupling_review.py。
已知legacy/raw历史RED、当前材料自由能和节点边界RED允许保持FAIL并披露，
但接触类、反转、静/动墙、多面执行器、热学、ID/标签/排序/重启链路不得新增FAIL。
每项回归timeout60s；不是将整体suite的exit1改成PASS。

使用新实测preflight_timing.json的dt、首步耗时估191步成本；旧失败路径26步用31.57s，
旧新接触首步测量曾预测约472s。预计>600s不启动；本次内部上限540s，外部570s。
只启动一次固定参数尝试；遇失败/超时保存最后有效checkpoint与REJECTED候选，
不在本轮反复改变CFL/接触/材料来追求跑通。输出每个接受步，避免丢失最后接受状态。

拟用命令（通过当前前置检查与预算后才执行）：

    .envs/pysph/python.exe -B scripts/run_afsd_debug.py --case configs/debug/starved.json --output outputs/engineering_e3_20261001/dx05_cfl01 --contact-variant event_aware_v1 --end-time 2e-6 --steps 1000 --output-every 1 --wall-seconds 540

外部subprocess.run(timeout=570)确保连初始化/输出也在600s以内。
时间/质量/模量不缩放；没有隐藏阻尼或自由粒子事后位移。

## 先验验收与停止

1. 无新增前置FAIL；同一固定问题达到2us才记IMPLEMENTATION_INTERVAL_COMPLETE。
2. 原运行保护保持：质量相对1e-12、密度5%、穿透/dx10%、速度20m/s、
   能量残差5%交换尺度+1e-7J、FKS温域293.15至855.15K内、全部NaN失败。
   任一触发即停止，不以用户放宽精度为由越过。
3. 每个有效末态通过diagnostics再输出；失败候选REJECTED且不得作为物理末态/重启起点。
4. 检查ID、标签、源配置hash、实际工具位姿、自由粒子速度/位移、塑性和温度变化。
   没发生的塑性/迁移必须报告零；不能因为想看沉积就预设阈值外的轨迹。
5. COMPLETE只证明所请求时间段完成。2us对应工具转角约6.28e-5rad，
   棒顶供料约2.12nm，远不足以展示宏观旋转、沉积或混合，更非论文复现。
6. 即使短探针完成，本轮不自动增加终时或启动CFL/三工况序列。

实际调用已读源码：drift半步 -> corrected梯度 -> 材料update及中点应力内力 ->
event-aware接触+执行器反力 -> 导热/热源更新 -> kick及第二次半步drift ->
诊断/有效帧与checkpoint。半步接触局部子区间不等于整套耦合隐式求解。

结果、真实命令、是否达到共同终时、成本和下一任务在本文件后续结果节记录。

## 实测结果（单次尝试，未调参重跑）

前置run_coupling_review：116项111PASS/5原FAIL/0ERROR；新接触/全链路检查通过。
新试跑FAILED/exit1：有效21步、22个VTK帧，终时2.2029619823669545e-7s，
未达到请求2e-6s。失败候选22、particle_id476、BORE_X_MINUS+BORE_Y_PLUS。
不是E2应力阈值阻断，也不是本轮能量/温域/穿透报警；接触局部求解未返回。
本轮没有重新跑第二个真实算例，没有改CFL/几何/材料/接触容限。
相较历史旧接触实现26步，本次21步不能声称推进更远；代码版本不同，历史都保留。

实际计算用时46.1507s，物理步中位1.56431s，峰值工作集180146176字节(约171.8MiB)；
有效dt约1.04901e-8s。491.8s只是启动预测，本次因失败提前结束，非完成成本。
旧runner的paper_cost_projection是乐观线性外推且自身提示层厚分辨率不足，
不能用其数值声称论文尺度可运行，也不能据21步预测最终工程交付工期。

### 已接受末态真实观察

- 657粒子，其中632个非执行器粒子；源/配置/22帧/checkpoint哈希及ID/标签读回PASS。
- 自由粒子最大位移4.81060e-8m，速度最大.369322m/s；基板位移最大1.83513e-10m。
- ep最大0，塑化粒子0，塑性热0，进给面穿越质量0：尚未形成塑性流动或沉积。
- 温度范围723.14999689至723.15004441K，热量按已有导热/接触/摩擦计算；
  摩擦入材料热2.96513e-7J，接触净热-4.59993e-8J，不把符号裁正。
- 累计外墙功3.29395e-8J、执行器功2.01572e-9J、向墙净散热3.42512e-7J。
- 能量残差2.85492e-7J，归一.217104%；密度最大相对变化.008265%；
  穿透1.97718um=.395436%dx；末态alarms为空，质量.0002217375kg。
- 已接受接触子区间存在-2.45212e-24J舍入尺度原始摩擦热，未裁零；
  这不等于旧negative midpoint friction物理尺度失败已经重现。
- 工具转角6.92081e-6rad、平移4.67028e-10m、棒顶规定进给2.33514e-10m，
  都远未达到宏观沉积/混合所需时间。不得生成或解释为成形动画。

最后有效checkpoint.npz对应step21。REJECTED_candidate.npz明确禁止恢复，
新读回程序也只采用series/checkpoint的有效时间，而非候选或性能计数。
REJECTED_contact_debug.json原有辅助trace可能来自上一成功callback；
当前身份从failed_local_context读取，不用旧trace冒充当前区间。

### 局部失败证据与下一最小范围

被冻结的cone coefficient=6.579703290960622e-23；对应隐式中点区间
h=2*coefficient=1.3159406581921244e-22s，并非整体SPH步长。
一面的base slip约5.12e-16m/s，另一面约.0145552m/s，数值尺度差异极大。
报告全局scaled_residual=3.3492872871342685e-14、signed_gap=2.842170943040401e-14，
超过求解器内部2e-14目标；原调用容限仍为2e-12，未在本轮放宽。

CONFIRMED：真正停止原因是该双面近停滑子问题未取得内部证书；
不证明大尺度负热/穿透/温度失控，也不能仅凭小的全局证书误差说每面力准确。
新增test_contact476_open_risk.py直接调用实际求解器，FAIL已保存；无xfail。
独立高精度参考只用于此冻结局部问题，绝不回填REJECTED候选或自动继续时间推进。

下一任务只限contact_cone.py/接触事件层的近零滑移、多面尺度处理：
先使用冻结476问题与既有432类级/439/454/反转/功恒等式回归，
量化每面力、证书、冲量和末速度误差；不要只加迭代或直接扩大容限。
如需改异常记录，domain_kick目前漏传inner interval_context，可补诊断并单测，
但本轮未改任何src/，不能声称该记录已存在。
不把E2精度修复重新作为此局部实现问题的前置条件，不换本构或后端。
局部修复经独立验证后才另行决定新的单次耦合尝试；不自动恢复拒绝状态。

### 独立冻结子问题参考（已运行，不是修复）

scripts/probe_contact476_reference.py仅读取被冻结的双摩擦圆盘凸二次问题；
80位Decimal、独立推导KKT乘子二分/块最小化，不调用被测求解器。
5次sweep、实测.155394s，原始及参考数据在contact476_reference.json；
可行性、乘子非负、互补性及驻点证书全部通过，exit0/REFERENCE_CERTIFIED。
参考scaled gap约3.08e-68；这里的1e-60仅是高精度参考证书要求，
没有作为生产容限，也没有更改生产内部2e-14目标或调用方2e-12容限。

参考两面的二维切向力(N)分别约：
(+.3783809173,-.08710952287)、(+.3899002403,-2.621836e-15)。
候选力向量差.4586616523N，不能说“力误差可忽略”；
候选独立高精度scaled gap仍为2.84141128e-14，说明不只是证书求和的舍入误报。
但由于此次局部h约1.32e-22s，在冻结局部算子下，
由力差对应的末速度差仅1.78836e-16m/s，目标函数差1.38691e-16W。
该条件估计不覆盖完整接触事件、事件时间误差或后续轨迹，绝不授权接受REJECTED。

SUPPORTED：近停滑事件细分与多面数值尺度问题值得优先修复；
尚未证明应修改哪条事件分支，亦未证明只放宽容限可以安全前进。
本轮只锁定可独立复现的数值阻塞，实际耦合状态仍为FAILED_UNCHANGED。

### 可重跑入口与状态

- scripts/run_coupling_review.py --output logs/engineering_e3_20261001/gates：exit1，111PASS/5旧FAIL。
- 真运行完整命令见logs/engineering_e3_20261001/PREFLIGHT.md、trial_driver.log：exit1/FAILED。
- scripts/inspect_engineering_integration.py --run outputs/engineering_e3_20261001/dx05_cfl01 --output 新JSON：
  数据完整性PASS，但因未到终时exit1，不把读回成功写成仿真成功。
- tests/review/test_engineering_integration_e3.py：5PASS；仅有效/拒绝分离、身份/单位/诚实分类回归。
- tests/review/test_contact476_open_risk.py：1FAIL；新实测失败回归，不删旧失败。
- scripts/probe_contact476_reference.py：exit0，独立局部参考证书PASS；真实算例仍FAIL。
- 真实2us完成、CFL序列、空间细化、塑性沉积、三工况生产：NOT_RUN/NOT_ACHIEVED。
- 本轮只新增诊断/回归与计划文档，所有已有物理源/参数/上游未变。

