# E3-3：真实0.5us短窗口已完成（2026-10-01）

## 结论

SHORT_WINDOW_COMPLETE / ENGINEERING_APPROXIMATE / NOT_PAPER_VALIDATION。
本轮只运行一次原物理问题，48个接受步恰到5e-7s，49帧；exit0。
用时122.2445721s，峰值工作集120254464字节(114.6836MiB)。
没有改src/、上游或物理配置，没有拒绝状态/超时/报警。
超过旧E3的21步/.220296us接受终点；不能声称重新经历完全相同的476事件。
塑性应变、塑化粒子、塑性热、供料面穿越质量仍为0，尚无沉积。
原2us目标、时间/空间收敛、三工况生产仍未完成。

## 运行前条件及可重跑命令

SHORT_WINDOW_E3_CONTRACT.md先于真实运行创建，范围/成本/停止条件未事后修改。
CONTACT476_SCALE_FREEZE.json验证830文件PASS；quick gates为3PASS(含64组)+1PASS。
旧116项111PASS/5FAIL及432组在源码未变情况下引用冻结证据，
本轮没有重复运行整套，也没有把旧FAIL改判。

    .envs/pysph/python.exe -B scripts/run_afsd_debug.py --case configs/debug/starved.json --output 新目录 --contact-variant event_aware_v1 --end-time 5e-7 --steps 80 --output-every 1 --wall-seconds 300

本次输出outputs/short_window_e3_20261001/dx05_cfl01_05us；
真实外层subprocess timeout330s，完整PowerShell/环境见logs/short_window_e3_20261001/COMMANDS.md。
从原初态开始，没有恢复旧REJECTED或旧源码checkpoint。
数据中的APPROXIMATE分类未更改；CLI的0.5us为显式实现窗口，
不是把较短物理时间伪装成2us完成。

## 真实状态与账目

| 量 | 实测 |
|---|---|
| 物理几何/分辨率 | 原dx=.5mm、CFL=.1、657粒子、层厚3dx不变 |
| 全时段最大密度相对变化 | .0157195% |
| 全时段最大穿透/dx | .400024% |
| 全时段最大速度 | .404609m/s |
| 最终最大自由粒子位移 | 9.48769e-8m=94.88nm |
| 最终基板最大位移 | 8.36601e-9m |
| 最终温度范围 | 723.149992288至723.150077391K |
| 最终/最大账面能量残差 | 5.9157883e-7J；当前归一尺度.4498698% |
| 最终外墙功/执行器功 | 6.6573545e-8 / 7.7665988e-9J |
| 最终向墙净热输出 | 7.4219972e-7J |
| 最终摩擦入材料热/接触净热 | 5.6489688e-7 / -1.7730284e-7J |
| 全时段最小原始摩擦子区间热 | -1.3925551e-23J，舍入尺度保留，未裁零 |
| 最大接触局部动能-冲量功残差 | 5.2409736e-24J |
| 塑性热、ep、穿越质量、新注入质量 | 全部0 |

质量始终.0002217375kg；所有49帧有限值、身份/被动标签保持。
能量账的独立重算与已存记录一致，但elastic_shear仍是当前Jaumann近似诊断，
旧闭循环非零功/E2应力/空间边界问题未因此被消除；不是严格自由能认证。
原方孔2um干涉/罚接触储能可主导启动瞬态，不能把运动等同真实送料铺展。

实际工具转角1.5707963e-5rad；走刀1.06nm，执行器轴向供料.53nm。
远不足以表现宏观旋转、塑性挤出、沉积和示踪混合。
没有生成动画，未以图形外观判定成功。

## 时序、性能与局限

实际顺序未变：上游半漂移 -> 材料预测/中点应力内力 ->
接触和执行器耦合kick -> 材料状态提交 -> 导热/热源 ->
上游kick/半漂移完成 -> 诊断 -> 接受输出/checkpoint。
不是接触后再强制覆盖自由材料轨迹。

前47步dt范围1.04900139192e-8至1.04906949246e-8s；
最后6.96108561504e-9s只为精确落在5e-7s，不能用它当稳定CFL步长。
物理步耗时中位2.2268681s，端到端平均2.5467619s/步。
48步累计305244个接触接受子区间、3791个内部拒绝细分；
这些是局部事件积分成本，不是隐含质量/时间缩放。

旧runner的paper_cost_projection把最小dt用于外推，含本次截短末步，
且其线性N假设/名义层厚未满足验证；不能用该输出预测论文尺度可行性。
本次完成不产生一般后续成本保证，尤其后续粘滑/塑性可能增加成本。

## 证据及实际测试

- trial_driver.log / performance.json：一次真实运行COMPLETE、exit0。
- scripts/inspect_engineering_integration.py：exit0；
  source/config/frame/checkpoint hash、实际终时、ID/标签/位姿读回PASS。
- tests/review/test_short_window_e3.py：最终7PASS；
  初版6PASS日志保留，随后增加逐步能量账/原温域限制的只读检查。
  检查覆盖全部帧、真实应力/温度/rho/ep/速度和checkpoint对应；
  不是通过调用同一物理算子生成参考。
- 原476测试1PASS，新局部3PASS/64组通过；没有新求解代码变更。
- 完整116回归、2us、CFL/空间敏感性、GPU/并行、三工况：本轮NOT_RUN。
- 旧失败数据、源码历史、拒绝状态全部保留；文件完整性最终由SHORT_WINDOW_E3_FREEZE.json验证。

## 下一单元：从已接受状态继续原2us目标，不再修476

可在下一轮核对源/配置hash和预算后，使用本次合法checkpoint继续到2e-6s。
剩余约143步，按本次端到端平均预测364.19s，另加20%约437.02s；
拟保持内部540/外部570s，仍一次，当前未启动。
这是有新数据支持的剩余区间估计，不覆盖/删除上轮从初态预测683s的记录。
若后续状态明显变慢或报警即停，不自动延时/改变物理。

建议下一轮命令（NOT_RUN）：

    .envs/pysph/python.exe -B scripts/run_afsd_debug.py --case configs/debug/starved.json --restart outputs/short_window_e3_20261001/dx05_cfl01_05us/checkpoint.npz --output 新目录 --contact-variant event_aware_v1 --end-time 2e-6 --steps 200 --output-every 1 --wall-seconds 540

需外层timeout570。source/config严格匹配才能恢复，不允许悄悄迁移。
本次0.5us处已截短一个步，续跑时间网格须明确记录；
不能把它冒充从零直达2us的完全相同时间网格，也不能据此声称dt收敛。
先完成原短程实现目标，再重新评估走向可观察塑性送料所需的物理时间/成本；
不以本轮0.5us通过授权宏观长跑或三工况。
