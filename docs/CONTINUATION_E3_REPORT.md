# E3-4：2微秒实现窗口完成，转向工艺时间尺度决策

日期2026-10-01。范围与预算先冻结于CONTINUATION_E3_CONTRACT.md。
PySPH固定1cd3dea9cfaac7d9f1bda2a2491ce4b42a31dbe5。
本轮没有修改src、物理配置、旧测试、上游或任何历史失败证据。
仅新增只读检查/观测脚本、独立回归、一次合法续跑及文档。

## 实际完成，不等于沉积复现

| 项目 | 实测结果/状态 |
|---|---|
| 原物理问题 | starved，dx=0.5mm、CFL=0.1、657粒子，event_aware_v1 |
| 输入 | short_window_e3_20261001/dx05_cfl01_05us，accepted48步/0.5us |
| 新运行 | continuation_e3_20261001/to_2us，143新增步，终点191步/2us，COMPLETE，exit0 |
| 单次预算 | 运行前预测364s，加20%为437s；实际413.0687763s，内540/外570s |
| 链总量 | 191步/192唯一接受帧；同一48步接缝验证后去重；累计wall535.3133484s |
| 内存/步速 | 本段峰值111.69MiB；物理步中位2.8213s；含输出约0.3462步/s |
| 全局拒绝/报警 | 无REJECTED candidate、无超时，143个接受诊断均无报警 |
| 首次自由材料塑性 | step63，t=0.657351465us；终态3粒子，max ep=1.691590689e-5 |
| 真正供料/沉积 | 穿越质量0；孔外层内自由棒料候选体积0；尚未证明持续沉积 |
| 复现层级 | ENGINEERING_APPROXIMATE短时三维热力实现探针；不是定量验证或70%完成 |

上表wall是runner报告；外层570s进程保护实际未触发。没有第二次真实试跑。
0.5us父段结尾存在截短dt；真实链时间网格保留，不冒充未中断直达2us或dt收敛。
新段dt范围1.043115910e-8至1.049032542e-8s；成本用两段诊断dt中位数，
不用父段终点的极小截短dt估稳定步长。性能文件原有paper_cost_projection仅历史乐观模型；
本轮决定采用mechanism_report_v2的当前657粒子实测外推，不据它声称论文全尺寸成本。

## 实际调用及数据核验

源码：src/afsd_backend/pysph/application.py、solver.py及其调用函数。

```text
PySPH host half drift
  -> material predictor + midpoint stress/internal force
  -> coupled event-aware contact/actuator kick (反力一并记账)
  -> commit material state
  -> conduction + plastic/friction/contact thermal sources
  -> host kick / remaining half drift
  -> diagnostics -> accepted VTK/checkpoint
```

本轮不改变该顺序。contacts内部可有试探子区间拒绝；这不等于全局REJECTED末态。
末步有213个局部拒绝子区间，不能把“全局无拒绝”解释成所有局部迭代一次通过。
原始舍入级负摩擦热仍在诊断中：本段最小-6.104353610e-24J，未取abs/裁零。

scripts/inspect_afsd_chain.py检查父checkpoint运行前receipt、每帧hash、
manifest/config/source一致性、每字段按持久ID的接缝相等和最终checkpoint。
只有accepted到达目标且runner/manifest均COMPLETE才给ACCEPTED_TARGET_REACHED。
192帧全部核验；首次/终末导出位置、速度、应力、温度、密度、ep、ID等是真实数据。
局限：子段初始内部scratch checkpoint未单独归档；不能由导出字段相等证明
未保存内部量或任意分段策略的restart-equivalence。现有NPZ只加载本项目可信产物。

## 物理观察与风险

三个塑性粒子均为自由棒料、原始中心x=y=0，原z为0.25/0.75/1.25mm：

| ID | 终态ep | 累计塑性功[J] |
|---|---:|---:|
| 541 | 7.738145165e-6 | 1.467572887e-9 |
| 542 | 2.935776512e-6 | 1.380939693e-9 |
| 543 | 1.691590689e-5 | 3.208388571e-9 |

它们本来就在1.5mm肩面以下，不能当作新送料或出口流动证据。初始2um过盈、
预结合SPH界面和启动瞬态可能触发塑性（SUPPORTED），未经因果对照不能唯一归因。
无基板塑性，原始负ep/负plastic_work计数0。规定进给总量仅2.12nm，工具仅转0.0036度。
自由粒子最大原点位移158.739nm，分段位移92.198nm；不能混淆。

末态温域723.149967760至723.150166786K，最大自由速度0.621644606m/s。
本段最大密度相对变化1.68119578e-4，最大接触穿透2.02662626um；
末态质量0.0002217375kg。未缩放质量/模量/物理时间，未改几何/层厚。

热账独立核验使用配置cp=900J/(kg K)、T0=723.15K、chi=1：
sum(m cp (T-T0))=5.879461459245761e-7J；
chi*sum(plastic_work)+friction_heat+contact_heat=5.879461414681838e-7J。
差4.456392272e-15J小于事前保守舍入预算9.843999517e-11J。
接触导热保持负号-8.420269121e-7J，摩擦入材料热1.423916152e-6J。
数值返回耗散1.748872227e-9J单列，未变成热；总能量残差也未变成热。
这是此snapshot的全局账目核对，不是热扩散精度或普适本构认证。

重要未解风险：总账面能量残差从0.5us的0.4499%升到2us的3.3109%，
末值4.353837491e-6J，虽未触发原5%+1e-7J保护，仍不能称高精度。
Jaumann变体积闭循环非零功、自由面/支承误差和E2应力FAIL全部保留。
有限棒不是保守持续入口；外圆/顶面接触不完整；预结合界面不证明真实焊合。

## 实际测试矩阵

本轮新增28项全部通过，另重跑父证据7项。证据目录logs/continuation_e3_20261001。

| 入口 | 结果 | 日志 |
|---|---|---|
| verify_scientific_freeze --manifest SHORT_WINDOW_E3_FREEZE.json（编辑前） | PASS，907项 | freeze_before.log |
| test_short_window_e3.py | 7 PASS | parent_checks.log |
| test_afsd_chain.py | 15 PASS，含失败状态不能升级/篡改/ID重排 | chain_tests_completion_gate.log |
| test_mechanism_observation.py | 7 PASS，独立体积/坐标/成本/标签/重排 | mechanism_tests_final.log |
| test_continuation_e3_physical_evidence.py | 6 PASS，独立热账/真实塑性位置 | physical_evidence_tests.log |
| inspect_afsd_chain.py / measure_afsd_mechanisms.py | PASS，读192帧；无推进 | chain_readback.log / mechanism_readback_v2.log |
| 观测器最终版本重算与v2报告全字段相等 | PASS | mechanism_final_reproducibility.log |
| 真实AFSD续跑 | COMPLETE，exit0 | trial_driver.log |
| 完整116项物理回归 | NOT_RUN；引用上一冻结111PASS/5旧FAIL，非本轮执行 | CONTACT476_SCALE_REPORT.md |
| CFL=.05/.025，dx=.25mm同终时 | NOT_RUN：当前预算不足；不声称收敛 | PROCESS_TIME_DECISION.md |
| GPU/并行、生产三工况、新积分变体 | NOT_RUN | 未启动 |

初版观测日志/JSON保留；v2仅澄清ROI身份和字段名，数值未改，见OUTPUT_SCHEMA。
完整命令见logs/continuation_e3_20261001/COMMANDS.md。
最终新旧文件保全记录CONTINUATION_E3_FREEZE.json与preservation_report.json；
完整性核验不等于科学验收。下一单元见PROCESS_TIME_DECISION.md，不再逐粒子修补或延长微秒探针。
