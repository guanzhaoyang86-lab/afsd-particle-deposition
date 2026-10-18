# 独立数组复核：4/8 增量（2026-10-02）

结论：`PASS_INDEPENDENT_ARRAY_CROSSCHECK_NOT_CONTINUUM_VALIDATION`。
两组均到同一制造小块终时0.2 s；4组40步、8组56步，17个共同规定时刻
逐一唯一匹配，未插值、未按相同步号或不同失败终时比较。
这不是AFSD沉积结果、连续体解认证或收敛阶证明。

## 范围与独立性

审查者先只读核对AGENTS/HANDOFF、实际冻结源调用和4组数组，审阅运行前合同
`docs/PROCESS_PREDICTOR_SENSITIVITY_CONTRACT.md`，后复核新比较器。
曾指出读回证据只比步数/终时可能错配；主代理随后补了final/events/warnings/
sources绑定，以及共同帧REQUESTED检查。没有改冻结求解器或物性。

本次重算脚本为`scripts/audit_process_predictor_sensitivity.py`，只用stdlib和NumPy，
不导入、不调用被测比较器、接触或材料更新器。先独立计算，再读取comparison.json
核对结果。输入：

- `outputs/process_predictor_onset_20261002/n3_i4`
- `outputs/process_predictor_sensitivity_20261002/n3_i8`
- `outputs/process_predictor_sensitivity_20261002/comparison.json`

检查每个report列出的artifact SHA256；两组配置去除increments_per_phase后相等，
ID/初始标签/参考坐标/参考体积/质量/支承分量/接触ID/feature相等。
来源清单和完整读回的更广泛核验仍以冻结reader/主比较器及其日志为据；
本独立脚本不声称再次完成材料/空间连续体验证。

## 实测差值

以下百分数是运行前冻结的双路径峰值归一化差，不是真值误差：
`D=max(peak4,peak8,0.01*nominal)`。末态仍用同一D。
所有10项都在预声明5%筛查内，不据此授予生产或定量复现认证。

| 场/量 | 共同时间最大绝对差 | 最大归一化差 | 末态归一化差 |
|---|---:|---:|---:|
| 粒子位移向量 | 1.2255944e-7 m | 0.270508% | 0.231477% |
| 接触合力向量 | 5.960414 N | 0.063315% | 0.026287% |
| 接触切向合力 | 3.066442 N | 0.357715% | 0.029409% |
| 应力张量Frobenius范数差 | 1.137703e6 Pa | 0.806827% | 0.587661% |
| 粒子等效塑性应变 | 1.1890010e-5 | 0.284315% | 0.266786% |
| 体积平均等效塑性应变 | 8.5547397e-7 | 0.030139% | 0.013828% |
| 累计塑性功P | 7.6992658e-5 J | 0.030139% | 0.013828% |
| 有符号摩擦功Q | 5.6140402e-5 J | 0.129941% | 0.044044% |
| 累计壁面功 | 0.002436586 J | 0.695137% | 0.693819% |
| 累计应力功 | 0.002418187 J | 0.786874% | 0.785165% |

接触合力由所有粒子的contact_force求和，不使用支承相抵后近零的reaction_N。
末态接触合力，4组为[-8214.442838,493.277185,656.845526] N，
8组为[-8211.981040,493.083353,656.684325] N。
末态P分别0.2554214501459417/0.25545677349555135 J，
Q分别0.043185381195100875/0.0432044099808658 J。
最大ep分别0.004181979456166672/0.0041738553463655445。
另以独立式`P_i=90e6*V_i*ep_i`核对完整保存路径。

## 事件与有符号功

两组各24个CONE_LIMIT事件组、30条接触边界事件；四阶段组数均2/4/9/9，
条目数8/4/9/9。按phase/持久ID/feature/occurrence核对，inventory、组内分组和
先后顺序相同。最大时间移位0.0001803998289688115 s，发生于第一阶段四个
角接触ID19/21/25/27，8组比4组早；低于预声明0.00125 s筛查。
这些只是新到锥边界候选，不能解释成全部黏滑转换或已发生滑移。

独立逐区间循环用assembled contact_force计算
`Q_i=-F_t_mid dot(du_i-dwall)_t`，没有使用保存的Q来生成参考。
重构Q与保存值最大差4组4.34e-19、8组2.17e-19 J；
`Wwall-Wmaterial-sum(Q)`最大区间残差分别1.17e-17/1.34e-17 J。
未用abs(Q)或max(Q,0)。原始负分量分别25/23个，最小值分别
-9.9680031e-20/-6.2643135e-20 J，符号保留且未越过原-1e-10 J停止阈值。
支承功两组均0。应力功与接触/支承功累计账差最大分别
2.03718e-11/4.91648e-11 J。

材料告警没有解除：末态`Ws-Es-P`从-0.004347974022177514 J变为
-0.002018196803874439 J，告警步分别2–40（39步）/3–56（54步）。
这里只有两级，不能断言趋零或已修正自由能；两组material_energy_certified均false，
没有把差额加成热。旧Jaumann闭循环及空间FAIL不因本轮比较而改变。

## 重跑和实际状态

在隔离环境、OPENBLAS/MKL/OMP线程数均1，用外部subprocess timeout60 s执行：

```text
.envs/pysph/python.exe -B scripts/audit_process_predictor_sensitivity.py --coarse outputs/process_predictor_onset_20261002/n3_i4 --fine outputs/process_predictor_sensitivity_20261002/n3_i8 --comparison outputs/process_predictor_sensitivity_20261002/comparison.json
```

真实日志：`logs/process_predictor_sensitivity_20261002/independent_audit.log`。
exit0，外部实测0.1783685 s；这是保存数据审计时间，不是物理试跑时间。
独立脚本对两组实际数据运行PASS；本脚本的专门单元测试NOT_RUN，主比较器的
独立fixture测试见本轮测试日志，不混为一套证据。

本审查允许结束这个固定法向小块的4/8筛查，转向预声明的旋转壁/局部送料组件
接口和独立基准设计；不建议为改善比例继续调这个已完成算例。
空间敏感性、法向开离/再接触、旋转工具、热耦合和AFSD生产本轮NOT_RUN。
