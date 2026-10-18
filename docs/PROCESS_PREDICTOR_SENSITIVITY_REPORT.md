# 同一物理问题4/8增量对照（2026-10-02）

本轮完成。冻结PredictorOnset两档均到共同终时.2s，读回PASS_COMPLETE_READBACK。
新8增量唯一尝试88.626309s、exit0，未重试、未调整物理或求解器参数。
17个共同请求时刻的10项预注册指标均WITHIN_SCREEN，最大归一化差.807%。
这是ENGINEERING_APPROXIMATE组件敏感性检查，**不是完整AFSD或收敛证明**。

## 范围与证据

基线outputs/process_predictor_onset_20261002/n3_i4保持40步/41帧。
新outputs/process_predictor_sensitivity_20261002/n3_i8有56接受步/57帧。
两档均27点、1cm立方体、9固定法向接触、54永久自由DOF，
PySPH 1cd3dea9cfaac7d9f1bda2a2491ce4b42a31dbe5不变。
E70GPa、nu0、rho2700kg/m3、制造J2屈服90MPa、mu.1、规定293.15K，
法向40um、切向20um、YZ方向(.6,.8)、4相各.05s全部不变。
仅increments_per_phase4变8；请求间隔.0125变.00625s，不是CFL/dx变化。
32请求经24个CONE_LIMIT事件组形成56接受步；889返回候选根不是物理步。
14共享修正、事件/法向/功检查不变；不裁raw负Q、温度、屈服，不加残差热或缩放。

事前估90–150s，外硬上限180s（用户600s内），实测外88.626s/advance87.847s；
旧4外76.738s。选定数组两档235764bytes，不是进程峰值RSS。
不据准静态制造程序的.2s外推显式AFSD成本；这里不是沉积时间。
新目录有history、57个accepted journal、checkpoint、VTK、23项来源、32次attempt。
checkpoint仍restart_allowed=false，非生产重启文件；温度未求解。

测量/单位/范围事前见PROCESS_PREDICTOR_SENSITIVITY_CONTRACT.md；
启动前hash见logs/process_predictor_sensitivity_20261002/prerun_hashes.log。
合同SHA256 b2a3e89aee5cfd21063c2169e180d9cc1547b0b36a2daf58c4f3a1d497578dbf。
比较脚本只用NumPy，不导入求解器，不调用contact生成参考。
两档23项来源一致，resolved config除增量数外相同，静态几何/ID/标签/质量/
约束/接触面相同。只匹配唯一实际accepted REQUESTED帧，不插值、不比不同失败终时。
旧源码/配置/结果/FAIL保留；旧i8/NOT_RUN是上轮历史，不覆盖。

## 实际测试

日志统一在logs/process_predictor_sensitivity_20261002/。

| 实际执行 | 结果 | 日志/外墙钟 |
| --- | --- | --- |
| 新后处理独立合成回归 | 27 PASS，不是27物理算例 | comparator_tests_final.log / .211s |
| 原组合兼容回归 | 38 PASS + 1历史RED | integration_compatibility.log / 7.966s |
| 唯一8增量完整尝试 | J2_CONE_PATH_COMPLETE，exit0 | n3_i8.log / 88.626s |
| 原reader审查实际新数据 | PASS_COMPLETE_READBACK，exit0 | n3_i8_readback.log / .547s |
| 同终时数组比较 | COMPARISON_COMPLETE，exit0 | comparison.log / .186s |
| 第二份独立纯数组审计 | PASS_INDEPENDENT_ARRAY_CROSSCHECK_NOT_CONTINUUM_VALIDATION，exit0 | independent_audit.log / .178s |
| 其他历史全套、第三时间档、空间、热、旋转、AFSD生产 | NOT_RUN | outputs同后缀/NOT_RUN.json |

新独立reference：向量(3,4,0)范数5、非对角张量Frobenius范数7；
制造账Q=3.25J、Wwall=4J、Wmaterial=.75J，均非调用被测contact所得。
覆盖缺帧/歧义/非单调/禁止插值、ID/phase事件及不同分组/缺失、5%筛查/近零分母、
负roundoff原样保留/明显负值拒绝、物性与静态篡改、partial拒绝。
tests_initial.log与最终日志重复同27项，不相加；expectedFailure不计PASS。
其余历史失败只保留，未重跑的不写PASS。

## 同终时测量

“归一化差”=共同17时刻峰值差 / max(两档峰值,.01*事前名义尺度)，不是实验误差。
向量Euclidean、应力Frobenius、标量差绝对值；物理Q/P/功保留原符号。
末态差也用同一全路径分母。力用contact_force_N合力，**不用**支座抵消后的reaction_N。

| 量 | 全共同路径最大绝对差 | 归一化差 | 末态归一化差 |
| --- | --- | --- | --- |
| 粒子位移 | 1.2255944e-7 m | .27051% | .23148% |
| 壁接触合力 | 5.960414 N | .063315% | .026287% |
| 壁切向合力 | 3.066442 N | .35771% | .029409% |
| 应力张量 | 1.1377035e6 Pa | .80683% | .58766% |
| 逐粒子ep | 1.189001e-5 | .28432% | .26679% |
| 体积平均ep | 8.554740e-7 | .030139% | .013828% |
| 累计塑性功P | 7.699266e-5 J | .030139% | .013828% |
| 累计摩擦功Q | 5.614040e-5 J | .12994% | .044044% |
| 累计壁功 | .002436586 J | .69514% | .69382% |
| 累计应力功 | .002418187 J | .78687% | .78516% |

逐共同时间/分母/原末态见outputs同后缀/comparison.json。
5%事前标准只作工程提示，不是新求解门槛；未事后放宽或拟合。
两档24事件组对应30条phase/particle_id/feature/occurrence边界到达条目，
ID顺序相同、无缺失；最大时间偏移.000180399829s（提示尺度.00125s内）。
它们是CONE_LIMIT新到摩擦边界候选，不是全部黏滑切换或任意曲线/旋转工具认证。

8档末态P=.255456773496J、Q=.043204409981J，壁功.350518560106J，
应力功.307314150077J、Es诊断.053875573385J，maxep=.004173855346，
体积均值ep=.002838408594，质量.0027kg。
独立重建raw Q最小-6.264313e-20J，23个负roundoff条目原样保留；
区间接触功恒等式最大偏差7.589415e-18J，reader累计偏差2.324529e-16J。

## 限制与下一项

Ws-Es-P由4档-.004347974022J变为8档-.002018196804J；两档均超过原1e-6J告警。
8档step3–56共54步warning，material_energy_certified=false。
两点减小不能证明趋零，不推翻旧Jaumann闭循环/空间FAIL。
return_dissipation_diagnostic .065441739971J不是热源。
三维自由度小块不等于三维沉积，真实AFSD历史仍2us/191步，无持续沉积。
没有温度/有限应变/正常开离/FKS/旋转方孔生产认证。

本项有界对照完成，下一推进**刚壁旋转运动学与接触功接口**：
先检查已有geometry.tool_pose/local_xy/surfaces的接触点、法向、壁速和开放方孔，
独立验证共转零滑、静壁/移动旋转壁功共轭，再设计coupled cone的每接触点接口。
不把当前统一wall(t)平移向量冒充旋转；不强制拖自由材料或让标签参与力学。
局部送料、正常开离、热和有限应变分开列为后续缺项，不本轮一次实现。
不再为此小块调参/补第三档追求全绿；不授权三工况生产、GPU或长跑。

完整命令见logs/process_predictor_sensitivity_20261002/COMMANDS.md。
第二份独立计算见PROCESS_PREDICTOR_SENSITIVITY_INDEPENDENT_AUDIT.md；
audit_process_predictor_sensitivity.py没有调用上述比较器，直接从数组重算后核对。
其专门单元测试NOT_RUN；真实数据审计与新27个比较器合成回归不能混为同一证据。
