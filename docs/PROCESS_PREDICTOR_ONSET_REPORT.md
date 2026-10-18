# PredictorOnset：首次贯通当前J2塑性—双切向摩擦完整小块程序（2026-10-02）

本轮结果：**J2_CONE_PATH_COMPLETE + PASS_COMPLETE_READBACK**。
在原固定几何/物性/边界下，从初态完成预压、切向正向、反向、返回四阶段，
共同终时.2s，16个请求区间经24次内部起滑事件分为40个接受步，41帧。
唯一完整试跑76.7377572s、exit0，未到外180s上限，也没有重试或数值预算增加。
这不是AFSD沉积时间；层次是ENGINEERING_APPROXIMATE组件耦合基准，
不是完整AFSD、有限应变能量认证或论文定量复现。

## 最小改动与冻结的边界

新scripts/probe_process_predictor_onset.py只有branch适配和原输出入口包装：
PredictorOnsetBlock继承OnsetBlock；branch调用上轮已验证的
solve_with_feasible_predictor。advance、boundary eligibility、严格内部二分、
commit、trial、tangent、state、wall全部复用原实现。
预测+cone共享14修正，line search8、events12/请求、bisection40、原接触容限不变。
config声明Q7_PREDICTOR_ONSET_V1，不修改旧源码、配置、失败和生产默认。
新增reader验证23项来源（parent16+Onset4+本轮wrapper/config/旧helper3）
及精确配置差异，然后实际调用原Q7物理输出审查，不调用求解器生成参考。

不变：1cm立方体、27粒子、9固定法向接触、54永久自由DOF；
E70GPa/nu0/rho2700kg/m³、制造J2常屈服90MPa、mu=.1、
法向压入40um、壁切向幅值20um、方向(.6,.8)、每相.05s。
293.15K仅规定温度，未求热；固定法向范围不支持正常开离/再接触。
自由材料位移由整体平衡求解，非规定自由粒子轨迹。
原压力/材料/paired corrected空间算子未替换，condM=3.385745。
没有裁负耗散/法向/温度、没有残差热、没有质量/时间/模量缩放。

## 实际测试与先验参考

先验容限/运行预算见CONTRACT，独立参考见REFERENCES；未试后改门槛。
线性独立弹簧mu10，只是制造基准，不是实际mu或论文参数：
起滑壁位移31/15m、路径时间62/75s，终态s31/15,q8/15,N7/15,tau14/3；
Q91/45J、Wwall1139/180J、Wmaterial155/36J。
实际被测组合使用原advance_onset，事件前黏着、事件后滑动，共两次原子提交；
超摩擦圆黏着根、整段强行释放跳过事件、负法向和损坏候选仍拒绝。

日志目录logs/process_predictor_onset_20261002/：
| 项目 | 结果 | 日志/墙钟 |
| --- | --- | --- |
| 独立组合解析/功/原子提交 | 8 PASS+1旧guard RED | independent_initial.log，.418s |
| 实际类/预压/标签/提交防护 | 5 PASS | actual_integration_tests.log，6.854s |
| 来源/预算/读回合成元数据 | 25 PASS；不是25物理算例 | readback_tests_final.log，2.333s |
| 新测试合并 | 38 PASS+1 RED（上述不重复计数） | integration_tests_final.log，7.928s |
| 原predictor兼容 | 13 PASS+1旧RED | predictor_compatibility.log，.372s |
| 原Onset兼容 | 27 PASS+1旧RED | onset_compatibility.log，7.353s |
| 唯一原4增量完整尝试 | COMPLETE，exit0 | n3_i4.log，76.738s |
| 实际输出读回 | PASS_COMPLETE_READBACK，exit0 | n3_i4_readback.log，.498s |
| 无求解器独立数组复核 | PASS；30字段旧物理前缀不变 | independent_audit.log |
| 8增量/空间/热/旋转工具/真实AFSD/生产 | NOT_RUN | n3_i8/NOT_RUN.json |

expectedFailure是被保留的原失败，不计PASS。其余历史套件本轮NOT_RUN；
全部旧FAIL仍在原路径和完整性清单中。测试通过不授予AFSD生产资格。

## 不再卡在旧.1125失败点

原.125目标的八个强制黏着接触超限，现在由已有事件处理器逐个定位：
新事件时间约.112888862/.113218568/.113974155/.115256401/
.115538107/.116358419/.117116186/.117429928s。
并非一口气释放八个接触或直接提交上轮的超限候选。
随后完成.125及后续反向/返回程序，最终达到.2s。
独立NumPy/JSON审计确认对应新起滑ID依次19/24/26/22/20/21/25/23，
feature均1。旧history全部30字段（8静态+22时序前17帧）逐字节相同；
旧accepted0–16物理数组不变，仅迭代诊断schema增加阶段/法向字段。
全程只有.125的TARGET分支用过一次predictor，1预测+4cone；不是隐式批量预移动。
实际attempt时间链：16请求、886返回的诊断候选根、40真正接受步；
候选数不能算成物理步。各候选与事件仍有独立trace和状态分类。

## 最终结果和仍在的科学限制

末态累计：
- 外壁功.348086595285J，接触材料功.304901214090J；
- 摩擦耗散.043185381195J，塑性功.255421450146J；
- stress work .304901214069J，Es诊断.053827737946J；
- 最大ep=.004181979456，体积平均ep=.002838016113；
- mass=.0027kg，固定法向gap和约束按原规则审查。

全历史最大接触law残差6.828055e-5N（原9e-5N内），
最大机械自由残差5.467811e-9（原1e-8内）。
独立累计接触功恒等式偏差最大1.110223e-16J；
raw摩擦增量最小-9.968003e-20J原样保留，不做abs或裁零。
这些说明接触/外功记录一致，不能证明材料自由能正确。

**Ws-Es-P=-.004347974022J告警仍存在**；
39个接受步带材料诊断warning，reader明确material_energy_certified=false。
旧Jaumann闭循环与空间失败仍未解决。return_dissipation_diagnostic不是热源。
FKS/温度依赖、导热/塑性功转热、工具旋转、正常开离与供料沉积均未被本次验证。
只有一档时间划分、单一空间分辨率，不称收敛或论文对照通过。

## 性能与下一项

运行前预计90–150s，实测76.738s（advance76.042s），新外部硬上限180s；
该上限事前声明，旧父60s保持历史来源，不影响数值或物理。
选定数组235764bytes，不是进程峰值RSS。
这里是准静态事件程序，不能拿.2/76.7直接预测显式AFSD论文尺度成本或CFL步速。

下一工作单元：保持物理几何/材料/边界完全相同，跑8增量到同一.2s，
检查反力、位移、塑性功、摩擦功与事件位置的时间划分敏感性；
先冻结测量定义和工程容限/告警分类，不能以两个不同终时声称收敛。
按本轮886 roots/24events成本，下一同规模粗估90–150s，仍需新轮确认预算；
不要无依据按粒子数外推论文大规模成本。
完成这个有界对照后，再讨论旋转壁/局部送料集成；不继续围绕已越过的初猜点调参。
正常开离、有限应变、热/FKS与旋转方孔适配仍是后续明确缺项。
本轮8增量按事前范围NOT_RUN，未启动三工况、GPU或生产。

复跑命令见本轮COMMANDS.md。新源码/日志/输出与交接落盘，不依赖聊天历史。
额外独立审计脚本audit_process_predictor_onset.py不导入求解器，计算signed功、
ID、前缀和预算；明细见PROCESS_PREDICTOR_ONSET_INDEPENDENT_AUDIT.md。
