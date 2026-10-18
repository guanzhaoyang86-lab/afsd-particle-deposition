# Q7：已有J2塑性＋双切向摩擦的最小耦合

2026-10-01。已实现并独立验证接口、原子提交和耗散分账；实际三维块完成预压段，
但完整加载程序 **BLOCKED**，不能称Q7整段通过或AFSD复现。未改生产src、旧配置、
旧tests/FAIL；PySPH仍1cd3dea9cfaac7d9f1bda2a2491ce4b42a31dbe5。

## 本轮实际实现

新scripts/probe_process_plastic_friction.py/PlasticFrictionBlock继承Q6的state/
branch/advance、圆形双切向自然映射和整体法切耦合；复用Q2的constitutive/trial/tangent
及原material.update，不改本构、不加硬化/阻尼/伪逆。显式选择已有Y90MPa常屈服J2，
不是FKS或AA6061校准。重写的commit先全部检查、后一次写入S/ep/P/R和摩擦状态；
不调用弹性commit后再补塑性、不在Newton/事件探针累计ep或功。

同Q6的1cm立方体/27三维点/9接触/54自由DOF、X0全夹持、固定法向刚壁、
E70GPa/nu0/rho2700、μ.1、293.15K规定温度；壁先压入40um，再沿YZ=(.6,.8)
按0→+20um→-20um→0程序走。没有设置自由粒子轨迹。Y90MPa取自Q2制造参数；
本轮不改变真实几何、物理参数或容限以绕过失败。

实际顺序：Q6.advance_cone纯探针/事件→Q6.branch约束/整体Newton→Q2.trial三维
梯度/实际J2/配对内力→Q6.state接触根/反力→Q7.commit重新trial并独立验接触与
功账→一次提交全部材料/接触→accepted日志/时序/VTK/diagnostic checkpoint。
这是准静态小应变制造问题，没有调用旧动态kick、温度求解或实际工具旋转。

重要新增行为：即使预压时墙面YZ位移为0，塑性横向膨胀仍会引起摩擦。
不能沿用Q6弹性nu0的“预压无摩擦”假设，也不能沿用Q6事件个数。

## 独立验收与证据

CONTRACT/REFERENCES在测试和真实试跑前落盘。核心参考不是用被测函数生成：
实际J2法向材料点的N决定独立二维弹性滑块的当前μN，求解双切向自由位移。
字面终态参考P=.004J、Q=.0025J、Wtotal=.0125J、ep=1/750、N5N；
随后法向卸载且切向壁静止：ΔW=-.0035J、ΔQ=.00042J、ep不变。
这是明确的一向耦合制造基准，不冒充同一材料点多轴塑性的解析解。
低模量toy只用于独立单位测试；真实块E/质量没有缩放。

| 检查 | 实际状态 | 证据 |
| --- | --- | --- |
| 7独立耦合参考，含current-capacity/反向能量/基底/kappa/return分账 | PASS | core_initial.log、q7_tests_final.log |
| 3实际塑性trial/tangent纯度、篡改候选回算/单次提交、塑性后标签不影响力学 | PASS | core_with_actual.log、q7_tests_final.log |
| 11独立reader对抗测试 | PASS，仅I/O/schema保护 | readback_tests_initial.log、q7_tests_final.log |
| 新测试共21项 | PASS，完整进程11.713s | logs/process_plastic_friction_q7_20261001/q7_tests_final.log |
| 原Q6兼容31项 | PASS，不计新测试 | q6_compatibility.log |
| .025s实际预压探针 | PASS，3接受/2事件，5.23s | actual_branch_probe.log |
| 原4增量完整程序 | BLOCKED，6接受至.05s，下一.0625s失败 | n3_i4.log、REJECTED_manifest/candidate |
| 实际已接受输出独立读回 | PASS_PARTIAL_READBACK，exit2，仍BLOCKED | n3_i4_readback.log、outputs下同名JSON |
| 8增量、同终时敏感性、空间n5 | NOT_RUN | n3_i8/NOT_RUN.json |
| 热、FKS摩擦、旋转工具、惯性AFSD、生产重启 | NOT_RUN | 未在本轮求解 |

本节未写目录的日志均在logs/process_plastic_friction_q7_20261001/，运行目录在
outputs/process_plastic_friction_q7_20261001/。11 reader tests用明确标记的
SCHEMA_FIXTURE_NOT_PHYSICAL；它们不是实际耦合力学证据。原Q4严格ep尾差FAIL、
有限应变闭循环/空间FAIL均未修改，本轮未重跑它们。

## 真实块结果与第一次阻塞

4每相增量（16请求＋额外事件）原计划到.2s，预算<50s、外部timeout60s；
实际首次进程6.026s返回exit2，非超时。最后接受step6/.05s，两个事件在
.0217140625944s和.0249001369264s；预压末8个接触滑移、中心1个粘着。
质量.0027kg，M最大条件数3.386；selected arrays235764bytes，不是峰值RSS。
本例无声学稳定dt/真实动态步速意义，不能外推论文尺度AFSD成本。

末态最大ep=.003716473803，P=.239646124388J、Q=.005284602077J；
壁功.306541366349J、材料接触功.301256764272J、应力功.301256764280J。
独立接触总功账最大残差4.35e-17J；内外配对功误差约7.63e-12J。
P_i=Y V_i ep_i与逐步R_i=1.5 G V_i (Δep_i)^2独立读回通过。
R累计.127418079068J是算法return诊断，不能当成额外物理热加入。

Es=.065895989870J，Ws-Es-P=-.004285349978J，超过事前1e-6J诊断线。
step2至6均保留MATERIAL_WORK_DIAGNOSTIC_WARNING；不裁温、不补热、不改判能量
认证。跨材料屈服点的中点功求积和Jaumann/小应变能量兼容性均有可能贡献，
此块体运行尚不能分离二者。单纯接触功账通过不能排除材料能量缺口。

失败候选：base_step6、time.0625s，最大12次Newton后root_relative=
1.15288355628e-6，目标1e-8；未被提交。最初残差.0907447，最后几次
.0066957→.000896383→.0000805904→.00000115288，迭代中求解子块条件数约15.5–27.8。
这次触发器不是negative friction Q、NaN或条件数门槛。原source/config hash、
accepted_0000至0006、完整history、REJECTED候选及Newton trace都已保存。

## 根因证据等级及后续边界

- CONFIRMED：默认12次迭代没有完成此.0625s候选的接触/平衡求根。
- CONFIRMED：预压阶段真实塑性膨胀已产生非零接触滑移，不能把它视为弹性预置状态。
- SUPPORTED：当前阻塞是加载拐点后局部非线性求解成本增加，而非全局切线失秩；
  条件数远低1e12且残差仍下降。是否只需要有界迭代恢复由下面独立诊断判定。
- POSSIBLE：塑性活跃集与径向预压摩擦方向转为壁面斜向驱动，加重非光滑Newton求解；
  本轮未逐点分离该贡献，不按函数名称断言。
- RULED_OUT（仅当前触发器）：并非旧negative midpoint friction报错，也非超时或GPU缺失。
- BLOCKED：共同.2s端点/后续卸载反向全路径、有限应变能量认证、生产AFSD。

只对同一最后accepted状态/同一released mask做固定12与14次纯branch对照，禁止提交
或推进路径。事前约束和结果见PROCESS_PLASTIC_FRICTION_Q7_BRANCH_DIAGNOSTIC.md、
diagnose_process_plastic_friction_q7.py及diagnostic.json/log。它不是生产restart支持，
不能把额外迭代候选移到accepted输出，也不能用它宣称原12次程序完成。

实际诊断已执行：branch_diagnostic_12_14_layout/diagnostic.json，完整进程.929s。
12次逐数组及完整trace精确复现原失败；固定14上限实际第13次即达
root_relative=1.637588335e-10。因而第一层预算不足为CONFIRMED，不只是猜测。
但该收敛分支仍有particle_id21/feature1的强制粘着力超圆约5.69123N，
应交回已有事件定位/分段，不能直接commit！已放开接触行根最大1.474e-6N，
方向检查通过；最小原始candidateQ约-6.4e-20J为舍入尾，不作abs/裁剪。
所有accepted粒子/参考/块体状态指纹未变，时间仍.05/step6，两个新candidate均
REJECTED_DIAGNOSTIC。因此不是“多迭代一下整条路径已通过”。
初次诊断在任何Newton前因恢复数组C/F布局造成合计4.5e-13N舍入差被严格拦下；
branch_diagnostic_12_14.log和branch_diagnostic_restore_inspection.log保留。恢复原力数组布局后才做
这唯一12/14对照，没有放宽record精确复现要求或追加第三个迭代预算。

下一步优先是加载拐点后的有界求解恢复/完整耦合路径，而非再换本构或添加新物理。
保留Q7原12次失败基线及原验收容限；若增加求解预算需显式配置、独立回归、源版本
和失败证据一起保存，不能以无上限迭代解决一切。完整共同终时通过后才做8增量对照。
此前不推进热/供料/转动工具，也不把本例.05s当作真实AFSD物理进展。

## 交付和授权

新增probe/inspect/freeze_process_plastic_friction及独立diagnose脚本、Q7配置、
核心/reader测试、合同/独立参考/本报告、原始日志和accepted/REJECTED分离输出。
原文件仅更新HANDOFF、IMPLEMENTATION_PLAN、OUTPUT_SCHEMA。
PROCESS_PLASTIC_FRICTION_Q7_FREEZE.json与preservation_report只验证证据完整性；
以freeze_final.log实际结果为准，不能代替科学验收。
当前级别：可审计工程近似组件＋真实预压耦合段；不是完整AFSD、不是论文定量验证。
原真实AFSD仍2us/191步、无持续沉积；三工况/长跑/GPU不放行。
