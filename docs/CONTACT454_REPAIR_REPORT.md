# 粒子454接触修复与后续439阻塞（2026-09-30）

## 结论与范围

**454局部问题已修复并由独立参考验证；完整耦合仍为APPROXIMATE/BLOCKED。**
本轮只修改一个物理源码文件contact_kick.py的数值表示/事件精度处理。
保持PySPH commit `1cd3dea9cfaac7d9f1bda2a2491ce4b42a31dbe5`，上游未改。
未改变材料、corrected梯度及伴随内力、热学、几何、参数、标签、NNPS或积分器。
没有GPU/并行优化、质量/时间/模量缩放、温度/屈服/负热裁剪。
以下路径均相对项目根。前次32/485原负热诊断仍见COUPLING_FAILURE_REVIEW.md；
该次历史报告的“3/454失败”由本报告更新，但原失败证据没有覆盖。

最终实跑107项：**101 PASS、6 FAIL、0 ERROR，退出码1**。
5项原方法反例继续FAIL，新增439独立回归FAIL；没有删除、xfail或改绿旧反例。
唯一一次新的固定真实试跑接受26步，第27候选在439双面接触锥迭代失败，
31.57秒停止，未达2µs。禁止生产、CFL/空间收敛结论和论文复现标签。

## 1. 预注册与调用位置

预注册：docs/CONTACT454_REPAIR_CONTRACT.md；
旧245项冻结在修改前验证PASS：logs/contact454_20260930/freeze_before.log。
修前源码备份：logs/contact454_20260930/source_before/。
旧contact_kick.py SHA256：
`7d4afd640bce07bcfe4d4caf6478cc1c83f9737d5216db8ba9f07d97f79b58fd`。
本轮修后/真实运行SHA256：
`14d12701ddfcfa60118ebeb2edf4fdbd34b1c07a17851dfeb581ec9c236f4b2d`。

实际调用未重排：

    PySPH half drift
      -> solver.evaluate: NNPS / corrected L / material candidate / internal force
      -> domain_kick: analytic contact faces + actuator projection
      -> kick: instantaneous stick/slip -> event interval -> midpoint cone_solve
              -> interval actual impulse/work -> local accepted endpoint
      -> material/thermal commit -> PySPH velocity kick + half drift
      -> global alarms -> accepted checkpoint/output

454和439的异常都在局部cone_solve内，未作为有效末态提交。
这不是先接触再任意覆盖执行器速度；局部求解沿用P/m约束后迁移率，反力功保留。
默认legacy不变，只有显式event_aware_v1使用本修复。

## 2. 454证据与最小修复

修前只读重放再次失败，日志probe_before.log，报告outputs/contact454_20260930/before/report.json。
输入为旧REJECTED接触局部上下文，不把旧REJECTED当作仿真重启态。

| 判定 | 证据 |
|---|---|
| CONFIRMED：旧3D投影有法向消去误差 | Decimal 70位独立平面投影；约1e-13m/s切向速度混入-7.51144e-18m/s法向分量 |
| CONFIRMED：旧迭代陷入两点循环 | 最后4次更新相互交替，幅值1.95263749e-11N，不随迭代降低 |
| SUPPORTED：重复投影/消去导致该循环 | 相同冻结输入改用2D正交切向坐标后通过，力误差6.05e-12N；没有放宽容限 |
| RULED_OUT：454材料/几何必须修改才能通过 | 同一内力、法向力、墙速度、质量和dt的局部重放已经PASS |
| BLOCKED：据454通过推断所有多面接触稳定 | 新439反例证明不能作此推断 |

具体补丁：

| 文件/函数 | 接口与原因 |
|---|---|
| src/afsd_backend/pysph/contact_kick.py::tangent_bases | 新增每面3x2正交基T，消除I-nnᵀ对近零切向量的反复消去 |
| 同文件::cone_solve | base从N×3改为已投影N×2；保留全部T_aᵀ(P/m)T_b耦合，返回每面3D力 |
| 同文件::kick | 先投影滑移/加速度再相加；记录近零事件但不覆盖速度；热仍为实际-Jt·u_mid |
| 同文件::kick的underflow门槛 | 仅把64eps×全kick的启发式改为t+h/2、t+h实际能否推进；新的制造反例触发旧门槛时h仍有124个时间ULP |
| tests/review/review_artifacts.py及3个旧测试 | 只改变新证据输出目录；旧输入/失败记录不覆盖；454测试额外检查独立停止时间 |
| scripts/run_coupling_review.py | 纳入新基准并按运行目录隔离产物；任何FAIL仍使总入口exit1 |

原力相对容限2e-12、锥上限2000、区间上限4096、rtol2e-7/atol1e-12、
128eps近零速度尺度均未放宽。真实库仑滑移反转继续允许。
最初新增制造测试的失败保留在precision_patch1.log、precision_underflow.log及对应输出；
记录原因后才修underflow，precision_patch2.log五项PASS。不是反复改物理参数拟合答案。

454独立常力向量库仑极坐标积分给出：

| 量 | 独立参考 | 实际 |
|---|---|---|
| 停止时间/s | 5.985479017657219e-9 | 5.98548031831218e-9 |
| 摩擦热/J | 8.164870224344708e-12 | 8.164870224516668e-12 |
| 最终切向速率/m/s | 0（连续解） | 2.0897160631081497e-15 |

停止时间相对误差2.17e-7，热相对误差2.11e-11；267个接受子区间。
速度未强制归零；事件记录残速2.08846e-15m/s、原零阈值2.41448e-15m/s。
证据：logs/contact454_20260930/final_suite/artifacts/candidate454.json。
坐标旋转和更紧rtol1e-10的结果也落盘，不依靠单个方位通过。

## 3. 实际测试矩阵

全部日志：logs/contact454_20260930/final_suite/；summary.json含命令/退出码/计数/源码hash。

| 组 | 实跑 | 结果与限制 |
|---|---:|---|
| contact454_precision | 5 | PASS：Decimal锥、约束有效质量、双面共享自由度、旋转/紧容限454、制造曲线停止 |
| event_contact | 12 | PASS：解析反转、静/动墙、共转、反力功、多面执行器、旧485、拒绝重启等 |
| new_candidate_RED（保留历史名称） | 1 | PASS：原454输入及独立热/停止时间；不是expectedFailure |
| closed_elastic_cycle | 1 | PASS仅指重复确认非零残功0.146111111111J；材料热力学仍FAIL |
| spatial_probes | 4 | PASS：仿射、畸变/自由面、功共轭、合力/矩相关基准、cond门槛及单云零能模探测 |
| event_solver_smoke | 3 | PASS：短链路/标签、初始成本等 |
| original_scientific | 38 | 35 PASS/3 FAIL：raw刚转、raw异密度功、本构Es反例未改 |
| legacy_reversal_RED | 1 | FAIL保留 |
| boundary_patch_RED | 1 | FAIL保留，常应力节点边界patch残差约32.6% |
| contact439_open_risk | 1 | 新FAIL：可行闭式KKT解存在，当前迭代仍超预算 |
| original_backend / original_scaffold | 21 / 19 | 全PASS |
| 真实固定2µs | 1次 | FAIL于候选27，仅26步有效 |
| CFL=.05/.025、dx=.25mm、生产、宏观对照 | 0 | NOT_RUN |

1D反转独立参考仍为t_stop=1/3s、u_final=-2/3m/s、Q=7/18J、
W_other=1/9J、ΔK=-5/18J。没有abs(Q)、max(Q,0)或残差补热。
闭循环16/32/64/128子步的残功趋向非零常数；本轮不换本构，
docs/FINITE_STRAIN_VARIANT_DESIGN.md仍为DESIGN_ONLY，不能因测试脚本PASS误判模型已一致。

## 4. 唯一真实试跑与新的439阻塞

事前门禁/预算：logs/contact454_20260930/preflight_gate.md。
N=657，初始dt=1.04906e-8s，约191步，首步0.906659s，线性预测173.17s；
预算600s。此估计不保证后续多面求解成本。本轮真实运行一次，失败即停。

实际输出outputs/contact454_20260930/fixed_dx05_cfl01/：

- 接受26步，t=2.727463701412469e-7s（0.272746µs）。
- 候选27拒绝：ID439，BORE_X_MINUS和BORE_Y_MINUS，两面共享粒子自由度。
- 实际墙钟31.5709408s，接受步中位耗时0.87482305s，峰值工作集229756928B（约219.1MiB）。
- dt范围[1.04900139192e-8,1.04906949246e-8]s；未缩放时间/质量/模量。
- checkpoint.npz为第26有效步；27帧含初态。REJECTED_candidate.npz与REJECTED_contact_debug.json独立存放，禁止作为物理末态/普通重启。
- 独立输出读回PASS（output_integrity.log），这只验证文件一致性，不把status=FAILED改为PASS。

末个接受步的诊断：质量2.217375e-4kg，密度最大相对变化1.12816e-4，
最大穿透1.97718e-6m，温度[723.14999603,723.15004652]K；累计摩擦热3.67084e-7J。
原始最小局部热-3.83599e-26J（舍入尺度）保留为负，没有裁零。
当前DEBUG总账残差4.01918e-7J、相对0.0030564；这不解决已知有限应变自由能问题。
ep_max=0，说明这段极短运行没有证明塑性沉积/混合已发生。

439的冻结局部锥可用**独立闭式凸KKT解**验证，不需真实长跑：

    lambda = coefficient/m
    coefficient = 2.6637601040312884e-15 s
    m = 3.375e-7 kg
    limits = [0.3895337242053182, 0.3879904620898042] N
    tangent base = [[-2.639703260799982e-9, 1.1107056659778313e-9],
                    [3.122135038184536e-9, 1.1107056659778313e-9]] m/s
    independent f0 = [-base00/lambda, -base01/lambda]
                  = [0.33445198355952604, -0.14072707286974026] N
    independent f1 = [-limit1, 0] = [-0.3879904620898042, 0] N

第一圆盘内部余量0.0266808086N，第二圆盘在边界；KKT梯度
[约0,0,5.98727755e-11,0]，第二约束乘子为正。因此不是无解。
CONFIRMED：当前投影梯度在2000次后更新仍1.58684206e-8N，
远高于7.79067448e-13N门槛；与454的两点舍入循环不同，它缓慢下降。
SUPPORTED：内点/滑移面转换与共享切向自由度造成慢收敛。
更详细条件数/收敛率归因尚未独立定量确认，不能称为已完成根因修复。
tests/review/test_contact439_open_risk.py保留FAIL；reference不调用cone_solve生成。
这项只验证被冻结锥子问题，尚不声称验证完整439局部kick轨迹。

下一最小任务：只在cone_solve内研究有残差验收的局部active-set/KKT求解，
先过439闭式参考，再过454/反转/多面约束/原始功测试；不得只增加迭代上限或放宽门槛。
本轮未实现该下一修复，也未继续启动其他真实试跑。

## 5. 可重跑命令与证据保护

在项目根执行；输出目录必须用新的名称。sandbox账户读取上游git时可使用进程级设置，
不写用户/系统git配置：

```powershell
$env:GIT_CONFIG_COUNT='1'
$env:GIT_CONFIG_KEY_0='safe.directory'
$env:GIT_CONFIG_VALUE_0='E:/Lagrangian Trajectory to Eulerian Field/Repo/AFSD_Codex_Starter/afsd_codex_starter/third_party/PySPH'
& .envs/pysph/python.exe -B scripts/probe_contact454.py --output outputs/contact454_recheck_01
& .envs/pysph/python.exe -B scripts/run_coupling_review.py --output logs/contact454_recheck_01
```

探针当前应exit0；总入口当前应exit1（107项含6FAIL），不要把exit1改为成功。
本轮实际真实命令如下，仅记录，不建议在439门禁未修复时重复启动：

```powershell
& .envs/pysph/python.exe -B scripts/run_afsd_debug.py --case configs/debug/starved.json --contact-variant event_aware_v1 --output outputs/contact454_20260930/fixed_dx05_cfl01 --end-time 0.000002 --steps 10000 --output-every 1 --wall-seconds 600
& .envs/pysph/python.exe -B scripts/validate_physics_outputs.py outputs/contact454_20260930/fixed_dx05_cfl01
```

原命令exit1，读回exit0。resolved_config SHA256：
`a50aa7974a18287ef0717639f7ff0f40f2d1ddb59389fe4bbea9c3ee1862a255`。
每次试验源码hash独立记录；真实运行后没有再改物理源码。
新证据冻结docs/CONTACT454_REPAIR_FREEZE.json；旧两个freeze文件及旧失败产物保留。
preservation_report.json列出允许的改动和未改核心/旧证据校验。
冻结PASS仅表示文件完整，绝非科学验收或生产许可。
