# 真实342粒子候选桥接：独立来源与调用审查

历史范围说明（2026-10-04续）：下文为四次评估之前的静态审查，NOT_RUN及50/60s
是当时状态/设计。随后在真实执行前明确修订为110/120s，四次已完成；现状见
PROCESS_BRIDGE_CONTRACT、REPORT、OUTPUT_REVIEW。保留原审查，不倒改试前记录。

范围说明：下文为试前版本，保留当时NOT_RUN及旧50/60s预算证据；不是当前状态。
实际四次结果和事前110/120s预算修订见PROCESS_BRIDGE_CONTRACT、PROCESS_BRIDGE_REPORT
及PROCESS_BRIDGE_RESPONSE_REVIEW；本历史正文不倒改为事后“预知”。

2026-10-03，READ_ONLY_SOURCE_AND_ARCHIVE_REVIEW。
本审查只读实际源码、JSON/gzip与PySPH checkpoint、重新计算选定SHA256；
**没有调用 evaluate，没有运行新轨迹，没有改物理或接受史。**
四次BE候选评估及新驱动回归在本审查范围均为NOT_RUN；其实际结果须另行登记。
约束依据为 `PROCESS_RESPONSE_NEXT_WORK_PACKAGE.md` 与
`PROCESS_GLOBAL_THERMAL_VARIANT.md`，不以既有27点PASS代替342点接通。

## 1. 结论与高风险检查

现有 `ProcessThermalTrial.evaluate` 能提供需要的同史342点塑性、接触、执行器、
BE热候选映射，不需先换后端/本构或移植HHT算法。本包仍只是未提交响应测量。

1. 固定C来自旧**显式热**接受阶段；新BE map不应被强迫等于旧热map，更不能
   强设C残差为零。旧接受本身的maximum_scaled_residual是0.8640199666842818。
2. `accepted_trial.json.gz` 内evaluation的accepted=false/restart_allowed=false
   是正确的提交前语义；公开接受依据是独立acceptance、history和accepted.npz。
   不能用文件名或evaluation.accepted=true作为候选来源认证。
3. `probe_updated_stage.inputs()`认证650ns父态，**不自动认证stage_001候选**；
   也不自动补入新BE模块的源清单。驱动必须补齐下节的闭包。
4. continuation的 `last_accepted_reference.json` 实际仍指650ns起始父态，
   不是最新850ns。最新态应由manifest.last_accepted_snapshot与history最后一项认证。
5. C±若改变塑性/FKS-floor/接触事件分支，只能如实报告响应非光滑；不调整epsilon
   追求相同分支。本切片不能认证整个Jacobian、宏步或输运。

## 2. 输入证据闭包：路径、字段、关系

下列路径均相对项目根。2026-10-03只读实查：continuation manifest登记的43个
source_sha256当前全部一致；下表选定文件也均与 `PROCESS_RESPONSE_FREEZE.json`
的files_sha256一致。这不是全部冻结文件的重审，也不是数值重启等价性测试。

### 2.1 两套配置和三个不同状态不能混合

原父目录 `outputs/coupled_afsd_20261002/short_07us/`：

- `plastic_accepted.npz`：step/count=62，t=6.503936329109041e-7s，342粒子。
  SHA256=`edb75c31c4b38f9d4e22142f29147efd7825fc7af6f51fe84cedc1a4b46da342`。
- 经 `assess_process_route.checked_inputs` 先核COUPLED_AFSD_FREEZE、manifest、
  原source及artifact，再由 `probe_updated_stage.inputs` 核state.accepted、
  count、N、state.source_hashes、state.parameters、state.configuration_hash。
- 原parent所在目录manifest是整段700ns轨迹；被选checkpoint是650ns。
  不能要求二者时间相等。
- parent粒子属性hash（旧evaluation.accepted_state_hash）是
  `ff94e2a7ae4c4a3861e2d6600f322cdaa44c42318fccba8b0a74cd5f997b11ea`。
  此hash不同于压缩NPZ文件hash；两者都应保存和检查。

历史续算目录 `outputs/nonlinear_continuation_20261003/continuation/`：

- manifest.classification=`UPDATED_STAGE_ANDERSON_CONTINUATION_APPROXIMATE`，
  status=`COMPLETE_REGISTERED_FOUR_STAGES_NOT_DEPOSITION`，accepted_steps=4。
- manifest.parent_reference.path/sha256/source_time_s/source_step必须匹配上述父态。
- `resolved_config.json.parameters`与parent参数同值；其configuration_sha256
  是fingerprint({classification,parameters,settings})，不是原parent的配置摘要。
  值为 `f8cb21fea743323fcbd65dbe69c17e958a23eb8a1db9d3e463157751682050ae`。
  不可要求两种不同配置载荷的hash相同。
- 新bridge另设自己的configuration hash，完整包含四调用顺序、dt、方向构造、
  epsilon、输入摘要、预算与分类；不能复用旧continuation hash冒充新配置。

### 2.2 关键文件摘要与接受链

以下短路径相对上述continuation目录：

| 文件 | 实际SHA256 |
|---|---|
| manifest.json | `6a49a0268d47ab437bd45f35d16c6d35c01e260fb74e74d4d06e64c358de390f` |
| accepted_history.json | `c63ef86e6650be4cfd96341ab36a2ccb022dbf65d0dcecc745a751bfdc157c0a` |
| resolved_config.json | `605ce3d50967f3b8ecce453ea5064969c26dff9230c0a071cfd37f34c8f5802d` |
| stage_001/accepted_trial.json.gz | `a327ef8fc4b724ec56243d5b08714423a485af638c3a5993feac3b4f53d0ee22` |
| stage_001/acceptance.json | `fa66e30d4b8ba19ba279dfa1657c3864ee6bdbf09804ecb12d8089baacf4f0a9` |
| stage_001/accepted.npz | `8edf7c62962d04bff26c447e64b7eb7bf4d7cac29de1997d6eefef456ce67c9b` |
| stage_004/accepted.npz | `18d97b64c5bcb960e5a12a7491a849287811551ae89c67ae0fdc0583507966b8` |

驱动需认证下列关系，而非只记录现有文件的新hash：

1. manifest自身匹配事前冻结；每个所用artifact实际摘要匹配manifest.artifact_sha256。
   路径解析必须留在对应目录内，不接受绝对路径或逃逸路径。
2. history[0]的stage_index=1、dt_s=5e-8、time_s=7.003936329109042e-7、
   snapshot与snapshot_sha256匹配stage_001/accepted.npz及manifest。
3. checkpoint state接受/重启语义、参数、source_sha256、configuration_sha256、
   parent_reference及t/stage_index匹配原manifest/history。particle_state_sha256
   为 `3c254fa7a3d5b3d223604e7dacfb0037016f2795f9465ca6ab795424fd344c21`。
4. gzip顶层正是candidate/evaluation；candidate只有登记的五字段，形状分别为
   N×3或N、有限；old evaluation的时间/ID/accepted_state_hash绑定650ns父态。
   old evaluation.classification应是旧UPDATED_CONFIGURATION版本，不是BE版本。
5. **已只读验证**：candidate.position/velocity/density/temperature与stage001
   checkpoint对应字段逐位相同；evaluation.material.S/ep与其checkpoint逐位相同；
   evaluation IDs、parent IDs、stage001 IDs顺序完全一致。mean_velocity没有作为
   独立checkpoint字段保存，须由已认证gzip绑定，不可用(v0+v1)/2代替事件均速。
6. acceptance.status为ACCEPTED_NUMERICAL_STAGE_APPROXIMATE_NOT_ENERGY_CERTIFIED；
   最大归一残差与history相等且<=1。它证明旧离散接受，不证明新BE候选已平衡。
7. manifest.last_accepted_snapshot=`stage_004/accepted.npz`；history[3]为
   t=8.50393632910904e-7、stage4，粒子hash=
   `8053dbda1773fb057115b322c92ed1d789029459207ee5f545c0d64b31feeeab`。
   这条公开最新史只读保存、前后hash一致，不被650ns重放替代。

`continuation_state.restore`可只读核原续算checkpoint，但其docstring也明确：
内部roundtrip还需history公开条目才有公开接受资格。新bridge不得调用commit/save。

## 3. 源码实际调用顺序及最小源清单

核心实读位置：`process_thermal_trial.py::ProcessThermalTrial.evaluate`；继承的
构造/域检查来自`updated_trial.py::UpdatedTrial`；原旧evaluate仅作语句对照。

```text
same parent(x0,v0,rho0,T0,S0,ep0,IDs) + five-field C
  -> validate_interval + candidate shape/domain + endpoint/before normal snapshots
  -> xmid/rhomid/Tmid; clone PA; v_for_gradient = C.mean_velocity
  -> mechanics.corrected_operators (real PySPH LinkedListNNPS dim=3)
  -> material.update(S0,ep0,L,Tmid,dt)                 [先材料]
  -> midpoint deviatoric stress + pressure(rhomid)
  -> paired mechanics.internal_force
  -> coupled_boundary.evaluate_boundary(xmid,v0,...,Fint,initial_z,IDs)
       -> normal_snapshot -> rigid_contact_kinematics.contact_rows -> geometry
       -> axial_actuator.control(initial_z,tmid,...)
       -> per-particle contact_kick.kick(P/m, all simultaneous faces)
            -> contact_cone.solve (function-internal import)
       -> actual event-interval mean velocity + raw impulse/work/reaction ledgers
  -> before/mid/after finite-feature guards (IntegratedSolver borrowed methods)
  -> actual plastic energy=(m/rhomid)*plastic_density
  -> certify_boundary_heat (all actual raw friction intervals, no clipping)
  -> implicit_thermal.advance_implicit_thermal (radial g, T0, raw sources, H, HTwall)
  -> replace thermal boundary ledger with BE Tmap; keep original Tmid diagnostic
  -> density map + x/v/mean-v/T map; map domain checks; residual=C-map
  -> verify private accepted properties hash unchanged; return UNCOMMITTED
```

约束在统一contact kick内参与projected mobility和reaction，不是contact后无账
覆盖速度。每次从同S0/ep0和v0出发；四次间不得把map回写parent或trial._accepted。
这里没有调用任何hht_path/curved/stop/response integrator；新HHT搜索修复不能
被记作该342接口已经使用的功能。

源记录须先保留已认证旧manifest的43项及输入freeze，随后显式增加：

- `src/afsd_backend/pysph/process_thermal_trial.py`：当前SHA
  `a9f385483b5f3225381fa87ba827cd4e83dd88561a38e128275d36ffcb86c232`。
- `src/afsd_backend/pysph/implicit_thermal.py`：当前SHA
  `5f6e77ce2b91587c033e6c9d6c7b1dcdcc1dcb242aa502500d07c7a8cb85233f`。
- 新bridge driver、实际调用的独立array checker、固定方向/预算 helper，
  `PROCESS_RESPONSE_NEXT_WORK_PACKAGE.md`及本次事前验收文档。
- 特别不能漏掉旧表已有的updated_trial、stage_iteration、mechanics、kernel、
  material、coupled_boundary、axial_actuator、contact_kick、**contact_cone**、
  rigid_contact_kinematics、geometry、coupled_solver、implicit_thermal所借用的
  coupled_thermal._vector/_scalar；输出restore/归档helper源也要保留。

source SHA是项目源码证据，不自动证明安装二进制来源。运行receipt还需记录
解释器/NumPy/PySPH/CyArray版本与实际加载路径、PySPH commit、线程环境；
若声称完整运行时闭包，应纳入实际加载的PySPH/CyArray扩展摘要，而非只哈希
third_party中的NNPS文本。`particle_hash`只覆盖PA properties（名、dtype、shape、bytes），
参数/solver_data/常量/模块身份须另外绑定，不能把它称为整个进程的状态hash。

前后都核同一个source/input清单；只在结束重新记录新hash而不比较，无法检测漂移。
原parent PA与trial._accepted各自前后核完整properties；candidate C也保持原始副本。
新调用输出accepted=false、restart_allowed=false、can_publish=false；
诊断用临时文件原子发表，不触碰旧accepted_history/manifest/checkpoint。

## 4. 旧显式热与新BE热：允许不同的字段

相同C、parent、dt下，力学语句顺序/参数保留：同中点NNPS、S/ep更新、压力、
内部力、冻结几何多面kick、执行器、raw摩擦与塑性源。BE变化显式命名，非静默替换。

| 项 | 旧UpdatedTrial | 新ProcessThermalTrial |
|---|---|---|
| FKS温度 | (T0+C.T)/2 | 不变 |
| pair导热温度 | T0 | Tmap（BE解） |
| wall导热温度 | (T0+C.T)/2 | Tmap（BE解） |
| 塑性/摩擦热输入 | chi*P、eta*rawQ | 不变，额外逐事件证书 |
| 热时间步界 | 显式quarter-capacity拒绝界 | former_explicit_dt_bound_s仅诊断 |
| boundary.contact_power_W | signed midpoint conduction | signed BE conduction |
| map.temperature、温度residual | 显式结果 | 新BE结果，不强令相等 |

BE在同中点几何构造对称非负pair conductance及每粒子H、B=sum(Hface*Twall)，
求 `[diag(m cp)+dt(L+diag(H))] ΔT = chi P+eta Q+dt[B-(L+diag(H))T0]`。
能量核对用Tmap=T0+ΔT和实际rawP/Q，不用候选C.T冒充已收敛温度。
旧法计算的是pair(T0)、wall(Tmid)；即使同物理参数也不是同热离散方程。

`_end_temperature_boundary_ledger`仅替换：contact_power_W/contact_conduction_W、
contact_conduction_into_material_J/wall_conduction_out_J及trace热通量；旧量保存在
`midpoint_temperature_conduction_diagnostic_W`。不能把两个层的通量相加重复加热。
normal potential、return_dissipation、kinetic/material残差都不是热源。

新thermal字段是pair_conductance_W_K、pair_heat_into_i_J、particle_energy_J、
ledger、linear_system_residual_J、former_explicit_dt_bound_s等；不能按旧
pair_power_W/thermal_dt_bound_s键硬读并填零。原始热源来源应逐项导出。

`certify_process_thermal(trial,C,e)`认证的是**candidate C.T**的BE余额，之后才
检查map相关账。当前旧C未经新根求解，不能强要求这个候选认证PASS；也不能
替换C.T为map.T后仍称原C通过。bridge应只复算映射内热账并单列完整残差。
`probe_updated_stage.acceptance`是旧显式热接受器，不适用于该新BE候选。

## 5. 成本依据、限制与响应解释

成本证据：`outputs/nonlinear_continuation_20261003/cost_probe/manifest.json`，
SHA=`8d20565e9cde4456a16f4dfb0acb0017926a9da13547403b3691deac7b6f9e78`。
其complete_trial_wall_seconds=4.04492670000036包括一次旧evaluate和gzip保存，
不包括新BE证书、独立读回检查或全部进程初始化。

只读实查原gzip：

| 指标 | 旧cost probe | 固定stage001 |
|---|---:|---:|
| 粒子/邻对 | 342 / 8131 | 342 / 8131 |
| 活跃contact / 接受子区间 | 104 / 6976 | 104 / 7581 |
| 零事件数 | 33 | 28 |
| 压缩候选字节 | 2375638 | 2530558 |

固定C子区间数已比旧probe多8.67%；新BE另外逐段认证并解342阶稠密热系统，
独立checker及四份归档也有成本。因此原 `1.5*4*4.0449267+5≈29.27s` 是
合理的**试前约30s预测**，不是BE实测/保证，也不是下一全局根成本。
维持已登记内部50s/外部60s；首调用的完整evaluate+证据保存+检查结束后，
必须按实际最大完整单次成本重算剩余额度。超过预算则后续NOT_RUN，
保留已经完成的非提交响应；不加方向、不重启第二进程补次数。

三机制是否覆盖必须由本次raw结果证明：真实接触、ep增量/P>0、执行器reaction/功。
旧stage001 PID269塑性功6.484598586830537e-9J和执行器功-5.595200466721313e-10J
只证明选择依据，不能抄作新评估结果。原总自由能仍诊断，不因局部账平衡升级。

四调用顺序C/C+/C-/C、固定epsilon与方向已在NEXT_WORK_PACKAGE事前冻结。
第一/第四相同历史结果可检确定性与不累加；中心差分只证明该完整映射的一个方向响应，
不自证导数准确性。分支更换时保留原始数据并标NONSMOOTH_NOT_CERTIFIED。
最新真实AFSD仍850.393632910904ns、crossings=0；此次650ns历史支路不得增加接受时间。

## 6. 本审查实际完成状态

- PASS（只读来源检查）：43个旧source当前digest一致；表列选定archive及新BE两源
  匹配PROCESS_RESPONSE_FREEZE；stage001四状态字段、S/ep与归档对应一致，ID顺序一致。
- PASS（静态调用审查）：上述调用/热时间层/源字段差别由实际函数确认，非按名字猜测。
- NOT_RUN：四次真实342点BE候选evaluate、其完整账目读回、新bridge接口测试。
- NOT_RUN：新接受步、真实宏步、输运、生产工况；本文件无任何此类验收授权。

本次唯一新增文件为本报告；无需也未修改任何被冻结的物理或历史记录。
