# 持续入口协议组件：NOT_INTEGRATED

2026-10-03。此文件先于本组件测试落盘。范围只有 PySPH 粒子数组的候选插入、
持久 ID 和显式外部源账，不接 NNPS、执行器、求解器、接触或热传导；不是物理供料验证。
不改现有任何物理、冻结文件或旧失败。新增字段必须由调用者完整给定，不能隐式补零。

## 事前验收与独立参考

- 原有粒子所有字段/常量逐位不变；ID、来源、排序对应及 JSON/checkpoint 往返精确相等。
- 独立制造入流：质量 [2,3] kg，速度 [(1,2,0),(-1,0,2)] m/s，
  温度 [310,280] K，cp=10 J/(kg K)，Tref=300 K。
  外源质量 5 kg、动量 (-1,4,6) kg m/s、显热 -400 J、动能 12.5 J；
  密度 [2,3] kg/m3 对应注入体积 2 m3。该常数不是 AA6061 参数。
- 源账与数组前后差值采用绝对容限 1e-12（相应 SI 单位）、相对容限 1e-12。
  这是双精度代数账测试，不是过程能量收支精度承诺。负显热不得裁剪。
- 累积显热必须使用同一 cp 与参考温度；不同约定不得相加冒充同一能量账。
- tracer、wholebar、AS/RS 标签替换不得改变源账；原始标签及坐标不可重置。
- 未解必需物理字段、缺失单位/来源、非有限值、非正质量/密度、schema 差异、
  ID 复用、重复事件、旧 JSON 配新数组、撤回注册表必须显式拒绝。
- 仅一次 <60 s 的测试进程；如 FAIL 最多两轮有证据修复，首 FAIL 保留。

## 接口与物理边界

`initialize_registry(pa, provenance, acknowledgement=...)` 只能用于无先前入口历史的新注册表。
`restore_registry(serialized_json, pa)` 校验注册表完整性及完整粒子数组绑定，排序允许。
`propose_insertion(pa, registry, incoming_pa, specification)` 返回新数组、候选注册表和事件源账。
输入数组不变；传入完整 incoming PA，ID 占位为 0，由注册表分配递增且不回收的 ID。
只有全部 schema、字段、来源和单位通过后才创建候选；不会推进时间或自由粒子轨迹。

字段逐项申明 `unit/source/status`；必需物理场必须是明确初始条件，scratch 可声明
UNRESOLVED 并保留存储值但不能出版为已解物理。schema 外扩不能丢弃旧历史。
相同 material_id 的标签无物性作用；cp 与参考温度均独立声明，不由 tracer 推断。
源账仅包含进入控制体的 m、m*v、m*cp*(T-Tref)、m*|v|²/2，体积为 m/rho。
压力流功、弹性储能、边界反力/外功及混合后的热力学更新为 NOT_COMPUTED；
总能量/焓闭合不是本组件 PASS 范围。允许显热为负，不使用绝对值或零截断。

调用者后续必须确定入口几何、插入触发/重叠检查、有限孔口接触覆盖、一致 NNPS
重建、物理接受条件、原子 checkpoint+注册表发表，以及包含外源的守恒/输运算子。
当前 `transport_observer` 的固定质量条件仍保留；本组件不得绕过它。
外部 checkpoint 只能加载已原子发表的 PA+注册表，孤立候选始终不是接受物理末态。
当前注册表绑定全部粒子历史：正常力学推进也会令原摘要失配。因此完成一次物理步后，
下一次入口需要基于已发表接受历史授权的 rebind/commit 接口；该接口 NOT_IMPLEMENTED。
绝不能通过重建注册表绕开校验，否则会失去此前分配记录。JSON/hash 只验证本地配对
完整性，不提供防恶意回滚的全局日志权威；PA 与旧注册表同时回滚需未来发表日志拒绝。
候选 ID 未经发表不能用作已存在物理粒子的证据；同一旧输入重复提案只是同一分支重算。

## 实测状态

代码实现：完成独立候选协议。最终独立代数/存储测试：14 PASS，非 AFSD。
求解器接入、实际边界供料、输运/沉积与入口稳定性：NOT_INTEGRATED / NOT_RUN。

实际运行（均有 55 s 外部超时，不启动真实算例）：

```powershell
.envs/pysph/python.exe -B -c 'import subprocess,sys,time,json; start=time.perf_counter(); p=subprocess.run([sys.executable,"-B","tests/review/test_inlet_protocol.py"],timeout=55); print(json.dumps({"outer_wall_s":time.perf_counter()-start,"exit_code":p.returncode})); sys.exit(p.returncode)'
```

- `logs/process_completion_20261003/inlet_tests_initial.log`：13 PASS；exit 0；外层 0.410322 s。
- 审查发现不同 cp/Tref 可混入累积显热，添加独立拒绝回归。
  `inlet_tests_thermal_reference_red.log`：14 tests，2 子项 FAIL；exit 1；外层 0.382154 s。
  首错原样保留；它不是旧物理或实测 AFSD 失败。
- 最小修复：注册表锁定首事件 cp/Tref，后事件改变即拒绝，不改热物理。
  `inlet_tests_repair1.log`：14 PASS；exit 0；外层 0.404161 s。
  修复轮数 1，不新增第二后端、不更改现有 solver/observer。

源代码：`src/afsd_backend/pysph/inlet_protocol.py`；测试：`tests/review/test_inlet_protocol.py`。
独立测试同时覆盖旧数组全字段不变、NPZ/JSON 和置换后持久 ID、标签无源账作用、
旧注册表配新数组被拒绝、非有限/未解必需量拒绝；它们不证明三维入口力学稳定。
