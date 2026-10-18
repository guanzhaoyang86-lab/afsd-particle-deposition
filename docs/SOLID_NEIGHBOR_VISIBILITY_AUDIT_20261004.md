# 有限工具实体遮挡审计（只读）

结论：所查初态、5 ms 与正式 22 ms 端点均未发现材料 NNPS 直线边穿过工具实体。22 ms 的 7,899 条保存候选边与从同一端点重新建立的真实 PySPH 邻域逐边一致。另查旧失败的 26.75 ms 候选，7,859 条边也无实体遮挡；该候选仍未接受、不可重启，未接入正式轨迹。

| 存档 | 边数 | 连续严格实体内部长度 > 1 nm 的边 | 最大实体内部长度 |
|---|---:|---:|---:|
| 独立热态初态 | 8,131 | 0 | 0 m |
| 合格 5 ms | 8,082 | 0 | 0 m |
| 正式 22 ms | 7,899 | 0 | 0 m |
| **未接受** 26.75 ms | 7,859 | 0 | 0 m |

因此上述端点没有遮挡边的参与 PID、组件或最严重实例。0/5 ms 用只读真实 PySPH NNPS 重建；22/26.75 ms 以保存的 `candidate_evidence.npz` 中原始行序 i/j 为依据，并独立重建核对。保存候选位置和粒子 ID 与对应存档数组逐值一致，工具位姿由保存值及参数/时间交叉核对。没有调用 trial.evaluate、材料更新、力计算、导热解或时间推进。

在工具局部坐标中，对每条线段用 z=0/H、x/y=±半孔宽及外圆交点作解析分段。每个开区间的中点判定 `0<z<H`、`x²+y²<R²`、`max(|x|,|y|)>半孔宽`。只接触边界不算穿实体；连续严格内部长度超过 1 nm 才触发拒绝。另对工具实体内缩 0.1 nm 作敏感性诊断，**该内缩不参与放行决策**。当前各存档连任意正长度严格实体区间都没有。计算仍采用浮点交点，不声称机器精度形式证明。

源码证据：`quasistatic_process.py` 在 x1 上调用 operators，直接将同一 i/j/g 传给修正梯度、配对内力和隐式导热；`mechanics.py`、`implicit_thermal.py` 未执行工具可见性过滤。批量变体只替换接触计算。所有处于非零核梯度支撑内的保存边因此进入这些算子；本次零遮挡不能消除后续铺展时出现穿实体边的风险，也不能验证自由表面核一致性。

新增纯函数位于 `scripts/audit_solid_occluded_neighbors.py`，接口 `audit_solid_edge_visibility(x,i,j,pose,p)`。输入必须为完整候选的 actual endpoint 原始粒子行序和实际图；调用方可用返回的 i/j 加入 PID。异常应停止发布；`occluded_edge_count>0` 应停止发布，不删除边或改变内力。接口只诊断端点图，不保证整个时间段的运动图可见性。

`scripts/check_solid_edge_visibility.py` 执行 11 项纯几何 fixture：孔内通路、穿墙、肩面相切、孔口绕棱、孔口切入、圆柱相切、沿孔壁、两侧墙穿越、刚体位姿变换及 0.5/1.5 nm 阈值。全部 PASS；不计作 AFSD 物理验证。22 ms 的 7 次 guard 测量最大耗时 0.0170891 s。

可重跑命令（输出已存在时按设计拒绝覆盖）：

```powershell
.\.envs\pysph\python.exe scripts/check_solid_edge_visibility.py
.\.envs\pysph\python.exe scripts/audit_solid_neighbor_snapshots.py
.\.envs\pysph\python.exe scripts/audit_solid_neighbor_failed_candidate.py
```

三条最终执行均 exit 0。主快照审计内部耗时 0.0849329 s；旧失败候选补充审计 0.0418916 s。主审计首次运行因报告试图从参数字典取不存在的 h 而退出，未写审计输出；改为读取实际粒子 h 后完成。曾尝试通过 restore 读取拒绝存档，按设计被拒；补充脚本明确使用只读 raw load，不绕过任何接受或重启门槛。

证据落盘：`outputs/coupled_cost_20261004/solid_edge_visibility_receipt.json`、`solid_occluded_neighbors_audit.json`、`solid_occluded_neighbors_unaccepted26_75ms.json`。三个 JSON 记录来源 SHA256、接口限制与时间成本。纯守卫冻结 SHA256 为 `e496a8756348844696643e465155a903c9ef68ae1d44c0cd5b8867207da69f21`。未修改既有数学源、已冻结 driver、接受存档或旧失败记录。
