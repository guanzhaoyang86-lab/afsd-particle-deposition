# 真实已存响应的搜索矩阵作用：试前范围与验收

2026-10-03，在本次代数检查运行前登记。只检查同一650.393633ns父态、dt=50ns、
已保存的真实plus/minus五场候选与真实map；绝不调用evaluate/contact kick/材料/热
更新，不再扫描epsilon，不构造零候选或伪物理mapping。500ns复算本包NOT_RUN。

输入为 `outputs/process_bridge_20261003/real342_initial/{report.json,
configuration.json,plus_candidate.json.gz,minus_candidate.json.gz,base_candidate.json.gz}`
及其已绑定 `outputs/coupled_afsd_20261002/short_07us/plastic_accepted.npz`。
源/归档摘要先核原report，结束再次比较同一清单。父态properties、参数、ID完全不变。

最多构造两个现成近似搜索矩阵：StagePreconditioner和ContinuityStagePreconditioner；
构造使用当前物理质量/模量、父态几何、原50ns、实际中点执行器，允许现有静态NNPS/
冻结弹性矩阵装配。每个矩阵仅对**已存真实(C+,F+)和(C-,F-)**各apply一次。
记录它们返回的correction，不利用proposal发动轨迹。

定义dC=C+-C-，dR=(C+-F+)-(C--F-)，dcorr=correction+-correction-，
E=dC+dcorr。冻结P线性时dcorr=-P dR；E描述这一实际方向上的搜索失配。
保存五场dC、dR、dcorr、E，单位分别m、m/s、kg/m³、K、m/s。
不除以epsilon以免掩盖原浮点分子；旧epsilon=1e-7m/s来源照录，不改变。

逐五场报告绝对inf/L2及原root缩放inf/L2：x/1e-12m，v及vmean/1e-9m/s，
rho/(parent.rho*1e-10)，T/5e-12K。分组：all、free_dof、constrained_dof、
plastic、no_contact、single_contact、multiple_contact、any_contact、zero_event、
no_zero_event。向量的free/constrained按分量掩码；标量按该粒子是否有相应自由度，
明确重叠，不能当质量分割。其它分组只从base的原trace/ep分类，不进入算子。

消减比为||E||/||dC||（同字段同组同范数）；dC=0时不造比值，标
ZERO_INPUT_CROSS_RESPONSE，并报告生成的E。组重叠/单位不同，不合并成物理误差百分比。
**不要求P误差小于任何新精度阈值**；P不是完整J逆，单方向且原分支不光滑，
不能据消减比授予导数、根、500ns精度或物理接受。只检查：

- 全部输入/输出有限、严格同源/同ID、父态和参数不变；proposal=c+correction逐位一致。
- 现有apply内部线性后向误差界1e-12不变；输出原测量值。
- 规定轴回代等于各候选原map-candidate的规定轴缺陷：1e-12m/s浮点一致界，
  不是缩小/扩大原物理接受门。raw five-field根门完全不动。
- 保存所有消减、放大或零输入交叉响应，不能以其大小筛选/重做方向。

独立event B/C全库存不存在，StageEvent对照=NOT_AVAILABLE；不调用旧有限差分kick
补齐，不移植27点HHT矩阵。内部55s、外部60s；单次调用，仅此一次。
成本预估：已测原冻结结构build约0.777s，两个装配加四次代数apply及归档预计<15s；
仅估计，非本次实测。超时停止，部分产物标未完成，不追加进程补做。

运行前状态：NOT_RUN。完成报告另见PROCESS_BRIDGE_ACTION_REVIEW。
