# Boundary eligibility v1（2026-10-02，运行前登记）

任务只修已诊断的起始接触分支分类。旧Q7/Q7R/src/配置/测试/失败输出保持不变。
已有1469项冻结校验PASS。新companion复用RecoveryBlock的14预算、实际J2更新、
整体法切切线、原始功和原子commit；不更换本构/后端，不接热/工具旋转/供料。

## 事前算法定义

从已接受态计算C=mu*N和margin=C-|tau|。任何margin<-ft的旧态拒绝；
ft仍9e-5N。只有abs(margin)<=ft且C>ft的接触作为BOUNDARY_ELIGIBILITY，
进入原自然映射联合求解。released表示可由方程选择黏滑，不代表指定滑动速度。
零载荷C<=ft不预释放：保留强制stick的预压/内部事件路径，不宣称一般无载碰撞已解。

用同一整体切线、同一初态求原大小mode_probe_fraction的探测分支；
探测后margin>ft的已候选接触是明确向内卸载，取消资格，保留随后的内部反向事件。
否则仍保留资格；原接受态/反力/耗散不投影、不裁剪。
真正内部的forced-stick行仍用原strict判据、40次二分、1e-12s事件时间定位。
每请求最多12事件、14次Newton、原线搜索/condition/force/slip/work限制全部不变。
原commit仍会拒绝未分段反向滑移、越圆、负原始耗散与不守功；不绕过它。
不能定位的反向或曲线路径允许BLOCKED，不以非负Q作为唯一成功条件。

## 验证先于真实小块

独立参考先登记于PROCESS_PLASTIC_FRICTION_ONSET_REFERENCES.md：
边界立即滑动、向内黏着、真正内部事件、共享耦合邻接激活、
带符号初态误差、二维方向变化、功账、分段和旧RED。标准物理量
参考atol1e-10 SI；小基准的更紧求根只用于解析验证，不修改实际27粒子预算。
容限带内初态不是精确连续体解；记录其早释放/中点求积误差，不授予精确事件时间认证。
另验实际零载荷起始、原commit拒绝不污染、标签无影响和原Q7/Q7R兼容。

## 真实运行预算

先4每相增量（16请求+事件）、27点/9接触/54自由DOF；预计20–50s，
外部subprocess timeout60s，单线程。失败保留accepted/REJECTED及逐请求分支/二分trace。
若完整共同终时.2s、独立读回通过且实际成本允许，才8增量；否则明确NOT_RUN。
一次失败后不继续调物性/容限/预算反复跑。实际AFSD/空间加密/生产仍NOT_RUN。
准静态制造时间不能当显式稳定dt或论文尺度时间；材料能量告警仍保存、不补热。

## 最小源码及接口

新scripts/probe_process_plastic_friction_onset.py：
boundary_eligibility(owner)、advance_onset(owner,t,on_accept=None,trace=None)，
OnsetBlock只替换advance入口；commit/trial/tangent等沿用旧实现。
原main序列化复用，额外onset_provenance.json记录新入口/配置及父恢复来源；
原父16来源不当完整入口清单。resolved_config.onset_variant显式记录改变。
逐请求诊断trace在尝试后保存，不是接受态/重启。所有新失败原样保留。
