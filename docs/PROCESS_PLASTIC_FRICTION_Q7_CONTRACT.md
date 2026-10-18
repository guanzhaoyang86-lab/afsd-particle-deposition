# Q7运行前合同：已有J2＋双切向摩擦，单一固定法向小块

本轮只新增companion，旧src/配置/tests/FAIL不改；PySPH固定版本继续复用。
先验父归档PROCESS_FRICTION_CONE_Q6_FREEZE.json1381项，再独立参考，再实际小块。
新类PlasticFrictionBlock继承Q6ConeBlock的实际三维trial/tangent、自然映射/整体约束解、
事件/输出设施。仅显式选Q2已有PerfectPlasticBenchmark(Y90MPa,T293.15K)，不是FKS。
保持Q6几何/质量/E70GPa/nu0、固定法向/μ.1/墙面程序；基底0度。
基准1cm块不是AFSD几何粗化，不以其.2s制造程序时间冒充AFSD物理时间。

## 提交、能量与容限

全部Newton/模式/事件探针从同一accepted S/ep试算，不能累加塑性；commit重新试算后
统一保存ep、逐粒子plastic_work、return_dissipation及逐接触friction_work。
保留Q6力/间隙/圆锥/自然根/逐接触rawQ验收；负法向、NaN、病态切线继续BLOCK。
所有旧数字原样：力相对1e-8、scale9000N、ft9e-5N、st1e-11m、功1e-6J、
逐接触Q<-1e-10J拒绝且原始小负尾差保留；Newton12/线搜8/事件12/二分40/1e-12s，
condM1e6、condK1e12，不pinv、不添刚度/硬化/阻尼、不裁ep或耗散。
方向误差检查沿用Q6量纲传播。material.update返回的负plastic/return量或ep下降拒绝。

分开报告 Wwall、Wmaterial、Wsupport、Ws、Es、P、Qfriction、return diagnostic。
Wwall-Wmaterial=Q，Wmaterial+Wsupport=Ws仍须在1e-6J以内。
Ws-Es-P在1e-6J以上明确MATERIAL_WORK_DIAGNOSTIC_WARNING，不转热、不改写字段；
沿用Q2的材料/时间求积诊断地位，不冒充热力学自由能验证。跨屈服点的中点应力功
有求积偏差；Jaumann有限应变旧风险未消除。return诊断不是新的物理耗散项再加一次。
独立基准必须在解析拐点分段，先证能量/原子提交；块体不保证每个材料点屈服界已分段。

独立reference见同前缀REFERENCES。字面单位基准参考默认位移/力/功atol1e-10SI、
ep1e-12；actualblock pure-trial/双重提交/标签不变要精确相等，已有Q4严格epFAIL不改。
对真实块只给可支持的离散算法/功账结论，不要求它等同于一维串联系统。

## 预算与依赖关卡

Q6完整27点路径约16s；J2迭代及预压横向滑动可能更贵或触发旧完美塑性奇异风险。
先独立测试、再4每相增量（16请求＋事件，估<50s；外subprocess timeout60）。
若完成且实测成本允许才8增量；否则8标NOT_RUN或保留真正失败，不更换物性躲失败。
这轮不为PASS反复改参数，最多针对明确新接口缺陷做小修；不重写本构/接触算法。
同终时才比较增量；不能承担或不适用的空间n5、热、FKS、转工具、惯性、AFSD均NOT_RUN。
任一Rejected候选单独存储并标REJECTED；已有接受态不能因失败探针被污染。
checkpoint诊断用途restart_allowed=false。最终旧文件只改HANDOFF/IMPLEMENTATION_PLAN/
OUTPUT_SCHEMA；不授予生产或论文复现认证。
