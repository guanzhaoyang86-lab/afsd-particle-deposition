# 联合候选隐式热变体：事前验收（2026-10-03）

分类 `PROCESS_GLOBAL_EVENT_AVERAGE_BE_THERMAL_APPROXIMATE`。
这是明确的时间离散变体，不是论文完整方法，不是已通过宏步/输运验收。
新 `ProcessThermalTrial` 继承 `UpdatedTrial` 构造、五场候选、域限制、同史材料与旧有限特征 guard；
没有修改或 monkeypatch 冻结旧源码。无新本构、接触几何、质量/模量/时间缩放。

## 调用顺序与唯一变化

旧 accepted → 候选中点 x/rho/T 与独立事件平均速度 → 原 corrected NNPS/梯度 →
原 Jaumann/J2/FKS → 原伴随内部力 → 原事件接触/执行器 → 原有限特征 guard →
逐事件摩擦证书 → BE 热图 → 五场 residual。每次重算从同一 accepted 历史出发。

`evaluate` 由于旧接口没有热注入点而局部复制；保留原力学语句的顺序/参数。
逐段差异清单：

1. 删除仅旧 `advance_thermal` 调用，替换为 `certify_boundary_heat` 与既有
   `advance_implicit_thermal`；输入同一旧 T0、候选 rhomid/xmid、原径向 g、实际塑性与摩擦区间能量。
2. 从真实活跃 face 构造 H=sum(hcontact*dx²)、B=sum(hcontact*dx²*Twall)。
   不把旧 Tmid 接触热功率误当 B。FKS 仍使用候选 Tmid，热输出进入候选五场闭合。
3. BE 传导全部使用同一 T1；boundary 的热功率/ledger/trace 热记录更新为 BE，
   原 Tmid 诊断单独保留。接触冲量、速度、几何、机械功不变。
4. 顶层/热层 metadata 显式命名新变体；不将旧 classification 冒充为新接受凭据。
   原显式 dt 界仅报告为 former，不限制无条件稳定的冻结 BE 热子系统。
   这不是机械大步认证，仍保留全部材料、接触、空间、温域和既有全局 residual 门槛。

## 独立证书与容限（先于运行）

- 每接触子区间独立复算 `Q=-Jt·u_t_mid`，重建壁参考系切向中点及旧 256 eps 舍入界；
  不允许跨事件正负抵消后才判断、不允许 abs/max 裁热、不允许外部任意放大 roundoff。
  热值一致容限 `1e-25 J+1e-12*abs(Q)`；证书系数相对 1e-12、绝对 1e-30 J。
  速度切向重建容限 `128 eps * (旧相对速度尺度+两端相对速度范数)`；原有负热判据不变。
- 温度解析绝对 2e-12 K，制造27点重排温度 2e-12 K；材料/机械比较沿旧试验 1e-6 Pa/1e-12状态量。
- BE 候选逐粒子热账 `1e-12 J+1e-10*sum(abs(Q))`，总账沿旧
  `1e-12 J+1e-9*sum(abs(Q))`；热对净和 `1e-15 J+1e-10*sum(abs(pairQ))`。
  温域保持原 room<=T<melt，mapped 与真实候选都检查，不裁温。
- 两节点 C=1、g=1、T0=[300,400]、dt=1，独立 BE 解 [1000/3,1100/3] K。
  旧显式 dt 界 1/4 s；旧实现对 dt=1 必须失败，新热解不能据此声称机械大步可接受。
- 单节点 C=2、H=3、Tw=500、T0=300、总源4 J、dt=.5，独立 T1=1354/3.5 K。
- 真实旧事件接触的一维反向制造参考 Q=7/18 J、v1=-2/3 m/s、事件均速=-1/18 m/s。
- 27点真实 PySPH 更新器：零平衡/非均温 BE 账、候选 T 对 FKS 反馈、标签重排无影响、
  accepted 不污染、重复调用一致、越温失败、旧/新非热物理一致、显式稳定界处区别。

`certify_process_thermal` 仅热/身份/域认证，不授予全局接受。调用者必须另外完整执行
旧空间/力矩/动量/执行器/事件均速/密度/材料历史/五场残差等接受检查。
同模块 `acceptance(pa,p,trial,c,e)` 已执行这些完整原机械检查；它不调用伪造旧 metadata
的 legacy acceptance。其相对旧检查的逐段变化仅为：具名 classification；旧热段替换为
`certify_process_thermal`；移除仅显式热 dt 判据；返回追加 variant/process_thermal_checks。
原机械阈值、Es/体积能/角动量诊断表达式均保留，依旧不是自由能认证。
全部真实 AFSD 新过程步、宏步时间敏感性、过孔输运：本组件范围 `NOT_RUN`。

## 实际运行登记

`thermal_variant_tests_initial.log`：16项中15 PASS、1 ERROR，外0.6337214 s。
错误是基底制造fixture原始150μm穿透，超过未变50μm guard，正确拒绝；不是热模型失败。
只将该制造fixture移至25μm穿透，不修改原真实配置/门槛/物理源码。
`thermal_variant_tests_repair1.log`：16项全部PASS，内0.254 s，外0.6367663 s；使用一轮修复。
独立2点/单点解析、27点非均温图方程、温度对原FKS返回反馈、真实反向事件、标签重排、
不污染历史、旧力学逐数组一致、负热/假舍入证书/越温/机械门槛拒绝均实际执行。
27点1ms静止非均温制造解仅证明联合接口中的热子问题，不代表AFSD宏步、持续变形或沉积。

重跑命令（设置 OPENBLAS/MKL/OMP_NUM_THREADS=1）：

    .envs/pysph/python.exe -B tests/review/test_process_thermal_trial.py

实际执行外包 `subprocess.run(..., timeout=60)`；日志保留命令退出码与墙钟。
未修改任何冻结旧源码；真实AFSD新试跑数量0。不得将这16项再叠加成过程物理步。
