# 500ns全局求根末候选：只读数值诊断

2026-10-03。基于已终止 `outputs/process_global_20261003/root_attempt`；
没有运行新full-trial、材料更新器或contact kick，没有改dt、阈值、模型、预算或旧哈希文件。
脚本 `scripts/assess_process_global_nonlinearity.py` 首次成功，外0.6171852s，exit0。
原始证据：`outputs/process_global_20261003/nonlinear_diagnosis.json`；
执行日志：`logs/process_global_20261003/nonlinear_diagnosis_initial.log`。

## 先区分预算未收敛与真正停滞

原运行32次达到事前预算，状态 `BLOCKED_NO_NEW_ACCEPTED_STATE`，新增接受步0。
最后五次scaled raw residual为528.42、236.45、153.98、107.02、72.553；
它仍下降，证据不支持将此称为发散或稳定停滞。
但最后candidate不是root，更不是可续算接受态，不能因误差“小”放宽标准。

| 原始字段 | 第32次残差 | 原标准 | 末候选字段判断 |
|---|---:|---:|---|
| x1 / m | 1.8138269e-14 | 1e-12 | 字段内 |
| v1 / m/s | 7.2553095e-8 | 1e-9 | 未满足 |
| rho相对 | 1.4989623e-14 | 1e-10 | 字段内 |
| T / K | 1.0231815e-12 | 5e-12 | 字段内 |
| vbar / m/s | 3.6276535e-8 | 1e-9 | 未满足 |

这些“字段内”不是物理接受PASS；能量/根/提交需全部条件。
线性预条件自由块backward error=1.283699e-16，不是线性代数精度耗尽。

## 源码＋末态数组支持的根因分级

**CONFIRMED：当前结构预条件没有实际事件接触速度响应导数。**
它用惯性自由弹性中点消元A，保留执行器规定轴，未加入事件平均速度对internal-force、
墙运动或接触几何变化的导数。旧接触kick仍物理驱动，没有被预条件器替换。

**SUPPORTED：漏掉多面近即时停滑的接触响应是剩余机械缺陷的主要位置。**
以末候选 `sum_i ||r_vbar,i||²` 为诊断权重，不是物理能量：

- 99.5357%集中在本kick出现zero-slip事件的粒子；72.6664%集中在多面接触粒子。
- 无接触粒子只占0.0001722%。因此先围绕全域材料率地板再调参不合证据。
- 最大四点PID317/341/197/221分别对应面(3,4)/(3,5)/(2,4)/(2,5)，均无执行器。
  每点两法向力约1.29N，两个停滑事件在3.60–3.62e-13s内；相对本次500ns间隔约7.2e-7。
- 对这些点，实际vbar与端点均值在x/y可相差约1.4e-7m/s，
  大于剩余vbar分量残差；不能改成端点均值来“消掉”它。

这些是定位证据，不是已计算完整Jacobian或证明唯一原因。
保存的zero-event标签之后是stick-feasibility检查，不能仅凭事件标签宣称整段严格粘着。

**RULED_OUT，限末候选：当前塑性分支与FKS率地板切换。**
最终342点plastic_work均0。独立从保存L、oldS做Rodrigues弹性试算和字面FKS公式重算，
所有点都在弹性侧，最近屈服差仍2.915884MPa；率重算与保存值逐位相合。
所有点距rate_floor=1/s至少1.294748/s；四个最大残差点rate约2.59–2.73/s。
这不能排除更早候选出现过其他分支，但末候选附近不应优先换本构或平滑FKS。

**RULED_OUT，限本次已保存边界跨度：有限feature切换。**
before→mid、mid→after无contact_id增删；此结论不能推广为孔口/肩缘已完整支持。

**POSSIBLE但未隔离：法向几何罚刚度、客观率/预应力和更新几何导数。**
每活跃面dt² kn/(4m)=0.121759，现预条件未含该项；不能只凭其小于1就证明可忽略。
空间cond(Mi)=2.21507，功伴随差1.11e-16W；这排除了此次空间矩阵退化，
不证明更新几何Jacobian可忽略。rho已满足原容限说明没有可见rho方程主导残差，
不是证明压力对机械残差完全无作用。

**NOT_RUN：完整切线分解/同分支方向导数与额外真实求根。**
本诊断只读；不以相关性代替因果试验，不在既有32次预算后续跑。

## 下一最小有效联合步骤

下一工作单元应把原kick的实际事件平均响应
`B=d(vbar_actual)/d(F_internal)` 接入候选Schur，并分别保存
`C=d(v1_actual)/d(F_internal)`；不是将物理kick的P/m替换成弹性W，
也不是把真实粘着判据写成永久删除材料自由度。
事件响应只用于非线性候选预条件/切线；物理map、原接受条件、独立vbar保留。

先核独立1D反向区间：m=1,u0=1,F=-2,muFn=1,dt=1，
令tau=1/(1-F)，u1=(F+1)(1-tau)，
vbar=tau/2+(F+1)(1-tau)²/2。
tau'=1/(1-F)²，故F=-2时

    tau=1/3, tau'=1/9
    d(u1)/dF=(1-tau)-(F+1)tau'=7/9
    d(vbar)/dF=tau'/2+(1-tau)²/2-(F+1)(1-tau)tau'=19/54

因此B既不等于自由中点1/2，也不等于C/2=7/18。
这些参考数应字面写入新独立回归，再测试多面＋执行器的受约束响应及共享交叉耦合。
局部数值差分若跨停滑/活跃面分支，必须明确reject或给广义导数范围，不能混成平滑切线。
仅当该联合基准通过且成本可负担，再另注册一次原500ns问题的求根；
不能自动扩充本次32次预算，也不降低raw门槛。

即使下一500ns根通过，也只是全局过程求解器接线的中间证据。
项目主要缺口仍是工艺时间尺度推进、有限开放孔口/肩缘及真实供料输运，
不是再累积几十个微小接受步或给现有0穿越结果补动画。

## 可重跑只读命令

```powershell
.envs/pysph/python.exe -B scripts/assess_process_global_nonlinearity.py --case outputs/process_global_20261003/root_attempt --output outputs/process_global_20261003/nonlinear_diagnosis_new.json
```

首次由外层 `subprocess.run(..., timeout=60)` 包装并保存stdout/stderr；
fresh输出目录文件规则防止覆盖初次诊断。脚本验证输入与源hash不漂移，不做真实advance。
