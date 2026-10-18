# 下一整体接口：更新构形候选与过程步（设计，未实现）

沿用同一PySPH pin、材料/NNPS/力/热/ID/输出。不是换后端或给自由材料规定路径。
本轮已实现的是process_residual.midpoint_residual，不是下文的全局求解器。
本设计不修改冻结配置/物理，不给长时、生产或准静态切换授权。

## 优先路线

以保留真实质量/惯性的隐式候选为第一路线；只有另有慢载、平衡和误差证据才研究去惯性近似。
不能把700ns末态速度清零，不能因名义工具Mach数很小就删除预压释放中的惯性。
隐式稳定性也不等于大时间步精度，特别是旋转、塑性、接触切换和热软化。
原Q7的27点0.2s制造准静态程序不能拿来作新方法成本证明。

## 一个可复用、无提交的入口

建议新增命名 `trial_increment(accepted_state, candidate, dt, parameters, rigid_program)`。
candidate包含x1/v1/rho1/T1、实际接触未知量与执行器反力；不是预设运动答案。

1. 克隆接受态到独立PySPH ParticleArray，不修改accepted数组、ID、累计P/Q或checkpoint。
2. 明确采用中点构形xmid、vmid、rhomid；重新调用mechanics.corrected_operators，不能沿用Q7冻结D/B/V。
3. material.update每个Newton候选都从同一S0/ep0调用；候选温度层明确，真实dt进入FKS率项。
4. 同一构形/rhomid/B/g下求pressure和伴随internal_force；记录功共轭身份式及cond。
5. rigid_contact_kinematics.contact_rows取同一工具位姿与候选几何；工具方孔开放，板棒均变形，原始z0仅选轴向执行器。
6. 组装本轮Rx/Rm/Rc，再加入连续方程 `Rrho=(rho1-rho0)/dt+rhomid*tr(L)`；不能继承Q7不更新密度。
7. 接触力来自同一候选的法向/摩擦方程，不冻结上一显式kick的冲量当全局解。
8. 最后才求/核热候选及所有分账，返回NOT_ACCEPTED；任意失败保留诊断且接受态不变。

该入口作为下一工作单元先实现/测试；本轮没有伪造其试跑结果。

## 不能遗漏的耦合工作

- 接触：全局约束后的effective response，不是每粒子P/m限幅照搬到静平衡；
  真实滑移反向、stick可行性、多面/法向开离/feature变化须分事件。
  本轮反例已证明总冲量满足动量也不能用端点速度造热；只对实际接受子区间计 `-Jt·u_mid`。
- 时间运动学：whole-step drift与event-resolved速度轨迹可能不同，需统一子区间的位移、应变和热功，
  不能仅修Rm而保留互相不一致的运动学。隐式方法不消除这一要求。
- 热：复用coupled_thermal的原始径向导热和源分配；若宏步超过真实热行和界，
  明确子循环或独立验证的隐式热解；原函数只拒绝过大dt，不已有自动子循环。
- 本构：Jaumann/Cauchy自由能缺陷保持风险；中点时间算法只能针对时间离散，不能修掉非零闭循环功。
  有限应变free-energy候选须另命名、独立测试，不在本接口顺手替换。
- 全局切线：Q7局部材料有限差分只能贡献本构部分；更新NNPS/B/几何/法向/约束后还需对应导数。
  不称其为完整一致切线，先做方向导数检查；病态/零模态失败，不pinv或偷偷增加材料夹具。
- 入口与孔口边缘：当前有限棒/预结合板棒不是持续沉积程序，真实供料的ID/质量/焓账、
  完整有限工具棱边和热启动依据仍单独欠缺，不能靠数值积分方法自动补齐。

## 独立基准与执行顺序

下一轮仅完成纯候选接口及以下少量联合基准，再决定一次全局试步；避免拆成许多同义单测刷进度。

| 基准 | 参考及要捕获的问题 |
|---|---|
| 更新三维仿射场 | 独立矩阵L与体积变化；确认NNPS/M/B来自候选构形，而非只改输出位置 |
| 非零初应力刚转 | 用独立sin/cos轴角给Q和Q S0 Q^T；显式区分连续旋转与有限步中点运动学误差，不调用被测rotation造答案 |
| 离散功/密度 | 同一候选的Σv·f+ΣVσ:L，Rrho与小步解析体积；不假定有限步中点连续方程对任意J恰好精确 |
| 材料trial隔离 | 多次相同/扰动/拒绝候选不能重复累计ep、塑性功、热或改变accepted hash |
| 接触事件功 | 原反向滑移手算、双面+轴向约束；验证真实分段运动学/热而非只验总冲量 |
| 热宏步 | 两点容量20/60J/K，C=1152/35W/K，初温400/300K；T1=325+75exp(-384t/175)，T2=325-25exp(-384t/175)，同终时折半子步比较 |

之后的全局试步必须先测一次trial/Jacobian成本，明确自由未知量、迭代/事件预算与预计墙钟；
超过600s不启动。首次global trial不必产生可视化，不从未平衡静态替换制造“沉积”。

## 精确复用位置

- src/afsd_backend/pysph/process_residual.py：本轮新增Rx/Rm/Rc/Rq纯组装。
- mechanics.py：corrected_operators/internal_force/pressure；condition_limit=1e6继续。
- material.py：现有FKS、Jaumann/J2返回；b1/率floor近似和独立return账不变。
- rigid_contact_kinematics.py、geometry.py、axial_actuator.py：真实边界/同点功/原始材料轴向控制。
- coupled_thermal.py：物理P/Q/接触/导热源；不得接受机械残差当热。
- scripts/probe_process_contact.py::reduced_correction：可借鉴完整受约束切线求解。
- scripts/probe_process_plastic.py::trial/tangent、probe_process_feasible_predictor.py：可借鉴纯trial/有界迭代思想，
  不能继承其制造Y、冻结构形或固定接触方向。

方法资料是MOOSE的[应力残差/构形一致性](https://mooseframework.inl.gov/moose/source/kernels/StressDivergenceTensors.html)、
[迭代材料返回](https://mooseframework.inl.gov/source/materials/RadialReturnStressUpdate.html)、
[独立惯性项](https://mooseframework.inl.gov/source/kernels/InertialForce.html)；
只作为设计依据，不安装该后端或声称它验证了本实现。
