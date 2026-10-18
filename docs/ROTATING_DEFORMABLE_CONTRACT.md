# 旋转肩面—三维可变形小块：事前合同

日期2026-10-02；ENGINEERING_APPROXIMATE，非AFSD沉积/论文验证。旧物理、FAIL及
ROTATING_CONTACT_FREEZE保持；只新增companion，不修改生产入口。

## 五步及边界

1. 冻结本合同/SI配置/真实调用映射。2. 复用原NNPS/材料/伴随内力/接触/积分接入。
3. 独立解析首步、三维算子、标签/ID/重启/拒绝回滚回归。4. 测成本后唯一粗档短程。
5. 独立读回粗档，通过后唯一同终时半步长对照；报告差异，不声明收敛阶。

固定64个中心点(4^3)，dx=.5mm，实体单元域2mm立方体；中心x=.01m、y=0，
z中心.00025/.00075/.00125/.00175m。每点m=2700*(.0005)^3=3.375e-7kg，
总质量2.16e-5kg。G=26.3GPa,nu=.33沿用旧配置，无质量/时间/模量缩放。
应力、速度初始零；三维邻域必须包含非共面邻居，内层32点xyz完全自由。
底层16点固定xyz(制造夹具，不代表整块刚体)；顶层16点仅z固定、xy受力自由。
肩面z=.00195m，中心gap=.0002m、检测overlap=.00005m；Fn=1N/顶点规定，
mu=.2、omega_z=20rad/s、走速0。Fn与法向夹具反力分存，非penalty/正常开离解。
方孔宽.008m仍开放但本块在孔外；不以肩面特例认证旋转孔壁功共轭。
初始接触面积对应平均压力4MPa；几何固定，不跑空间加密。

293.15K为规定等温输入；未解温度不输出为温度场。原FKS/Jaumann更新不换；
b1量纲与rate_floor等旧未决保留，短程预计弹性但实际ep仍记录，不伪造为零。
Q为原始摩擦耗散账项，未沉积成热；Es仅旧应力平方/体积诊断，非已认证自由能。

## 实际调用链与独立核对

HostLeapfrog半漂移→PySPH三维NNPS/corrected_operators→旧material.update(L_before)
→rho_mid及中点应力→配对internal_force→contact_rows/逐粒子kick(含约束reaction)
→速度/密度kick→半漂移。材料位于接触前，不把期望顺序写成实际顺序。
逐步另用同一冻结邻域/B和实际kick中点速度计算L_mid，不用它替换本构输入。
保存sum Fint·v_mid 与 sum V sigma:L_mid的伴随功、旧L本构功及两者差，
显式暴露时间层分裂。内力净力/净矩/所有3x3应力/连续密度全部记录。

首步独立答案：初始L/S/Fint=0，顶点壁速=(-20y,20x,0)，e=wall/|wall|，
Jt=.2*dt*e，v1=Jt/m，x1=x0+.5*dt*v1；其它点v1=0。
Fn向-z冲量与法向夹具reaction相反；Q=Jt·wall-|Jt|^2/(2m)，
W=Jt·wall，K=|Jt|^2/(2m)。首步参考不调用被测kick/几何/算子。
非仿射随机三维场核对伴随功；独立仿射A参考不调用被测函数生成。
没有此自由表面小块的连续体解析解，禁止把自算结果称独立精确轨迹。

## 事前容限、预算与失败规则

仅允许dt=5e-9/2.5e-9s到共同2e-7s(40/80步)。先2步成本探针，
声速c=sqrt((K+4G/3)/rho0)，dt<=.1*h/(c+vmax)；此为必要声学界，
不是纠正自由表面算子的谱稳定证明。实际NaN/密度/速度守卫仍生效。
每次外120s/内110s，预测有余量才启动；超限BLOCKED，禁止无限调参或重试。

独立首步：速度abs1e-12m/s、位置2e-17m、功1e-18J；仿射梯度1e-10/s。
伴随功abs1e-12W+1e-10*sum|Fint||v|；净力abs1e-12N+1e-10*sum|Fint|；
净矩abs1e-14Nm+1e-10*sum|x cross Fint|。这些是代数舍入界，不是实验误差。
每步接触/动能账abs1e-15J+1e-9*sum各功绝对值；摩擦负值仍用旧kick舍入界，
不得abs(Q)/裁零/残差热。cond(M)>1e6失败，不pinv。质量relative1e-12；
密度偏差>.01、速度>1m/s、位移>1e-5m、非有限立即拒绝/回滚。
固定DOF速度abs1e-12m/s、位置2e-17m；标签交换物理结果逐位一致；
排序/重启对齐ID后浮点atol1e-12(应力1e-6Pa)，质量/ID/标签精确。

全局K+Es+P+Q-Wwall-Wconstraint及其旧时间层分裂持续报告；这里Es含剪切与EOS体积能
诊断，P为累计塑性功(不是压力/法向弹簧能)；return算法耗散另记，不补热。
当前Fn/固定z夹具没有法向弹簧储能，不能虚构penalty势能。
若abs残差>1e-12J+.01*(|Wwall|+|Wconstraint|+Q)，明确ENERGY_WARNING，
不因完成即授能量认证。旧Jaumann闭循环FAIL不被此次小应变结果覆盖。
两档差异只比较共同终时和同一ID：position1e-7m、velocity1e-2m/s、
stress1e5Pa、density1kg/m3、累计Q/W差1e-8J；这是短程工程筛查，非收敛阶。
粗档必须完整且独立代数/输出读回通过才跑细档；能量WARN可保留作APPROXIMATE证据。
新checkpoint仅此固定组件允许重启，生产false；失败candidate=REJECTED不可重启。
所有旧FAIL不删改；三工况/热/供料/长时/空间加密NOT_RUN。
