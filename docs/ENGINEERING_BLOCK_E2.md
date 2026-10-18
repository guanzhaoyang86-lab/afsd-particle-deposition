# E2：三维受载小块（2026-10-01）

## 运行前合同

ENGINEERING_APPROXIMATE / SMALL_STRAIN_REFERENCE_CONFIGURATION_BENCHMARK；非AFSD，非AA6061。
保留现有生产源码；脚本复用PySPH粒子/三维NNPS、HostLeapfrog、material.update、
corrected梯度及其成对伴随内力、上游checkpoint。不使用此前投影SBP原型。
参考几何固定算子是明确的小变形基准变体，不是有限变形生产solver.py复测。

立方体[0,L]^3，L=.01m，rho0=2700kg/m³，E=70GPa，nu=0，
G=35GPa、K=E/3；全部是制造参数，nu=0允许三维连续体有独立单轴解析解。
规定等温293.15K，不求热；测试律高屈服1e12Pa，必须实际ep=0，不冒充合金模型。
alpha=1e-5，c=sqrt(E/rho0)，k=pi/(2L)，omega=c*k，B=alpha*c²*k。
初始位移/速度/应力0。自t=0起恒定参考体加速度b_x=B*sin(k*X)，b_y=b_z=0。
唯一运动约束X=0面v_x=0（双向无摩擦刚性roller，非单边接触）；切向自由。
其余面外牵引0，其他全部自由度由力更新；不规定自由点运动，不加阻尼/稳定化。
体力只依原始坐标定义，不取自内力/参考解误差；反力记录、不做功。

独立线性三维连续体参考：
u_x=(B/omega²)*(1-cos(omega*t))*sin(kX)，
v_x=(B/omega)*sin(omega*t)*sin(kX)，
sigma_xx=E*alpha*(1-cos(omega*t))*cos(kX)，其余分量0。
R_x=-E*alpha*L²*(1-cos(omega*t))；体力合力=E*alpha*L²。
Wext=.5*E*alpha²*L³*(1-cos(omega*t))；
Kinetic=.25*E*alpha²*L³*sin(omega*t)^2；
Elastic=.25*E*alpha²*L³*(1-cos(omega*t))^2。
共同终时tf=pi/(2omega)=L/c≈1.964us：W=3.5e-6J，K=U=1.75e-6J，R=-70N。
这些值由连续体波动方程推导，不使用被测SPH算子生成。

粒子格点含物理边界；张量梯形参考体积，边界半权，m=rho0*V，
每档体积均L³、质量均rho0*L³。不是质量缩放，是同一域积分权重。
参考坐标上提取现有corrected梯度矩阵，原始NNPS、h=1.3dx、cond<=1e6。
每步真实位置由未改HostLeapfrog推进；参考算子不随微小位移移动。
材料在半步/整步几何增量上调用实际update；体积应力用线性K*trace(grad u)。
生产rho/EOS/thermal耦合未调用，输出只能标reference_density或linearized_density。
原Jaumann更新保留客观转动，和线性参考差是O(alpha²)的小量，但要实测能量误差。
应力/位置不同时间层必须明确；外功=体力·实际位移，约束反力功=0。

运行安排：先只做n=5的5步计时探针，输出PILOT_NOT_ACCEPTANCE。
名义dt=CFL*dx/c；另取参考弹性离散谱上界的显式稳定保护，不改质量/模量。
谱上界：nu=0时sym为正交投影，omega_max²<=c²||sqrt(V)D/sqrt(V)||_2²，
用||A||1*||A||inf上界；dt<=1/omega_bound（不超过中心差分稳定界的一半）。
随后时间序列n=5,CFL=.1/.05/.025；空间n=5/7/9，CFL=.05，均到tf。
预计名义40/80/160步，空间80/120/160步；真实步数取ceil(tf/dt)，不越终时。
依据pilot估墙钟；单命令timeout60s。预测超过60s不启动，标NOT_RUN；绝不启动>600s。
同一解的误差趋势只在相同终时比较，不以失败提前终时声称收敛。

预设容限（工程筛查，非论文标准，不随结果调整）：
- 代数/实现：质量相对1e-12、ID精确；梯度/力接口相对1e-11；力/参考坐标力矩归一1e-10。
- 材料/数值：ep、塑性功、返回耗散严格0；NaN失败；cond(M)<=1e6；最大位移/L<1e-4。
- 实际累计功与K+U差/Wref<=1%；支承约束位移为0，支承功为0。
- 最细位移体积加权L2相对<=10%，全应力张量L2<=20%，反力相对<=20%。
- 侧向寄生位移L2/独立轴向位移L2<=5%。这些是宽于舍入的预先工程门槛。
- 每次空间细化位移/应力误差至少下降10%；否则趋势FAIL，不能仅凭误差较小通过。
- n5相邻dt解差在最细<1%，细化差比<=.4（预期二阶）；若无可辨差<=1e-10注明舍入。
容限用于暴露受载误差，不用于将原常应力节点0.326残差改判PASS。

输出真实粒子ID、参考/当前坐标、速度、应力、ep、载荷/反力、能量时序、checkpoint，
参考字段用reference前缀；未求热量/接触/有限应变密度不伪造。初始及末态可写小型VTK。
PASS只能覆盖本小变形受载问题；塑性接触、热耦合、真实AFSD、GPU和生产均NOT_RUN。

## 实际调用与接口范围（数值序列前核验）

HostLeapfrog.advance -> 固定上游LeapFrogIntegrator.one_timestep：
stage1半步位置漂移 -> Block.evaluate -> sync(半步材料状态) ->
实际伴随internal_force + 独立体力 -> 法向roller约束/反力冲量 ->
stage2速度kick及第二次半步位置漂移 -> sync(整步材料状态) -> 记录。
生产SolidSolver.evaluate没有调用；新基准不证明生产应力/热学分裂正确。
提取矩阵只调用probe_boundary_sbp.extract_original，未调用其project/SBP构造。
独立接口测试检查该矩阵与实际corrected_operators、配对internal_force一致。

报告区分actual kick反力/累计冲量与整步瞬时反力；解析参考只在对照使用。
张量梯形的体力合力参考为70*q*cot(q)N，q=pi/[4(n-1)]；连续体参考70N。
不将离散体力重新归一到70N；其积分误差属于空间近似误差，必须报告。
初始自由横向平移/绕法向刚转模态未额外钉住，检查实际寄生横向运动。
checkpoint只作诊断快照，restart_allowed=false；重启此变体NOT_IMPLEMENTED。

## 实际结果：共同终时完成，但工程应力验收FAIL

全部5档真实力驱动小块到达同一tf=1.9639610121239317e-6s。
没有物理告警中断；每档exit1是误差验收FAIL，不是未跑完或NaN。
没有调整几何、材料、载荷、容限、h/dx或加阻尼来追求通过。
本轮src/、旧配置/测试/失败证据、固定PySPH上游均未改。

| 粒子数/CFL | 实际步数 | dt/ns | 位移L2误差 | 全应力L2误差 | 反力误差 | 能量账残差/Wref | 推进墙钟/s |
|---|---:|---:|---:|---:|---:|---:|---:|
| 125/.1 | 41 | 47.9015 | 12.359% | 28.659% | 19.138% | 0.01728% | 0.127 |
| 125/.05 | 81 | 24.2464 | 12.356% | 28.656% | 19.144% | 0.004425% | 0.252 |
| 125/.025 | 161 | 12.1985 | 12.355% | 28.655% | 19.146% | 0.001120% | 0.507 |
| 343/.05 | 121 | 16.2311 | 10.195% | 30.708% | 12.561% | 0.002059% | 1.162 |
| 729/.05 | 161 | 12.1985 | 7.620% | 27.530% | 8.863% | 0.001177% | 5.301 |

表中是对独立制造连续体解的数值误差，不是实验/论文误差。
ceil浮点向上得到41/81/121/161步，报告准确dt；没有缩放物理时间。
时间差异比较：相邻解位移差/参考范数6.33248e-5、1.62887e-5，比例0.257224，
时间敏感性PASS。不能把这一点扩大为n9全部模态的独立时间收敛证明。
空间位移误差比.8251/.7475通过预设趋势；应力比1.0716/.8965不通过。
最细应力27.53%>20%门槛，整体FAIL_ENGINEERING_TARGETS；不放宽标准。
粗网格也保留逐档误差FAIL；允许完成事先规定的诊断细化，不代表其物理验收通过。

所有档质量/ID、自由度约束、原始ep=0、力/力矩、能量账、有限值检查通过。
最细侧向寄生位移1.037%<5%门槛；约束功严格0。
独立全局冲量读回最大相对缺口5.99e-16。
最细体力功3.552117e-6J与连续参考3.5e-6J仍差1.489%；
动能1.915712e-6J、弹性诊断1.636446e-6J。
“能量账闭合”不等于“应力/动能等场与连续体一致”，这里已实际区分。

实际数组内存（不是进程峰值）：125/343/729粒子状态数组35500/97412/207036字节，
参考梯度375000/2823576/12754584字节；最后一档setup0.0443s、推进5.301s。
本基准没有接触/热学/更新NNPS几何，禁止把这速度外推论文尺度AFSD成本。
进程峰值RAM NOT_MEASURED；无GPU/并行/生产任务。

### 误差定位（事后只读诊断，不是改验收）

证据logs/engineering_e2_20261001/artifact_check_final.log：
n9的95.927%应力平方误差来自sigma_xx；
X=0、1.25、2.5mm三层合计占92.739%全域应力平方误差。
支承面加权sigma_xx=1.066407MPa，独立参考0.700000MPa；
下一层0.479959MPa，参考0.686550MPa，不是统一材料倍率偏差。
这支持支承邻域/离散受载响应问题，但尚未证明某一边界处理就是唯一根因。

证据等级：
- CONFIRMED：n5时间截断不是约28.66%应力偏差的主导来源；三档dt结果几乎不变。
- CONFIRMED：受力、约束、质量/ID、能量账与工况时间均真实；应力精度尚未过关。
- SUPPORTED：主要空间误差集中支承附近的轴向应力，符合边界/低能模态响应嫌疑。
- POSSIBLE：梯度应力恢复误差、边界力—体积积分不相容及弱刚度模态分别贡献，未分离。
- NOT_RUN：n9额外dt细化、塑性接触、有限变形、热学、新AFSD、GPU/生产。
旧Jaumann变体积闭循环问题仍在，E2低应变等温通过部分检查不解决它。

### 实际测试与证据

- tests/review/test_engineering_block_e2.py：7 PASS，0.041s。
  独立解析数值、质量/载荷、真实3D梯度、伴随功/力矩、被动标签、
  roller切向自由、checkpoint读回ID且禁止不完整重启。
- test_engineering_material_e1.py：6 PASS；test_spatial_probes.py：4 PASS。
- test_spatial_probes.py --boundary-red：原FAIL(.326157)保留且复跑。
- test_coupling_benchmarks.py --cycle：诊断PASS确认旧0.146111J问题，不是物理验收PASS。
- 5次共同终时小块：全部完成/误差FAIL；sensitivity.json：时间PASS、空间/整体FAIL。
- check_engineering_block_outputs.py：5组NPZ/VTK/checkpoint/源hash/ID读回PASS。
- 完整耦合回归、E2完整重启：NOT_RUN/NOT_IMPLEMENTED；没有静默认证。

首次5步pilot因NumPy bool JSON写入失败；日志和已写原始文件全部保留。
仅新脚本序列化改为Python bool，未改物理；重跑pilot成功。
之后加固定制造参数guard与实际中点反力/动量字段再运行共同终时序列。
两个旧脚本快照对应首次失败和最终pilot的源版本。
容限配置/预注册原文在logs/engineering_e2_20261001/，没有后改。

## 下一最小工作单元

先在同一小块上分离误差来源，不进入AFSD：
1. 将独立解析位移采样送入现有gradient，仅测应力恢复误差；不推进粒子、不调参数。
2. 以独立解析应力/加速度检查实际成对内力+体力的节点动力学残差，
   逐层报告支承附近与内部差异；约束反力单独去除。
3. 检查sigma_actual-E*sym(D*u_actual)，独立分离材料/客观率更新偏差。
这三个离线诊断使用原n5/7/9状态，按固定物理带X<=L/4和余域报告原始节点残差；
共同终时独立a_ref=0。不新增时间算例、不滤掉振荡节点、不调参数。
分解后才决定是否研究低刚度模态或有出处的边界改法；禁止重启通用SBP求解器或同时改核/材料/接触。
即使开源近似路线也不能把27.53%的应力错误用“能跑/守恒”掩盖。
核心实现路径仍mechanics.py::corrected_operators/internal_force、integrator.py、
material.py::update；本轮只调用，下一轮未获证据前不修改生产源。
详细命令见logs/engineering_e2_20261001/COMMANDS.md；冻结ENGINEERING_E2_FREEZE.json。

