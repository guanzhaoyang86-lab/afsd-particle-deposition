# Finite-strain free-energy variant：只做设计，不接生产

状态：DESIGN_ONLY / IMPLEMENTATION_NOT_RUN。本轮没有改material.py、G、EOS、
FKS参数或塑性热；不拿一个新本构覆盖当前失败。触发依据是实际材料更新器闭循环
四档细化后均为0.146111111111J，见COUPLING_FAILURE_REVIEW.md第4节。

## 为什么不能仅改Es

当前无旋转对角路径中，S=2G dev(logF)；同样的偏加载/卸载发生在V0和Jc*V0。
故闭循环偏功=(1-Jc)*G*(A:A)*e²*V0，与任意dt细化无关。
S和F闭合、当前Es归零，但功不零；完整反走路径会翻转这个功的符号。
给残差加一项热、改变Cp/G、放宽能量阈值或只重命名Es，都不修复这个性质。

## 一个明确独立的弹性候选（推导，不声称来自论文）

先做isothermal、无塑性的material-point实现，生产不开启。
状态F、J=detF>0、b=F F^T、bbar=J^(-2/3)b。
按参考体积的自由能[J/m^3]：
Psi(F)=G/2*(tr(bbar)-3)+K*(J-1-lnJ)。

导数给Cauchy应力：
sigma=G/J*dev(bbar)+K*(1-1/J)*I。
因此p=K*(1/J-1)=K*(rho/rho0-1)，保留当前显式线性rho压力闭合。
总弹性能E=V0*Psi，功率应满足dE/dt=V*sigma:D。
这是新的compressible isochoric neo-Hookean弹性变体，不是原Jaumann率的重命名。
G、K不缩放；F=I的小应变剪切/体积模量仍为G、K。

拟新增独立模块，不动旧update：
material_energy_variant.py:
- elastic_response(F, G, K) -> Cauchy stress, reference_energy_density, J；
- prescribed_path_step(state, L, dt) -> F更新及能量/应力诊断；
- J<=0或不有限直接失败，禁止pinv/应力或温度裁剪；
- provenance.variant="FINITE_STRAIN_ENERGY_DESIGN"，默认不可用于run_afsd_debug。

暂不设计一个伪完整的热塑性生产返回映射。
下一阶段若塑性扩展：明确F=Fe Fp、等体积塑性、与自由能共轭的应力度量、
FKS流动应力的率/温依赖及耗散不等式，保存Fp/必要硬化变量及checkpoint版本。
只有“把旧S替换成新sigma”不够；不能保留不共轭的旧塑性功公式而宣称一致。
FKS b1、参考率、温度材料参数等来源未决继续BLOCKED。

## 独立material-point tests（规格，全部NOT_RUN）

必须在实现前冻结容限。建议double精度基准量纲归一后：
代数参考1e-10相对+1e-12绝对；差分导数需做步长扫描，
不能挑选唯一有利步长。闭循环须显示dt细化趋势，不能只比末应力。

1. 纯旋转F=Q：J=1、Psi=0、sigma=0，参考Q用独立轴角式。
2. 非零预变形/预应力叠加刚转：F'=QF，sigma'=Qsigma Q^T，Psi不变。
3. 小应变模量：独立线弹性极限测G和K，不调用响应函数构造期望。
4. 纯体积F=lambda I：独立J=lambda^3；
   p=K*(1/J-1)，Psi=K*(J-1-lnJ)，与当前EOS原函数对应。
5. 简单剪切F=I+gamma ex tensor ey、J=1：
   Psi=G*gamma²/2；sigma_xy=G*gamma，
   sigma_xx=2G*gamma²/3，sigma_yy=sigma_zz=-G*gamma²/3。
   三个平面和反向剪切均检查。
6. 真正有限变形闭循环：复制本轮四段路径（不改几何），正向/反向，
   最终F/S闭合，积分sigma:D*Vdt随dt细化趋于0，
   同时与Delta(V0*Psi)比对；不得加energy residual heat。
7. 方向导数：独立中央差分能量 vs first Piola P=J*sigma*F^(-T)；
   另核对P:Fdot与J*sigma:D，非对称L含spin。
8. checkpoint/F状态、拒绝J<=0、随机旋转/随机可逆形变；
   不得把未实现塑性/热字段填零冒充已求解。
9. 塑性扩展才新增：detFp=1、独立标量屈服路径、sigma与Lp功共轭、
   非负物理塑性耗散、热分配/材料源与总能量闭合。当前NOT_RUN。

准入：独立模块通过这些基准并经设计确认后，才讨论接入SPH。
即使该弹性候选通过，也不解决边界patch、接触近零事件、工具外面、
论文材料标定或2us稳定性。当前生产本构保持不变。

