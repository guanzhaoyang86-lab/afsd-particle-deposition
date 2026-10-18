# VALIDATION_CONTRACT v1 — independent review pre-registration

2026-09-30方法类验收补充见CONTACT_CLASS_CONTRACT.md/METHOD_CLASS_REVIEW.md。
固定坐标、物理参数、ROI、N01–N10均不改。多解接触按独立KKT/合力/功验收，
不强认任意面力分配唯一；舍入预算变更及失败证据单列，原始负热不得裁零。
局部方法认证不等于全耦合或论文验证；新METHOD_CLASS_FREEZE.json不覆盖旧冻结。

2026-09-30授权的首发诊断/最小接触变体修订见COUPLING_DIAGNOSTIC_CONTRACT.md，
结果见COUPLING_FAILURE_REVIEW.md。本文件坐标、物理配置、ROI和原容限不变；
新变体必须显式选event_aware_v1并标APPROXIMATE/BLOCKED，不能用可视化切换。
旧SCIENTIFIC_FREEZE.json保留为旧版本证据，新修订另有COUPLING_REVIEW_FREEZE.json。

冻结日期：2026-09-24。本轮在运行新数值测试之前写入。原始副本及SHA256保存于
logs/review/contract_v1_pretest.md、contract_v1_pretest.sha256；后续补充只增加说明，
不得按测试结果放宽本表门限或删除失败。阈值是数值审查门槛，非论文实验误差。

## 1. 预先冻结的数值验收

| 编号 | 检查 / 独立参考 | PASS门槛与依据 |
|---|---|---|
| N01 | 3D核积分 / 分段多项式解析积分；独立导数 | 归一化绝对1e-7，导数相对1e-6；数值求积/差分误差，不是材料拟合 |
| N02 | 随机非共面邻域 vs 暴力欧氏距离；z向/剪切自由度 | 邻居集合完全相等；三轴均参与，禁止厚度复制 |
| N03 | 非零初始偏应力刚体旋转，独立轴角矩阵或解析sin/cos | 应力/不变量相对1e-10，零塑性；不可调用被测rotation作oracle |
| N04 | FKS独立math标量式、三个剪切面/组合J2及独立根求解 | 标量应力/返回残差相对1e-9；允许绝对1e-4Pa的舍入下限；温域/非法值必须失败 |
| N05 | 独立两体导热及塑性功，完整对称张量缩并 | 两侧热守恒相对1e-12；解析值相对1e-9（离散时间衰减另允许1e-3）；热流从热到冷 |
| N06 | 接触一步功-热与实际速度增量、静/动壁符号；深穿透检测 | 独立离散能差相对1e-10+绝对1e-12J；法向/滑动符号严格；深陷刚体不得被无声漏掉 |
| N07 | 物理状态排序/标签交换与移除/完整重启 | ID映射完全一致；标签变化同排序时物理字段逐位相同；重排还原后相对1e-12；代码hash变化禁止旧重启 |
| N08 | 小耦合安全/收支 | 质量1e-12相对；有限值；密度偏差<=5%；穿透/dx<=10%；瞬时速度<=20m/s；残差<=1%交换能尺度+1e-7J。较旧5%门槛收紧，仅作为审查判定，不靠改运行容限掩盖失败 |
| N09 | 同一物理终时dt、dt/2、dt/4 | 相邻两个最细解：动能/总塑性功/峰值速度/体积平均ep差<=5%，峰值ep差<=10%；小分母用固定绝对下限ep=1e-8、能量=1e-10J、速度=1e-6m/s。误差应呈下降趋势；只有三点且未渐近则NOT_CONVERGED |
| N10 | 同物理几何/边界的dx=.5mm与.25mm | 上述整体量差<=10%、峰值ep<=20%；两网格只能检验敏感性，不能授予渐近空间收敛。绝对下限同N09；变化超限FAIL，不反复调参 |

独立检查必须以解析式/手算/不同算法生成参考，不能以被测模块回填答案。
旧测试会重跑但只作回归，旧日志不计本轮证据。不能承担的检查写NOT_RUN。

## 2. 物理、坐标与敏感性冻结

仅PySPH 1cd3dea9cfaac7d9f1bda2a2491ce4b42a31dbe5；保留原代码/日志。
论文共同量：300rpm=10pi rad/s，+x走速.00212m/s，层高.0015m；进给
.00106/.00212/.00424m/s沿-z。P0 p6误称rotational velocities而使用mm/s，采用p4量纲一致定义。
实验室z=0为初始基板上表面，+z朝工具；omega沿+z、俯视逆时针。
AS由(omega cross r).Vtool>0确定，当前为y<0；RS为y>0。工具参考系只平移，
x'=x-Vtool*t，y'=y，z'=z；共转方孔检测坐标不得拿来代替轨迹参考系。
可视化只允许刚性视图变换/显示单位转换，不得换AS/RS、改粒子位置/温度/力学参数。

本轮定量数值敏感性只用明确APPROXIMATE/DEBUG物理变体，不升级PAPER_COMPARISON。
继承PARAMETER_SOURCES D01–D10物理假设，除已证明错误的最小修复外不调材料或摩擦
以追求任何论文比值。b1=12.44无量纲是近似解释，Table2 MPa矛盾未解决；
Troom/Tmelt/rate0/低率floor、常数k/Cp、预热/过盈/预结合界面仍为假设。

同终时目标先定为2.0e-6s；时间三档CFL=.1/.05/.025，固定dx=.5mm；空间对照
dx=.5/.25mm，固定所有长宽厚、孔径、过盈、温度、速度、摩擦系数、材料、物理执行器
厚度(.5mm)、接触单位面积罚刚度。必须记录h随dx的正常离散变化，不把罚刚度或
执行器厚度随dx改变误称“仅网格细化”。若旧实现边界绑定dx，先记录并修复该缺陷。
每次最多600s，超时/告警保留不完整结果，不用不同终时冒充敏感性比较。

## 3. 实验测量与后处理冻结（页图复核完成，不产生参考数据）

Fig8的实验量是氧化物截面面积比A_AS/A_RS，不是3D点数比或体积比。
以工具轴线为中线，两侧横向分析带各15.88mm；实验五区沿程边界未知，不等分后称原文。
Region1启动、Region5终止，Regions2–4仅为原文名义稳态区；数值稳态需另给判定，不以
时间很长或比值接近论文判定。精确沿程端点/切面方向/切片厚度在取得来源前为null/BLOCKED。

可执行近似后处理必须提供：截面原点、法向、厚度、沿程ROI端点、工具轴线位置、
AS/RS分析带、材料/标记选择、起止排除区及时间。缺任一关键量拒绝“实验面积比”。
指定分辨率的SPH体积分配/截面重建面积算子可作APPROXIMATE，需独立几何金标准、
切片厚度/像素敏感性；未实现时标NOT_RUN。R_count、R_volume单列，禁止重命名为R_area。
分母0返回null/UNDEFINED，不加epsilon；没有标签或局部带未到ROI不等于实验零氧化物。

tracer_type=初始z0在3–4mm的局部外表面/中心代理；wholebar_tag仅额外全棒诊断。
二者必须分别选择且写进测量结果，不能拼接/互换来补齐数据；均AA6061-T6，不入受力。
Fig10弯曲富集/中心偏AS、Fig8外层偏RS、Fig11/12轨迹变化都只是待检验假设。
轨迹ID只按初始几何/固定随机种子选择，记录全部选中ID；不得事后挑“好看”的路径。

校准/验证集合要在调参前写明。使用某图校准后，该图不能再声称独立验证。
本轮没有实验校准、没有捏造CSV、没有图数字化比值。若以后读图，保存原图hash、
坐标标定/像素点、读数者、重复读数和线宽/分辨率导致的区间；误差栏不得虚构为实验误差。

## 4. 必须保存的结果与冻结机制

结果含classification、status、contract_version/hash、config/source/upstream hashes、
exact argv/exit、物理终时、实际步数/dt范围、峰值进程内存、步耗时、质量/外功/热流/
总能量/数值耗散、密度/温度/应力/塑性/穿透极值、NaN及任何报警。字段未求解则缺省并说明。
敏感性比较同时保存绝对值、归一化分母及相对差，不只报PASS。
后续可视化必须核对冻结config/source签名，只能读取副本；任何物理改动须新variant、
新契约版本和重新验收，不得覆盖本轮证据。

## 5. 页图复核补充（不改第1节预注册容限）

本轮重新渲染并人工查看P0 PDF第5、6、7、9、10、11、12、14页完整图像，
见logs/review/paper/page_{5,6,7,9,10,11,12,14}.jpg。P0 SHA256：
0d82cd71c53f1005931581750b5a3eff582f6e5a65c5f61e883bbef15d63c119。

- 第5页核为15/(62*pi*h^3)乘分段q^3-6q+6或(2-q)^3，支持半径2h；
  Table1包含3D D/spin/Cauchy应力。原始体积加权梯度与rho^-2应力力在不同密度时
  并非离散功共轭。本轮启用的完整性修正/伴随力是显式APPROXIMATE方法偏离，
  不声称它逐字实现Table1。原公式仍保留供审查，不作生产路径。
- 第6页FKS热项为1-(1+exp(-c1*Tstar))^(-1/c2)，1/c2作用于完整括号。
  第14页b1确实标MPa，而第6页Eq6要求无量纲；目前不能判定作者勘误。
  b1=12.44无量纲和参考温度/率仍ASSUMED，绝不机械乘1e6。
- 第7页Fig6：相对于工具轴线各侧15.88mm；ImageJ阈值面积范围1–1000 µm²。
  缺切面位置/方向、五区精确沿程边界和图像分割细节，不能用等分区间冒充原文。
  该阈值不能直接应用于粒子体积或数目。实验图像和XRCT原始数据未提供。
- 第9页Fig8的方/圆实验点与模拟线不是同一算子；有局部实验点超过1，不能以
  “所有区必须外层偏RS”验收。Region1/5为开始/结束，2–4为名义稳态，非精确时间窗。
- 第10页Fig10的铜示踪/弯曲富集和结束聚集只作待检验形貌；图标沉积中心不能
  替换工具轴线重新居中。第11/12页仅展示所选轨迹，不能作为全群体统计。
- Fig11/12的显示映射可取(horizontal,vertical)=(y_lab,-(x_lab-V*t))，
  因而AS在左、相对材料沿图上方流动；这是刚性视图旋转，不交换AS/RS或改变omega。
  Fig11初始尺寸标9.53mm，Fig12标7.63mm，后者与论文棒宽有疑点，保持UNRESOLVED，
  不据此改物理几何。本轮没有图像读数或误差百分比。

## 6. 可执行后处理与审查冻结文件

`afsd_validation.measurement.measure_proxies(points, scalar_fields, spec)`只读现有
VTK读回数据。spec的字段定义和严格校验在src/afsd_validation/measurement.py；
所有几何量SI，所有字段必须显式给出（包括ROI来源与ASSUMED/DERIVED/REFERENCE_CONFIRMED
状态），禁止null、自动补零、混合标签或事后重选中心线。无真实ROI时不执行论文测量。

算子：粒子中心位于半开x区间[x0,x1)、闭合切片|n·(x-origin)|<=thickness/2、
|y|<=.01588m，且material_id=1、显式单一标记值/范围时才选入。
(-omega*y)*Vx>0归AS，<0归RS，y=0单独排除并报告。
R_count=N_AS/N_RS；R_volume=sum_AS(m/rho)/sum_RS(m/rho)。两者均非实验面积，
experimental_area_ratio始终null/BLOCKED；即使有粒子，STEADY_ASSUMED也不变成已证稳态。
结果含完整spec及其SHA256、选择的全部ID、原始分子/分母、分类、单位与未定义原因。
9项独立几何测试见tests/review/test_measurement_review.py（不是仿真/实验验证）。

审查物理配置冻结在configs/debug/starved.json及configs/review/{dt_half,dt_quarter,space_half}.json；
resolved_config.json是实际运行的完整参数/来源/状态快照。各运行source_hashes.json
冻结实际求解源码；改源禁止旧checkpoint直接续跑。docs/SCIENTIFIC_FREEZE.json
另固定本轮源码、配置、协议及关键证据hash，运行scripts/verify_scientific_freeze.py
检验。该文件不是防恶意篡改的外部签名；更改物理必须显式新版本重审，不能重写hash
掩盖变化。后处理/可视化先校验冻结，再校验series中VTK SHA；不写入原始输出目录。

本轮参数不作实验校准，Fig8/10/11/12全部保留为未来验证集合；材料/热学/接触假设
见PARAMETER_SOURCES及SCIENTIFIC_REVIEW。后续拟用它们校准须先改变用途记录，再另留
独立验证数据。当前没有获得生产授权。
