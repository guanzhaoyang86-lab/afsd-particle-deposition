# 方程核验、参数来源与 DEBUG 偏离（2026-09-24）

本文件定义可执行的 HOST_FKS_DIMENSIONLESS_B1_HOTSTART_SCALED_GEOMETRY 变体，
不是作者程序，也不是论文参数已经补齐的声明。所有计算 SI，三档 case 的具体值、单位、
来源、状态见 configs/debug/*.json；configs/cases/* 的缺项仍保持 null。

## 原文图像核验

已实际查看用户原 PDF **第5、6、14页完整渲染图像**，不是只读提取文本。
证据：logs/paper_verification/P0_page_{5,6,14}.jpg（另保留 PNG）。
原文 SHA256：0d82cd71c53f1005931581750b5a3eff582f6e5a65c5f61e883bbef15d63c119。
渲染命令：pdftoppm -f PAGE -l PAGE -r 140 -jpeg -jpegopt quality=88 -singlefile papers/sph_simulation.pdf logs/paper_verification/P0_page_PAGE。

- p5 Eqs2–3：三维系数 15/(62*pi*h^3)，q=|xi-xj|/h；q<1 时 q^3-6q+6，1<=q<2 时 (2-q)^3。原印刷矢量距离表达按范数解释。自项不求方向；不同粒子重合直接失败。没有改成 cubic spline。
- p5 Table1：rho_i sum(m_j/rho_j)(vi-vj).gradW；对称 Cauchy 应力动量；sigma=S-pI；三维 D、spin 和 Jaumann；导热谐均值 4ki*kj/(ki+kj)。实际离散路径在 mechanics.py/thermal.py。
- p6 Eqs4–8：H=a1+a2*atan(a3*ep)；Lambda=1+b1*Tstar^b2*ln(rate/rate0)；Theta=1-(1+exp(-c1*Tstar))^(-1/c2)。**1/c2 作用于整个分母括号**，不是仅对 exp 或 Tstar。
- p14 Table2：a1=317.6 MPa，a2=1 MPa，a3=70，b1=12.44（表列 MPa），b2=5.65，c1=13，c2=0.01。Eq6 要求 b1 无量纲，矛盾未被消除。仅 D08 变体用无量纲12.44，绝不乘10^6。
- p6 的三档 mm/s 被称为 rotational velocities；按 p4 的300 rpm旋转和1.06/2.12/4.24 mm/s轴向进给实施。矛盾保留，不称作者已勘误。

## 参考文献取得情况与采用范围

| 源 | 实际取得/访问状态 | 可采用内容 / 未解决项 |
|---|---|---|
| P40, Stubblefield et al., Materials & Design 202:109514, [DOI](https://doi.org/10.1016/j.matdes.2021.109514) | DOI落地HTML已存；出版商全文页403；本次未取得可核验全文 | k(T)、Cp(T)、精确几何/启动程序仍 UNRESOLVED，不从摘要编造 |
| P41, Stubblefield et al., JMPT 305:117593, [DOI](https://doi.org/10.1016/j.jmatprotec.2022.117593) | 检索到出版商摘要/开放稿入口文字；打开页面403；未取得全文 | P0这套FKS标定的参考率/温度和b1矛盾仍 UNRESOLVED |
| P44, Fraser (2017), [UQAC原始学位论文](https://constellation.uqac.ca/id/eprint/4246/1/Fraser_uqac_0862D_10345.pdf) | 已取得359页原PDF，papers/references/P44_Fraser_2017.pdf；选读PDF页87–99、111–144 | 线性EOS、塑性功、摩擦热/分配、FKS形式和三角面接触依据；不能把其旧参数表当P0新标定 |

P44 SHA256：b30565aa83862b1474ad819b90dbdc24d48d595615102344565a7d9cc13fc9d3。
提取可重跑：python -X utf8 scripts/inspect_reference_pdf.py papers/references/P44_Fraser_2017.pdf 87:99,111:144 --output logs/paper_verification/P44_selected.txt。
保留的 P40/P41_landing.html **不是论文全文**。访问失败是本环境本次结果，不断言文献永久不可获得。

P44逐项：

- PDF92页 Eq3-37/38：p=c²(rho-rho0)，文中c式含 K+4G/3。本 DEBUG 采用体积模量 K=2G(1+nu)/(3(1-2nu)) 闭合 p=K(rho/rho0-1)，波速限步另用 sqrt((K+4G/3)/rho0)。这是显式热弹塑小体积变形选择，**不是逐字等同P44的EOS系数**，仍需波速/压力敏感性核验。
- PDF112页 Eq3-106：塑性功热 chi*sigma:deps_p，讨论 chi=1；本变体 chi=1，q_new*dep 是隐式返回的塑性功，1.5G*dep² 数值返回耗散单列，不偷偷变成热。
- PDF112–113页 Eqs3-109–111：摩擦相对切向功和两侧热分配；分配基于 sqrt(k*Cp*rho) 的热感受率。DEBUG eta=0.5 只对应等热感受率假设；工具真实材料参数尚缺。
- PDF124–128页 Eqs3-127–133/Table3-2：原FKS含额外 b3，旧表 b1 无量纲、值1.7，参数与P0不同；没有替换P0参数。提取 Eq3-129 的分子与该文“室温为零”描述相矛盾，代码采用已看图的P0室温归一化式，不复制可疑旧式。
- PDF129起 §3.11：刚性零厚度三角面、SPH从面接触及罚力。后续 Eq3-160 含 stick/slip 与剪切屈服限制。本轮**只复用PySPH平壁罚力并组装有限解析面**；完整三角面/粘滑/角边一致性仍需开发验证。

## 可运行 DEBUG 假设表

| ID | 采用内容 | 来源/适用性及下一步 |
|---|---|---|
| D01 | 方棒2.5×2.5×4.5mm；基板6×6×1.5mm；肩半径2.7mm；工具顶z=4mm；dx=.5mm | 明确缩小几何，不是9.53mm论文棒。层高仍1.5mm=3dx，最低调试解析要求，不代表网格收敛；432基板+225棒=657可变形粒子 |
| D02 | h=1.3dx；无人工粘性/应力/XSPH/核修正 | 稳定化未校准，只允许微秒探针。仿射内部误差容限2.5%，边界完备性尚缺 |
| D03 | 孔宽2.496mm，即每面2微米初始过盈 | 主动设置可测接触载荷，初始罚弹簧能计入总能；不是预设金属运动。棒/基板从t0共享SPH邻域，相当于预结合界面；没有预测冶金结合或分离 |
| D04 | 所有金属初温723.15K、工具恒温733.15K、底座703.15K；直接以300rpm启动短探针；无停留/抬升历史 | 人工热启动，不是已计算预热。dwell=0/lift=0是已声明程序选择，不是缺失值补零；非零设置拒绝。到步数/墙钟上限即终止 |
| D05 | k=167 W/(m K)，Cp=900 J/(kg K)，常数 | DEBUG假设，不是P40温变表；G/nu取P0常数。无热膨胀/相变/热弹耦合源，温度仍影响FKS屈服 |
| D06 | kn=.05*G*dx；mu=.3；法向无阻尼；切向min(mu*Fn,m*|vt|/dt) | 实际调用上游RigidBodyWallCollision.loop，en=1避开其硬编码有效质量阻尼；每步无记忆滑动正则化，**不声称静摩擦粘着**；有限面角点可双约束，未做完整边缘接触 |
| D07 | 接触h_c=1000 W/(m² K)，面积dx²/活跃面，eta=.5；自由面绝热 | 恒温刚体为无限热库，散/入热单独计账；不解工具/底座温度。未计对流/辐射是显式DEBUG热边界，不能用于论文热稳态 |
| D08 | b1无量纲12.44；Troom293.15K、Tmelt855.15K、rate0=1/s；rate_eff=max(rate,1/s) | 未解析P41；该floor为可见正则化且累计计数，不是静默裁剪屈服/温度。rate=sqrt(2Ddev:Ddev/3)，不是以期望结果选率；越温域/非正屈服立即失败 |
| D09 | CFL=.1，dt-half=.05；波速、罚接触和热显式限步取最小；全量真实质量 | 未缩放时间/密度/模量。复用上游drift-kick-drift，热/材料分裂一阶；非整套二阶热塑耦合。已做dt敏感性，尚未收敛 |
| D10 | 密度5%；穿透/dx10%；速度20m/s；能量5%相对+1e-7J绝对；质量1e-12相对；所有NaN/Inf失败 | 调试告警门限，不是实验误差容限。相对能量尺度=初始机械/接触能或累计外部交换绝对量之较大者，绝不用巨大初温热能掩盖误差 |

## 运动、供料、示踪与坐标

初始自由粒子速度为零。仅方棒最上层执行器规定 x 走速和 z 轴向进给，y自由；
没有强制整棒随工具自转，也没有对自由粒子施加“旋转+径向+平移”速度场。
旋转由方孔接触传递；工具/底座刚性，基板和棒均有3D应力/密度/温度更新。
执行器反力功按约束冲量乘该自由度中点速度累计。

有限棒在封闭材料集合内被压入；外部新增质量/焓严格为零（已知边界，不是缺值）。
通过z=层高平面的首次进给质量按持久ID/标记累计；无限长入口尚未开启。
append_feedstock仅已测试身份接口，不代表已实现守恒供料边界。

外表面tag=1、中心代理tag=2均按初始位置保存；tracer_type仅在z0=3–4mm局部带非零，
wholebar_tag是额外全棒诊断。全部material_id=1、AA6061-T6；不解铜/氧化层异材。
中心代理半径dx/2是缩小几何调试选择，不声称解析论文1mm铜丝。
实验室：+x走速，+z基板法向，omega沿+z（俯视逆时针）。仅平移工具系x'=x-Xtool，
不随工具转动；方孔检测另用共转局部坐标。AS由(-omega*y)*Vx>0判定，当前AS为y<0。
从没用AS/RS或示踪标签更改受力；反转旋向测试通过。具体论文画面投影仍未定量配准。
