# 独立审查补充（2026-09-24）

当前APPROXIMATE变体采用完整性修正梯度及伴随功共轭力，非原Table1逐字公式；
原式刚体转动和异密度功反例失败。执行器厚度固定0.5mm，接触单位面积罚刚度
固定2.63e12 N/m³；细化仅改变离散dx/h，不改物理BC。见VALIDATION_CONTRACT。
离散热改实际kick中点功，负耗散拒绝；多力耦合失败、本构储能不一致、外圆/顶面
接触缺口均未解决。所有原材料/启动假设继续保留，不因unit PASS消失。
Fig12初始图7.63mm与棒宽9.53mm差异UNRESOLVED；不据此改几何。
Fig8精确ROI/切面和图像原始数据缺失，count/volume不可冒充实验面积。

# ASSUMPTIONS / UNRESOLVED

实施更新（2026-09-24）：所有论文级缺项仍保留，不能由DEBUG使用值推断已经解决。
可运行变体的D01–D10、PDF第5/6/14页图像核验、P40/P41访问阻塞、P44已取得内容及
EOS/接触/FKS偏离集中在PARAMETER_SOURCES.md；实测范围见PHYSICS_STATUS.md。
A01执行解释为p4的300rpm与三档轴向进给；原p6措辞矛盾没有删掉。
A11已固定本实现+z旋转/+x走速、AS y<0和仅平移工具系，论文画面投影仍待配准。
DEBUG不静默升级为APPROXIMATE/PAPER_COMPARISON；未完成宏观定性对照或定量验证。

状态约定：UNRESOLVED / PAPER_CONFIRMED / REFERENCE_CONFIRMED / DERIVED / ASSUMED / CALIBRATED。
每次解决需补：采用值、单位、来源页/式/源码位置、理由、适用配置、敏感性测试。没有数据不是将值设为零的依据。

| ID | 事项 | 初始状态 | 操作 |
|---|---|---|---|
| A01 | p.6 把三档 mm/s 写成 rotational velocities | UNRESOLVED（有明确候选解释） | 用 p.4 300 rpm + 三档进给作配置依据，并留下冲突记录 |
| A02 | FKS b1 的 MPa 表头与 Eq.6 量纲 | UNRESOLVED | 查 P41；原始数值12.44保留；禁止盲目乘10^6 |
| A03 | Troom、Tmelt、参考应变率及 FKS 正则化 | UNRESOLVED | 查原始模型与合金适用范围；测试边界 |
| A04 | k(T)、Cp(T) | UNRESOLVED | 查 P40；常数近似单列，不能称论文原值 |
| A05 | 几何外径、基板尺寸、棒的模拟长度、孔隙 | UNRESOLVED | 查 P40/补充材料；图估计单列 |
| A06 | 起始加热/驻留/抬升/停止过程、路径总长 | UNRESOLVED | 区分实验描述与模拟初始条件；短时间demo不等于整段沉积 |
| A07 | 摩擦律/参数、接触导热、功热分配及散热 | UNRESOLVED | 查 P40/P44，并分别列出上游默认和论文设置 |
| A08 | EOS、h/间距、稳定化、CFL、质量或时间缩放 | UNRESOLVED | 记录来源与所有数值加速；先做未缩放小测试 |
| A09 | 标签带到五分析区/完整沉积的采样协议 | UNRESOLVED | 对照 Fig.5/6/10/11/12；局部与全棒标签分开 |
| A10 | 原始实验数据和误差棒数值 | UNRESOLVED | 论文发表时不共享；读图值标为图数字化，不虚构原始数据 |
| A11 | AS/RS、工具系与实验室系映射 | UNRESOLVED（实施约定） | 在配置中固定，再做旋向反转对照；不能追着结果改标签 |
| A12 | 3D 自由表面/氧化物面积比的测量算子 | UNRESOLVED | 点数、体积与截面积指标严格分名；ROI须可重跑 |

允许在显式命名的 approximate/debug 配置中使用有依据的假设推进工程，但论文对照配置不得静默继承这些假设。模型审查者须决定哪些缺失阻断定量复现，哪些可经敏感性分析支持定性比较。

