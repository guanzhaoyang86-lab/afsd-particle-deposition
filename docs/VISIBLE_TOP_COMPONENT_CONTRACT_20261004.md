# TOP 接触资格扩展的先行组件合同

本合同在组件运行前登记。只验证现有 ChainRuleDiscreteContact / BatchedChainRuleContact
在有限 TOP 面及 TOP–BORE 棱边的数值实现；不修改物理公式、几何或任何旧冻结文件。
不运行真实 AFSD 候选，不增加接受时间，也不替代同一真实父态的 coarse/half 事件比较。

已存在的 TOP 几何及静态势能梯度测试不等于 TOP 有限步摩擦、创建/释放和旋转链功
已经覆盖。新增 fixtures 使用 R=2、孔宽=1、H=1、粒子半径=0.1、kn=100、kt=200、
mu=0.3 的明确制造单位制，并附 1e-3 几何尺度检查。它们不是实际生产尺寸模型。

独立参考：平 TOP 距离 d=z-H，势能 U=kn/2 max(r-d,0)^2；TOP–X+ 孔口棱边
TOP 距离为 sqrt((x-a)^2+(z-H)^2)，孔壁距离为
sqrt((x-a)^2+max(z-H,0)^2)（夹具取 x<a、y 严格在孔面内、靠近上口）。
重复面势能按现有模型显式相加，不改成最近单面。法向有限步相对功等于 -ΔU；
滑动摩擦用独立一维 Coulomb/BE 弹簧算式核验，算法耗散和弹簧释放能单列，不并入热。

必要范围：TOP 压缩、创建/释放、粘着/滑移、TOP+BORE 入口/出口、棱边法势能链功、
有限旋转/平移与任意共同三维坐标旋转、完整 scalar/batch 接触及几何字段、输入/history
纯性、开放孔口。另核验 TOP 两端不活跃而中间接触的实际世界弦路径：预算耗尽必须
UNRESOLVED，足够查询应 FOUND；有限面 clear-path 证书包含 TOP。

沿用旧组件容限：解析力/能量标量绝对误差 1e-12、rtol=0；几何量 1e-12。
scalar/batch 数值数组逐字段 rtol=3e-13、atol=3e-13，布尔、标识和非数组元数据一致。
一般坐标旋转力向量最多沿用旧同类组件的 2e-12 绝对误差，功仍 1e-12。
内部有符号法向/切向功误差必须满足求解器原有逐项 roundoff certificate；不修改
证书、负摩擦判断、吸引法力拒绝或任何过程门槛。尺度夹具用解析同次缩放参考。

runner 保存实际命令、退出码、测试名/结果、源码/合同 SHA、原始 stdout/stderr；
已有目录拒绝覆盖。成功标记仅为 PASS_TOP_CONTACT_COMPONENTS_NOT_PROCESS_VALIDATION。
真实 step_022 中 shell / visible-path 尚未执行，不能引用本组件 receipt 替代它们。
