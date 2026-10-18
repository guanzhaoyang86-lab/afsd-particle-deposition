# 有限开放方孔工具组件（2026-10-03）

## 运行前冻结范围与验收

分类：`FINITE_OPEN_SQUARE_DISTANCE_APPROXIMATE_COMPONENT_NOT_AFSD`。
这是独立几何/保守法向势能组件，不接现有求解器，不删除或放宽任何旧位移、
feature transition 或穿透 guard；没有摩擦、热、材料或自由粒子推进。
所有尺寸用 SI。制造测试尺寸不是论文参数。

工具实体采用局部集合 `0 <= z <= H, x*x+y*y <= R*R,
max(abs(x),abs(y)) >= a`，其中 `2a` 为贯穿开放方孔宽，`R > sqrt(2)*a`。
枚举下平肩、顶面、有限外圆柱、四个有限孔壁的精确最近点；平面环域投影
枚举孔四边。取边界最小欧氏距离，实体内部为负、自由空间为正。
这不是用 `max/min` 布尔 SDF 代替真实欧氏棱边距离，不圆滑或封闭方孔。
外侧有唯一最近点时法向为该距离的梯度；同一棱边最近点可同时拥有两个
feature 名，feature 只是诊断身份。距离为零的尖棱以及不同最近点等距的
中轴有 `NONSMOOTH` 状态；活跃法向力明确拒绝，不平均法向、不 pinv。

新的**单实体距离惩罚变体**为 `delta=max(radius-d,0)`、
`U=0.5*k*delta^2`、`F=k*delta*grad(d)`。正部仅定义单边几何重叠，
不是裁功/热/温度。势能由单实体距离定义，不随有限 face 名字切换重新归零。
此变体不等价于旧逐面叠加弹簧；非等距的多个孔壁同时穿透时只表示最近实体
边界，不能据此声称旧多面摩擦/执行器已获验证。集成必须单独声明接触模型。

运行前容限：独立解析距离/法向/力/势能绝对 `1e-12`（制造尺寸 0.1–2 m）；
毫米尺寸几何绝对 `1e-15 m`。中心差分 `h=1e-7 m`，力/势能梯度容限
`2e-7 N`、刚体时间导数功容限 `2e-7 W`。这些容限只检验双精度几何公式，
不是 AFSD 离散误差。连续棱边转换前后距离/势能极限按解析 `epsilon=1e-8`
检查，独立路径功积分误差 `2e-6 J`；无在测试后调整容限的授权。
旋转协变、刚体点速度、开放孔道、非光滑拒绝、非法输入亦须测试。
单次测试进程外部 timeout 60 s；最多两轮证据修复，首错日志保留。

## 接口与集成风险

- `FiniteOpenSquareTool(radius, bore_width, height)`：只接受有限正尺寸。
- `query(point, origin, rotation, translation_velocity, angular_velocity)`：
  返回 signed distance、最近点集、feature 集、唯一法向或缺失状态，以及同位置
  刚体速度。一般正交 3D 姿态可用；转动/平移由调用者提供，不驱动材料。
- `normal_penalty(query_result, particle_radius, stiffness)`：独立单实体法向
  力/势能；非光滑且活跃时抛出带查询证据的异常。

不接受 tracer/AS/RS。没有基座、摩擦历史、温度、接触事件时间积分、执行器、
插入粒子或 restart 适配。一般最近点速度与粒子中心速度位置必须功共轭；
法向力沿两点连线所以两位置的法向刚体功相等，但摩擦不能这样互换。
尖棱切换需事件定位/多约束可行性；不得接上组件后直接取消旧 guard。
活跃非光滑点拒绝仍是完整通过孔口的集成阻塞；组件 PASS 不能当作有限棱边
时间推进已解决。真实 AFSD 集成、输运、旧多面摩擦回归和三工况均为 `NOT_RUN`。

## 来源

集合和欧氏投影公式为本项目几何推导；原项目有限平肩/方孔定义来自
`geometry.py`。旧 `IntegratedSolver._feature_transition_guard` 明确拒绝带势能
face 消失，说明该 guard 不是过时错误。刚体 `v=V+omega cross (x-origin)`
与现有 `rigid_contact_kinematics.py` 一致。本组件不复制第三方代码、不改上游。

## 实际结果

首轮 `finite_tool_tests_initial.log`：21 tests，19 PASS、2 ERROR，外 0.400375 s。
一项为全局距离容限将距棱边 `1e-8 m` 的平滑面误判为两个最近点；平方距离差
应以 `(qa-qb) dot ((qa-x)+(qb-x))` 稳定比较，舍入界 `64 eps sum(abs(terms))`。
另一项为当前 NumPy 无 `np.trapz`，测试改为显式梯形和。首错保留。
修复轮 1：上述比较修正并新增回归；不改原验收容限或物理。
`finite_tool_tests_repair1.log`：22 tests 全部 PASS，内 0.457 s、外
0.610479 s、exit 0；尚余一轮允许修复但无证据需要使用。

实际命令为 `.envs/pysph/python.exe -B tests/review/test_finite_tool_geometry.py`，
由独立 Python `subprocess.run(..., timeout=60)` 包装，stdout/stderr 合并保存。
日志目录：`logs/process_completion_20261003/`。未启动真实 AFSD 或长任务。

独立基准包括：3-4-5 孔口/肩外缘距离 `0.05 m`、法向分量 `0.6/0.8`；
`k=100 N/m, radius=0.1 m` 时 `U=0.125 J`、力模 `5 N`；
尖棱半径变化路径 `r:0.04 -> 0.06 m` 的法向功 `0.1 J`（端点解析势差）；
刚体平移/旋转势能时间导数与壁面法向功一致；毫米缩放、一般三维姿态、
平肩/孔壁/棱边连续最近点转换、非光滑显式拒绝。

结论：**有限几何与独立单实体保守法向分支 PASS；真实接触时间积分/输运
NOT_RUN，完整孔口通过仍 BLOCKED**。已有 guard 不动，不能据此授权生产。
