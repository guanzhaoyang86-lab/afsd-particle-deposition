# 独立两次热库存出生诊断合同 V1

本合同只授权一个从 8 mm 独立 initial_v2、t=0 开始的两步诊断，不延长任何现有接受分支。
模型分类为 `HOT_QS_TWO_INVENTORY_BIRTH_DIAGNOSTIC_V1`。任何阶段和 checkpoint 都必须
`accepted=false`、`restart_allowed=false`；V2 scope 固定为
`DIAGNOSTIC_COUPLED_INLET_BENCHMARK_NOT_PROCESS_ACCEPTANCE`，不能升级过程接受。
与预装长棒的独立对照、时间收敛、持续质量流率入口及完整沉积均为 NOT_RUN/NOT_MODELED。

## 初始条件与接口单位

唯一原点是 `outputs/expanded_hot_domain_20261004/initial_8mm_visible_v2/INDEPENDENT_INITIAL_NOT_RESTARTABLE.npz`，
SHA256 `31a8ae1910f8189ab3f437f26f663b523fbf0d8195bd41e1a78cc000696fcfa3`。
保留其实际初始分类、零时刻来源和完整数组；新分类明确记录为本独立诊断的模型范围，
不伪造旧 spatial transition 或过程父 checkpoint。物理参数使用 expanded_hot_domain.parameters(.008)，
feed/traverse/omega/rho/G 不缩放；粒子间距 0.5 mm，所有计算采用 SI。

在 t=0 和 t=3.90625 us 各添加 25 个 5×5 方层，中心 z 分别为 3.25/3.75 mm，
x/y 为 −1,−0.5,0,0.5,1 mm。出生温度 723.15 K，rho=rho0，S=0，ep=0，
塑性功/回归耗散/接触史/弹性储能均为零，速度为 (0,0,−feed) m/s。
速度是显式 reservoir 初值；QS 保存的后续速度为 dx/dt，出生动量和动能仅观察、不得转成热。
这是外部库存出生，不用 25m/dt 冒充规定物理质量流率。

入口控制区域为 component=1 且 world-z>=2.5 mm，每个机械步准备时重选，仅 z 规定 −feed*dt；
x/y 始终自由，原基板底层三轴固定。出生后 N=943/968、free DOF=2011/2061、受控点50/75。
激活/释放 IDs 保留；两步没有释放时结果必须写 NOT_EXERCISED，不能说释放已验证。
标签不参与力热。内部轨迹由实际求根产生，不做预制运动。

`scripts/two_insertion_diagnostic.py::InletTrial` 直接初始化现有 visible 完整 map 所需属性，
复用其 evaluate/tangent，无旧分类伪装。固定库存的**一步内部**用 V2 observer 的 legacy view
调用原 discrete_process.accepted_copy，只获得 detached 候选及原 fresh/raw/all-prefix 检查。
V2 seal 根据真实出生 registry 和坐标重新计算 observer；最终落盘永远保存 V2 birth observer。
第一步封存后读取真实 NPZ，反转 PA 行序，再按 persistent ID 恢复 contact/free；第二步必须
从该读回状态准备。这个诊断读回不能调用过程 restore 或声称 checkpoint 可以续过程。

## 每一步必须具备的独立证据

两个根均重新创建 trial，出生之间不复用 secant。完整 map 的实际可见 i/j/g 同时用于
密度/梯度/成对力/热及切线。原 force/heat/penetration/negative work/raw energy/全前缀门槛不变。
每根再做 fresh 完整 map、scalar/batched 完整场比较、实际图文件、修正矩阵条件、连续 shell
和 retained-edge path（32层/256查询）、原接触扫掠与 SHOULDER/TOP 事件限制。
每个实际候选图另组当前几何的**无载弹性参考**自由刚度，核 eigen/rank/least-mode 能量恒等式；
它只排除图/约束的弹性零模态，不证明实际非线性本构切线正定或未来流动稳定。
新生粒子与既有/同层中心不得重合或 nominal cells 初始交叠，出生接触势能必须严格为零。

source-work 使用已冻结 inlet_source_work 的真实 parent/prepared/actual、实际反力位移、
接触两端储能和力热账。外部出生边界明确牵引为零，才允许 pressure-flow=0；内部控制面
不另加压力流功，以免与实际 actuator 反力功重复。统一 Tref=parameters.room。
逐步开放库存余额保存并累计，残差保留而不补热；source-work 不替代任何求根或热门槛。
V2 外部 verifier 必须读回实际八 gate artifact、实际 source artifact、绑定 PA/history/time
及文件 SHA，不能用无条件 True。显著原始能量缺陷仍由原 gates 拒绝。

## 成本、输出和失败行为

CLI：`python scripts/run_two_insertion_diagnostic.py --stage cost --output <全新项目内目录>`，
之后同目录 `--stage diagnostic`。`--preflight-only` 只读核参数、输入、来源及成本资格。
所有输出文件新建，旧 FAIL 不改。每个 worker 540 s硬停、外部600 s；两根加验证的预测必须
<=450 s，成本不足直接拒绝。cost 只在第一出生后的完整 map 测量 full-map、实际 full tangent+
linear 及实际图参考 rank 成本；第二库存成本按自由 DOF 三次比显式估计，并乘1.5安全因子。
求根最多24更新/100 map/两次 Jacobian刷新，实际更新上限由保守成本公式收紧，无额外 fallback。
预估公式及各项测量写入 receipt；成本预测不等同实际壁时保证，任何后续截止仍 fail closed。

只有两次真实根及所有 gates/V2/source/读回链均成功，才标
`PASS_TWO_INSERTION_DIAGNOSTIC_NOT_JOINT_OR_PROCESS_ACCEPTANCE`。
即使成功，accepted_steps=0、accepted_time_advanced_s=0、joint_benchmark=NOT_RUN；
诊断时间另列 7.8125 us，最终库存968。失败保留已完成诊断前缀和最后未接受候选，不重置历史。

实现准备和语法检查不是执行证据。主代理排期前不运行 cost、NNPS、真实根或物理测试。
