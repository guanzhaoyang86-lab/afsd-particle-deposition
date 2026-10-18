# 已建立 TOP 接触后的独立时间步诊断合同

本合同及新 controller 在任何本合同真实候选运行之前登记；不修改冻结物理、旧
controller、旧合同或旧 FAIL。方法仍为
`HOT_QS_VISIBLE_SUPPORT_CHAIN_RULE_CONTACT_APPROXIMATE_V1`，PySPH commit、参数、
空间算子、摩擦/材料/热学和现有门槛不变，不做质量、时间、模量缩放。

## 要回答的问题与范围

首次 TOP 创建区间 42.875–44.875 ms 的 250 us 诊断已经在首个共享终点出现局部
应力差约 17.74% / 11.54%，超过既定 10% 门槛；原
`outputs/visible_top_timestep_20261004/mult_04/manifest.json` 保留 FAIL。
更大的首次事件 500/1000 us 分支未继续，不能视为通过。

本轮另问：TOP 接触已经存在时，250 us 是否仍与 62.5 us 的实际接受轨迹一致？
新入口 `screen_visible_established_contact_timestep.py --window N --multiplier M`
仅允许 N=2..100、M=4/8/16。第一项实际拟运行是 window 2、M=4：从 **window 2
的真实 44.875 ms 父态**重放 2 ms，与 window 2 的 32 个 62.5 us 接受步比较。
绝不从 46.875 ms 终态开始，也不沿用首次事件失败候选。

M=8/16 对应 500/1000 us，只保留显式独立入口，本轮初始状态均为 **NOT_RUN**。
入口一次只执行一个诊断，不自动串跑、不自动升级时间步。后续若运行必须分别报告
实际命令、预算及结果；250 us 失败后，不自动执行更大时间步。

本诊断成功只能限定为当前已建立接触窗口与一种细步参考一致，不能洗白首次事件
FAIL，不能据此全程采用大步，也不能证明下一次接触创建、释放、支撑边切换或入孔
事件可靠。事件发生时可能仍需要小步；接触 history 的时间收敛仍记未建立。

## 前置证据及保持条件

- 该 N 号窗口自身的独立 readback 必须先完整 PASS，32 个实际细步、完整 TOP
  事件资格及连续能量前缀均成立；部分接受前缀不得供本诊断使用。
- 使用 `top.preflight(N, 'coarse')` 获取同一真实父态和冻结来源；不 monkeypatch
  旧脚本的全局变量或函数。
- 从原始已审计父态接触 journal 读取 TOP active，父态必须已有实际 TOP 接触；
  journal 与其已通过 readback 的 SHA 绑定。不为该判定重算几何或力。
- 原 15 项 TOP 组件收据继续由 top preflight 核验，不重复声称新组件或物理认证。
- 所有参考 accepted checkpoint 的 SHA、父链、dt、完整父能量前缀以及独立
  readback 中对应 checkpoint 路径/SHA 必须一致。所有候选/参考来源运行期间
  不变；旧首次事件失败 manifest 也作为只读引用绑定。

## 每个候选及每个共享终点的接受条件

候选仍调用冻结 `visible_top_candidate.candidate`。保留 fresh 完整残差、材料/热学、
质量/密度/温度/穿透与能量门槛、scalar/batch 完整候选比较、完整能量前缀、有限
接触壳证书以及保留支撑边整步路径证书。非孔面只允许 SHOULDER 与已资格 TOP。

每个大步终点均比较实际 62.5 us 参考的相同时间状态，不只检查整个窗口末态：

- 全局与局部位置差 ≤5 um、温度差 ≤0.5 K；应力相对差 ≤10%，塑性相对差 ≤20%。
  局部应力/塑性分母分别保留 1 MPa / 1e-5 下限。调用参考 fine 状态作为
  `post.compare` 和 `event.local_compare` 的第一个参数，不倒置归一化。
- 累计种子包括候选及覆盖区间内**每个**细参考子步的接触事件、blocked IDs、
  crossing IDs；每个子步的原始和过滤图均加入一跳邻域，不能只取区间末图。
- 每个共享终点的严格 crossing ID 列表必须一致；保留逐粒局部表、双方列表、
  参考细步编号、来源 SHA 与候选全部原始诊断。
- TOP 实际活动必须出现在整个完成的候选分支内。任何失败立即停止，保存 FAIL，
  不放宽门槛、不改分类、不从失败候选继续另一个分支。

全部候选 checkpoint 保持 `restart_allowed=False`，不写 accepted.npz；manifest
的接受步数和接受时间增量永远为零。只有诊断成功状态
`PASS_ESTABLISHED_CONTACT_TIMESTEP_REFERENCE_SCREEN_NOT_RESTARTABLE`，没有过程
接受、全局增 dt、完整沉积或实验验证授权。

## 预算与运行

参考成本为该 N 号完整 fine stage 的 wall_seconds / 32，同时取与其中最大实际
候选成本的较大者。单个大步预测 = 上述成本 × M；候选数 = 32/M；完整预测为
`30 + 1.5 * 候选数 * 单个大步预测` 秒，必须 ≤450 s。

单次内部 540 s、外部 600 s，求根留在 525 s 预算内；OMP/OPENBLAS/MKL 各 1
线程，无 GPU。每个窗口/倍数使用独立输出和监督日志，拒绝覆盖旧文件。

初始运行命令（仅在完整 window 2 readback PASS 后）：

```powershell
.\.envs\pysph\python.exe -B scripts/screen_visible_established_contact_timestep.py --window 2 --multiplier 4
```

controller 创建阶段只作 AST/语法和静态结构核对，不调用 preflight 或物理求解。
源文件和合同一旦用于运行将由其证据冻结，不再就地修改。
