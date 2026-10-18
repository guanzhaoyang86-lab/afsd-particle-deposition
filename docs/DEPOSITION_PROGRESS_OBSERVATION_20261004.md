# 44.875 ms 实际接受链的沉积阶段观察

此后处理没有新增接受时间，没有调用 PySPH、NNPS、接触、材料或热求解。最终采用
519 步、520 个真实接受时刻；原动态分支与此独立 450°C 热态准静态分支不相加。

`readback_deposition_progress.py` 先读取原 432 帧和 VISIBLE_SUPPORT 新增 16 帧，精确
重放已签名的 448 帧观察结果。随后依次追加 40 个及 32 个真实帧；只删除精确相同的
共享父帧。每帧绑定实际 accepted checkpoint SHA、ID、质量、位置与塑性应变；新增
接触活动来自 PASS readback 已绑定的原始完整接触 journal，并与 checkpoint 接触
history 逐项相等。初始接触未计算仍标不可用。所有输入在运行前后哈希一致。

## 实际观察

| 量 | 42.875 ms | 44.875 ms |
|---|---:|---:|
| 实际接受帧数 | 488 | 520 |
| eligible 中已严格穿过方孔出口的质量 | 0.675 mg | 1.0125 mg |
| 新穿越材料中心在孔投影外、肩面下的质量 | 0 | 0 |
| 新穿越材料半个 dx 横向分离的质量 | 0 | 0 |
| 新穿越材料的实际肩面接触质量 | 0 | 0 |

PID274 新增穿越事件的记录时间为 43.9926225133 ms，所在真实接受区间为
43.9375–44.0000 ms。44.000 ms 是首个真实穿越后状态，不为穿越瞬间插值塑性应变。
最终 PID238/262/274 分别在出口面下 29.356/63.386/6.711 µm，但仍距方孔投影边缘
790.503/759.640/801.293 µm；三者均尚未横向铺展到肩下环形区域。

从首个真实穿越后样本开始，三者累计塑性增量分别为
0.03906197345 / 0.03397830641 / 0.004214229664。
实际相邻样本的正塑性增量分别覆盖 113/121/14 个连续区间，即
7.0625/7.5625/0.875 ms；未发现超过 roundoff 容限的负增量。
这是接受样本之间的增量观察，不证明任意子步时刻的塑性率连续。

44.875 ms 时初始就位于肩面下、当前孔外肩下的棒料有 12 粒，实际肩接触有 2 粒。
它们不具备最初上方 eligible 资格，不能计入新送料沉积。

## 342 粒有限域可以回答多长过程的问题

当前域可以研究经过各步门槛审核的有限库存、局部出口及随后铺展事件。它不能给出
完整肩部稳定沉积的有效时间窗：基板名义宽 4 mm，而肩径 5.4 mm，初始时就未覆盖
整个肩部投影。增加时间步数不能消除这一几何限制。

初始 eligible 库存只有 75 粒、25.3125 mg，等效于名义棒截面的 1.5 mm 长度；无持续
补料。按规定进给 2.12 mm/s，一个 dx=0.5 mm 对应 235.849 ms；1.5 mm 对应
707.547 ms；300 rpm 的一转为 200 ms。这些仅是执行器运动和长度换算，不是自由
材料通过时间、可验证过程持续时间或计算耗时预测。当前规定进给和工具走过各
95.135 µm；最近尚未穿越的 PID298 中心距出口 17.857 µm，也不能据此承诺下一步必过。

要称完整持续沉积，需要另行处理足够宽的基板、入口库存/补料及相应边界验证；若改变
初始化，应另立分支范围，不能把新状态冒充当前 519 步的自然连续结果。

## 重跑方式与结果

```powershell
.\.envs\pysph\python.exe -B scripts/readback_deposition_progress.py --readback outputs/visible_deposition_20261004/window_001/readback --readback outputs/visible_top_event_20261004/window_001/readback --output outputs/deposition_observer_20261004/NEW_UNIQUE_DIRECTORY
D:/py/anaocnda/python.exe -B scripts/plot_deposition_progress.py --observation outputs/deposition_observer_20261004/readback_44_875ms_v1
```

观察 driver 可继续添加 `--readback`，必须按真实父态顺序排列，并使用未存在的输出
目录。现有图及收据已存在，绘图入口将拒绝覆盖；当前绘图文字专门描述此三粒仍在孔内
的状态，不能直接拿来描述未来已有肩外输运的状态。

- 42.875 ms 观察：exit 0，3.853 s。
- 44.875 ms 观察：exit 0，5.269 s。
- PNG 绘图：exit 0；已通过 `view_image` 人工视检，标注 `NOT FULL DEPOSITION`。
- 观察 summary：`outputs/deposition_observer_20261004/readback_44_875ms_v1/summary.json`，
  SHA `b9f8be71115bf87a9463f9cd21cff0647bd77a2ca7a844effc21005ca384b5be`。
- 图：同目录 `deposition_progress.png`；真实输入、逐帧来源及采样持续证据分别为
  `observation_inputs.npz`、`actual_frame_bindings.json`、`plot_receipt.json`。

没有改动冻结物理源码、旧观察模块、旧 checkpoint、旧 FAIL 或已有报告；上述 PASS
仅指只读观察和可视化，不是全量测试、模型物理验证或完整沉积达成。
