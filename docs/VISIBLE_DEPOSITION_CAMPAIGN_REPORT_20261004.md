# 真实输运推进至 49.9375 ms：过程、失败及当前边界

当前正式结果为 **49.9375 ms / 655 个接受步**。这是 342 粒、独立初始化为
450°C 的工程近似热态准静态分支；不能与原惯性分支时间相加，不能称完整沉积或论文
定量复现。从本次连续推进起点 40.375 ms / 447 步起，共增加 **9.5625 ms / 208 步**。

实际严格穿越粒子由 238/262 增至 **226/238/262/274/298**，总质量 **1.6875 mg**。
末态五者均仍在方孔投影内。穿孔后继续塑性变形不等于已形成肩外铺展。

最终被动观察核对了 **656 个真实接受帧**，跨孔后肩外质量、半粒距分离质量及
实际肩面接触质量均为 **0**，没有将初始化就在肩外的材料算作新沉积。
按 PID226/238/262/274/298 顺序，跨孔首个真实样本后的连续采样塑性增长时长为
0.625 / 12.125 / 12.625 / 5.9375 / 1.625 ms；这不是对子步连续率的证明。
观察与真实动态图位于 `outputs/deposition_observer_20261004/readback_49_9375ms_v1/`，
其中 summary SHA256 为 `b881396672665e1d6932d27eead695d46a8554e798680ce81f7f852144415d98`。
`deposition_progress.png` 已实际渲染、查看，明确标注 NOT FULL DEPOSITION。

机器可读的整轮索引为 `outputs/visible_deposition_campaign_20261004/report_v1/summary.json`，
SHA256 `42e4421516b3e582f6185b913c1c17731d14ebb9e646a13fb73a15e45ca9dfd9`；它逐窗绑定
正式父态、独立摘要和真实监督收据，物理求值次数为0。

## 改动与验证

1. 新增有来源与预算限制的续跑 controller，保留真实父态、全部能量历史、材料与
   接触 history、被动身份和原空间模型转接。粗支永不成为 restart；细支仅在一整对
   全部物理/路径门槛和同父全局/局部比较通过后发布。失败前缀单独记账。
2. 通过 15 项独立 TOP 组件检查后，新增明确的 SHOULDER+TOP 事件资格。原残差已
   包含这些有限面；实际力、热、材料、摩擦及势能公式没有更换。新的真实 TOP 窗口
   两支均出现 198/222/318/342，并完成逐步 shell 与支持边整步证书及粗细比较。
3. 实测大步失败后保持门槛不变；另从真实 46.500 ms 父态采用 62.5/31.25 us
   局部细化。没有把 250 us 失败候选、旧 TOP 失败候选或未通过比较的两半步发布。
4. 新增独立 saved-evidence reader 与被动沉积阶段观察入口，保留所有真实帧；
   无 NNPS/force 的读回只核对证据与算术，不能替代物理验证。

TOP 组件实际命令为 `.envs/pysph/python.exe -B scripts/run_visible_top_contact_gate.py --name v1`，
exit 0、15 tests。收据 `outputs/visible_top_components_20261004/v1/summary.json`。
本轮没有全量回归；历史 1246 项的 10 FAIL / 1 ERROR / 3 expected failures 未重跑。

## 实际物理运行

下表耗时为各 stage 的真实外部监督时间。失败状态与合格前缀同时保留，reader 的
PASS 仅代表保存证据一致性，不能把原物理 stage FAIL 改成 PASS。

| 系列 / 窗口 | 新接受步 | 正式终时 ms | 粗支 exit / s | 细支 exit / s | 窗口结果 |
|---|---:|---:|---:|---:|---|
| visible_deposition / 001 | 40 | 42.875 | 2 / 65.333 | 2 / 155.008 | TOP 原白名单阻断，保留前缀 |
| visible_top_event / 001 | 32 | 44.875 | 0 / 53.144 | 0 / 126.494 | 完整 2 ms 比较通过 |
| visible_top_event / 002 | 26 | 46.500 | 0 / 54.002 | 2 / 105.680 | PID270 局部应力差失败 |
| visible_local_refinement / 001 | 32 | 47.500 | 0 / 51.419 | 0 / 146.225 | 完整 1 ms 比较通过 |
| visible_local_refinement / 002 | 32 | 48.500 | 0 / 49.658 | 0 / 121.765 | 完整 1 ms 比较通过 |
| visible_local_refinement / 003 | 32 | 49.500 | 0 / 52.166 | 0 / 130.549 | 完整 1 ms 比较通过 |
| visible_local_refinement / 004 | 14 | 49.9375 | 2 / 29.027 | 0 / 59.081 | 仅前 7 对合格，完整窗未完成 |

14 个 stage 合计约 1,199.55 s；另一次 250 us 非接受诊断耗时 9.702 s。
每个真实 stage 均有 <=450 s 事前预测、540 s 内部 / 600 s 外部限时，实际最大
155.008 s。保持单线程 CPU、原质量/模量/真实物理时间、300 rpm 与原走速/进给；
没有 GPU、三工况生产或缩放。这些实耗不能线性外推为完整沉积完成承诺。

执行入口：

```powershell
.envs/pysph/python.exe -B scripts/run_visible_deposition.py --window 1 --stage coarse
.envs/pysph/python.exe -B scripts/run_visible_deposition.py --window 1 --stage fine
.envs/pysph/python.exe -B scripts/run_visible_top_event.py --window 1 --stage coarse
.envs/pysph/python.exe -B scripts/run_visible_top_event.py --window 1 --stage fine
.envs/pysph/python.exe -B scripts/screen_visible_top_timestep.py --multiplier 4
.envs/pysph/python.exe -B scripts/run_visible_top_event.py --window 2 --stage coarse
.envs/pysph/python.exe -B scripts/run_visible_top_event.py --window 2 --stage fine
.envs/pysph/python.exe -B scripts/run_visible_local_refinement.py --window N --stage coarse
.envs/pysph/python.exe -B scripts/run_visible_local_refinement.py --window N --stage fine
```

最后两条实际 N=1,2,3,4；每窗后运行对应独立 reader。所有入口拒绝覆盖旧目录。
具体命令、退出码、stdout 与监督记录在各系列 `logs/`，独立 reader 在对应
`outputs/<series>/window_NNN/readback/` 保存 summary、trajectory、完整能量账与执行收据。
**不能自动运行 refined window5**：window4 未满足完整窗口资格。

## 不能被改判的失败

- 250 us 相对真实 62.5 us 的首次 TOP 创建试验：首共享终点局部应力差
  17.7404% / 11.5400%，超过 10%；0 接受步，500/1000 us 未继续。
- TOP window2：PID270 局部应力差逐步积累至 10.0721%。其自身无接触，参考应力
  范数超过 1 MPa，下限未触发。后来的 1 ms 局部细化通过不证明旧全历史时间收敛。
- refined window4 coarse step008 的 50.000 ms 候选：支持边整步证书用完
  256 次查询，仍有 PID339–341 在 fraction [0.875,1] 未认证。实际最深细分8，
  上限32；没有 found_edges、退化边或数值域失败。不能称已发现穿透或释放边。

随后独立纯几何诊断分别对两个剩余子区间使用原 API，64 / 96 次查询得到证书，
未发现交集。父、端态都保留该边，不支持 start-blocked/end-visible 的解释。
这说明原全局包络查询预算是明确阻断；诊断没有追认原候选或新增接受时间。
结果 `outputs/visible_path_diagnostic_20261004/quarter_turn_v1/summary.json`。
仍须另行登记、更充分验证几何包络效率或总预算方案，再从真实父态重新求解与比较。

一个尚未实现的改进是避免现界 `||v|| + theta*max||r||` 丢失共转抵消。
对于线性世界弦和精确 SO(3)，若 `q(s)=r(s)R(s)`，可由
`||q'(1/2)|| + (2*theta*||v|| + theta^2*max||r||)/2` 给出全步速度上界。
实际浮点实现还必须处理 R0 正交缺陷，以及 q' 抵消时按两个大项各自尺度计入的
绝对舍入误差；不能只对小的差值加相对误差。小角/近 pi 的回退、边界夹具和真实
保存路径仍需独立资格。本轮仅审查数学可行性，没有实现、运行或宣称新界通过。

细节见 `VISIBLE_TIMESTEP_FAILURE_ANALYSIS_20261004.md`。后事件大步入口仅完成静态
实现，其所需完整参考窗口没有建立，250/500/1000 us 均为 NOT_RUN。

## 正式终态与账本

- 末态：`outputs/visible_local_refinement_20261004/window_004/fine/step_014/accepted.npz`
- SHA256：`396080e1fdd6d9f211fb521d2c44e47ecb19c17c0d327d781b3bb5f61156d4f0`
- 独立摘要：同 window004 下 `readback/summary.json`
- 摘要 SHA256：`6659f3459c989641885d5786df76d26600abae1f9b5e7b67d23ba43b118d7ce1`
- 总质量 0.000115425 kg；末温 723.339255–723.985764 K；最大塑性应变
  0.4845456662；最大等效应力 11.8945543 MPa。
- 全 655 步 raw residual 有符号累计 6.01298648 uJ；逐步绝对值累计
  **6.02386653 uJ**；absolute work scale 0.08529587245 J。
- 最后合格 0.4375 ms：raw residual 绝对累计 2.282565e-10 J，外功
  4.07388865e-4 J；最大局部应力差 1.61299%。

这些是明确工程近似下的内部账与时间比较，不是实验误差。旧 5% work-scale
能量门槛没有放宽，也不能当作物理精度保证；未裁温度、负耗散或材料响应。

## 完整沉积仍缺什么

当前 4 mm 基板从初始就小于 5.4 mm 肩径，只有有限料棒库存，没有持续入口事务。
继续累加当前小域时间不能解决完整肩覆盖或补料。现成入口组件仍缺接受后 registry
同步、开口力学功/对流储能、执行器层释放和可变粒子身份/观察/重启接入。

具体可复用函数、6/8 mm 独立域与 dense 成本、最小入口联合验证见
`FULL_DEPOSITION_INTEGRATION_NEXT_20261004.md`。扩大域必须独立初始化，不能在当前
655 步 checkpoint 上补粒子后冒充同一封闭过程。FKS 低率与参数、Jaumann 闭循环、
空间/边界敏感性及实验数据问题仍未解决；当前项目未完成。
