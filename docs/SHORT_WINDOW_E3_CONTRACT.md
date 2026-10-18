# E3-3：一次0.5us集成窗口，运行前冻结（2026-10-01）

用户要求继续上一轮后的预算内集成检查。本轮不改src/、旧配置或接触算法。
工程近似目标及旧FAIL不变；这不是原2us目标完成、收敛或沉积验证。
只运行一次；失败或超时后记录，不在本轮继续修接触/换参数/重新发起。

## 固定问题、窗口与成本

- 原configs/debug/starved.json、dx=.5mm、CFL=.1、657粒子、层厚1.5mm；
  材料/热学/接触刚度/执行器厚度/温度/速度/几何/约束全部保持。
- PySPH 1cd3dea9cfaac7d9f1bda2a2491ce4b42a31dbe5；
  event_aware_v1含上一轮已验证的contact_cone逐盘数值缩放。
- 从原初态重新开始，不迁移旧代码checkpoint，不恢复任何REJECTED。
- 本轮明确的短窗口：5e-7s=.5us；最多80步，每个接受步输出。
  选取理由：覆盖此前候选22及历史27附近，并给后续一段推进空间；
  不是为了截取已知成功终时，结果尚未知。不是2us替代验收或dt敏感性。
- 上轮测得dt=1.0490609943676047e-8s；预计ceil(5e-7/dt)=48步。
  首步3.5764233s给出171.6683184s线性估计，乘1.3约223.17s。
  内部硬上限300s，外部整个进程timeout330s，低于600s。
  不以较短窗口结果回写原2us预估，也不在结束后自动延长。

## 前置证据及接受条件

1. 先验证CONTACT476_SCALE_FREEZE.json，旧830文件、上游pin/clean。
2. 源/配置未改变，复用已冻结的116项111PASS/5旧FAIL以及432组证据；
   本轮只快速重跑新64组/独立参考/完整476kick和原476测试，不重复整套。
3. 质量相对1e-12、密度5%、穿透/dx10%、速度20m/s、
   能量残差5%交换尺度+1e-7J、材料温域及全部有限值原保护不放宽。
4. 唯有runner COMPLETE、accepted series/checkpoint实际到5e-7s、
   所有接受诊断无报警、hash/ID/标签/字段一致才称SHORT_WINDOW_COMPLETE。
   即使成功，paper/deposition/mixing/production仍false。
5. 与原第22步门槛比较仅报告是否跨过“步号/物理时间”，
   不声称重新经历完全相同的粒子状态或相同事件顺序。
6. 逐帧读回真实坐标/速度/温度/密度/应力/ep、ID/被动标签与工具位姿；
   已知Es不是严格自由能，只把账面残差标为当前近似模型诊断。
7. 若失败：checkpoint仍是最后有效状态；候选REJECTED不可恢复。
   辅助last_completed_contact_callback_trace可能陈旧，只用failed_local_context定位。
   HARD_TIMEOUT时manifest可能仍INCOMPLETE、performance缺失；不得伪造补齐。

## 计划命令与本轮停止

    .envs/pysph/python.exe -B scripts/run_afsd_debug.py --case configs/debug/starved.json --output outputs/short_window_e3_20261001/dx05_cfl01_05us --contact-variant event_aware_v1 --end-time 5e-7 --steps 80 --output-every 1 --wall-seconds 300

外层subprocess.run(timeout=330)，线程环境OPENBLAS/MKL/OMP均1。
详细shell命令见logs/short_window_e3_20261001/COMMANDS.md。
跑后只做读回/验证/交接；不启动第二算例或GPU/并行/生产。
