# 两次局部时间步失败及后事件大步诊断状态

本报告读取已保存的真实候选、接受 checkpoint、原始接触 journal、时间比较和监督
日志。没有 force/NNPS/材料/热求值，没有新增接受时间。未采用任何不完整的 reader
v1 JSON；下列结论直接绑定原运行证据。

## 当前结果

1. 首次 TOP 创建区间的 250 us 诊断从真实 42.875 ms 父态出发。在首个共享终点，
   全局应力差 0.7044%，但 PID318/324 局部应力差为 **17.7404% / 11.5400%**，超过
   10% 门槛。保留 `FAILED_NONRESTARTABLE_TIMESTEP_SCREEN`；exit 2、外部 9.702 s，
   0 接受步。首次事件的 500/1000 us 未继续。
2. 随后的 TOP window 2 使用既有 125/62.5 us 比较，细支实际接受 26 步至
   **46.500 ms / 全局 545 步**。pair14 的 46.625 ms 比较失败，fine27/28 两候选均
   未发布。原细支状态为 `FAILED_WITH_ACCEPTED_PREFIX`，exit 2、外部 105.680 s。
   此次全局应力差 0.3530%，PID270 局部应力差 **10.0721%**，严格 crossing IDs 一致。

## PID270 的已存证据

PID270 为杆顶执行器层点，actuator=1、free=[True,True,False]。fine26/27/28 和
coarse14 中，它的全部接触 active/sliding、接触力、摩擦耗散、释放能和算法耗散均为
零。没有该点自身 TOP 接触创建/释放证据。

其局部应力差从 pair1 的 0.4559% 逐步增至 pair13 的 9.4183%，再到 pair14 的
10.0721%。末对位置差 17.076 nm，温差 1.49318e-5 K；细/粗张量范数为
1.228152/1.228165 MPa，范数相对差仅 1.07318e-5，而张量差范数为 0.123701 MPa。
差别主要落在张量分量/方向；展开为九维分量向量的夹角是 5.7733°，该量不是材料
刚体转角。细/粗 strain_rate 为 0.0287283/0.0369761 s^-1，ep_rate 为
0.0284486/0.0368427 s^-1。

末端两支图相同，PID270 邻居相同；pair11 的全图曾暂不相同。观察支持进一步检查
累积时间离散和材料/空间耦合响应，尚不能从这些已存量确定因果来源。没有依据直接
放宽 10% 门槛，也不能把本次失败解释成重复的局部 TOP 创建事件。

## 后事件大步入口：本轮 NOT_RUN

已新增 `screen_visible_established_contact_timestep.py` 和
`VISIBLE_ESTABLISHED_CONTACT_TIMESTEP_CONTRACT_20261004.md`；AST compile 和静态
结构核对 exit 0，没有调用 preflight、没有真实求解、没有重复组件测试。

该入口要求所选窗口自身具有完整 32 个细参考步及完整独立 readback PASS，并且
父态确实已有 TOP 接触。**window 2 只有 26 个接受步，因此后事件的
250/500/1000 us 全部 NOT_RUN。** 不裁成部分参考以满足该入口，不从未接受的
46.625 ms 候选重启，不把脚本已写完算作大步已经验证。

从真实 46.500 ms 父态另行登记 62.5/31.25 us 同父局部细化，可以判断该新短段的
时间差异；即使通过，也不会证明 0–46.500 ms 整条旧历史全局收敛，更不会洗白上述
两个 FAIL。完整沉积仍受方法和有限域共同限制：当前 4 mm 基板初始即小于 5.4 mm
肩径，eligible 库存有限且没有持续入口补料。

## 精确证据绑定

以下 SHA256 为本报告读取时的真实文件值；路径均相对仓库根。

```json
{
  "scripts/screen_visible_top_timestep.py": "cfcb4bbc15c95138abdfba0ac8654e754c231256eeca501c859cdc5252b29bd8",
  "scripts/screen_visible_established_contact_timestep.py": "b1ff0c3e7927ee8262850d9f05a8a9f7c4aa409b01a3d623a575a8d4decf1ce2",
  "docs/VISIBLE_ESTABLISHED_CONTACT_TIMESTEP_CONTRACT_20261004.md": "04e26cb0af0ecf1ef976464c91ec9effccd5dfa3ac460e95157609fa458c6996",
  "scripts/audit_saved_particle_temporal_difference.py": "4add33064507ba55f6f8fe647f7b848a666dc3707da002ad0385c9215d376f3d",
  "outputs/visible_top_timestep_20261004/mult_04/manifest.json": "a3fce98a626ce893cc58679e953e5c988aa0c128f35ae93901a81a091bd99caf",
  "outputs/visible_top_timestep_20261004/mult_04/step_001/temporal_screen.json": "e03030bad706c0f9efb4f36feaac93e2f2b6efda5a451261288b4c9a35af5b5e",
  "logs/visible_top_timestep_20261004/mult_04.log": "d56b07a62bd21997621b1bd7b4058899b057ec7fb391f07e16330735ff143afb",
  "logs/visible_top_timestep_20261004/mult_04_supervision.json": "5b749e3eb27096ea78f6cac8400c678363b2d6b3fe5e8f047f24b125154da770",
  "outputs/visible_top_event_20261004/window_002/fine/manifest.json": "bb3c2c900fce889d9e4ede003d071cb19dbd77b55f20787df145ad000cbb7fa6",
  "outputs/visible_top_event_20261004/window_002/fine/pair_014_temporal_screen.json": "49eb747b2bc18fd69aed48602c578e10fd356844062ced27d8fd040b3a9baae6",
  "outputs/visible_top_event_20261004/window_002/fine/step_026/accepted.npz": "b395158b765378f4dcd64f2fa7a98d236efb8a260be79a385d3acdf71d4be3ce",
  "logs/visible_top_event_20261004/window_002/fine.log": "e95abdaa2080c915666ceffa42a49d1b2bd024c521e83a841ec743728648fcba",
  "logs/visible_top_event_20261004/window_002/fine_supervision.json": "50eb53be810862b4d4bd12f74123b78094114e0c81e72e66377a046498aa7356",
  "outputs/visible_established_contact_timestep_20261004/window_002_reference_blocker_pid270.json": "59951793c2972ac0e0217217adb0dfe8cd3d71993efcfff10a5632a4e29a5c69"
}
```

PID270 诊断 JSON 另绑定全部 14 个 pair screen、4 个实际候选 checkpoint、接触 journal
和保存的邻域图 SHA；其 PASS 仅表示已存量读取与算术完成，不是物理/因果验证。
