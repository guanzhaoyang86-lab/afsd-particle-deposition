# Predictor + Onset 全路径输出独立复核（2026-10-02）

状态：PASS_INDEPENDENT_SAVED_ARRAY_AUDIT_NOT_AFSD_VALIDATION。
这是现有输出的独立 NumPy/JSON 后处理，不导入或运行求解器、接触分支、材料更新器。
没有再跑真实路径、材料点或测试。使用 scripts/audit_process_predictor_onset.py；
完整机器结果在 logs/process_predictor_onset_20261002/independent_audit.log。

## 输入、可重跑命令与完整性

新目录 outputs/process_predictor_onset_20261002/n3_i4；
旧目录 outputs/process_plastic_friction_onset_20261002/n3_i4。
按各自 report 清单重算新63、旧34个输出文件哈希，原基线来源哈希及新
predictor_provenance 的23个来源条目全部一致。

```powershell
.envs/pysph/python.exe -B scripts/audit_process_predictor_onset.py --run outputs/process_predictor_onset_20261002/n3_i4 --baseline outputs/process_plastic_friction_onset_20261002/n3_i4
```

实际上述命令 stdout/stderr 通过 Tee-Object 保存 independent_audit.log，exit 0。
只读，生产/物理求解次数为0。新输出40接受步、41帧、24事件、终时0.2s；
旧输出16接受步、17帧、终时0.1125s。

## 旧物理前缀逐字节不变：PASS

history.npz 全30字段逐数组检查 dtype、shape、C-order 数据字节：
8个静态字段整体相同，22个时序字段的前17帧完全相同。
包括位置/位移、应力、塑性应变与功、法向/切向力、相对滑移、signed Q、
冲量、运动学速率、模式、ID和标签；不是仅 allclose。

accepted_0000 至 accepted_0016 的全部物理数组也完全相同。
初始文件所有字段一致；step1–16仅 record_json 的 iteration_trace 不同：
新轨迹明确标 CONE_ROOT、法向力清单等诊断信息。解析 JSON 后，所有其他
标量物理记录一致。压缩 npz 整文件或 record_json 字符串哈希差异因此不是
物理前缀变化，不将元数据不同掩盖成整文件 byte-identical。

## 接触过程、初猜使用与预算：PASS

独立遍历16个 requested attempt，886个返回分支根，40个连续接受记录。
全路径仅1个分支实际使用 predictor，仅1次预测修正：

- 旧accepted step16 / 0.1125s，TARGET 0.125s。
- 初始最小 N=-6.972110549861 N。
- 预测1＋原cone4＝总5修正，返回根最小 N=14.599258321040 N。
- 该目标根没有直接提交；原事件二分随后逐个定位八个接触起滑。
- 886个分支的共享修正最大值14，未出现14预测＋14cone的双预算。

首次越过旧末态的事件是 step17 / 0.11288886198308318s，ID19、feature1；
ID27是已释放接触。以下是0.125s请求内新接触到达摩擦边界的顺序：

| 接受step | 时间s | 新到达边界粒子ID | feature |
| --- | --- | --- | --- |
| 17 | 0.11288886198308318 | 19 | 1 |
| 18 | 0.11321856825970292 | 24 | 1 |
| 19 | 0.11397415501276222 | 26 | 1 |
| 20 | 0.11525640095653489 | 22 | 1 |
| 21 | 0.11553810685938190 | 20 | 1 |
| 22 | 0.11635841864121188 | 21 | 1 |
| 23 | 0.11711618621028780 | 25 | 1 |
| 24 | 0.11742992799129748 | 23 | 1 |

这里“到达边界”由 event 中未释放且 abs(mu*N-|tau|)<=原9e-5N判据识别，
不等同于该结束区间已发生有限滑移。原metadata的模式和released原样保留。
后续step29、30、32–38分别对应ID19、27、22、24、20、26、21、25、23。
新增17个事件 bracket 宽度均<=1e-12s；逐事件完整ID/已释放ID在日志。

## 从数组独立重算 signed 功：PASS；材料能量闭合：仍警告

以持久ID连接接触行与粒子位置，独立使用
Q=-0.5*(tau_old+tau_new)·(du_contact-dw)，
Wwall=sum(mid_contact_force)·dw，Wmaterial=sum(mid_contact_force·du)。
没有调用 vector_contact_work，没有 abs(Q) 或 max(Q,0)。

- 重算累计 Q=0.043185381195100875 J。
- Wwall=0.34808659528485164 J；Wmaterial=0.3049012140897508 J。
- 支承功0；最大逐区间 Wwall-Wmaterial-sum(Q)=1.605e-17 J。
- 重算Q与保存interval Q最大差4.337e-19 J。
- 存在25个负的原始Q分项，最低-9.968003103622417e-20 J（step33、ID24）；
  仍原样保留，没有剪成0；该尺度远低于原负耗散1e-10J门槛。
- 重算接触/支承累计功与保存应力功最大差2.0372e-11 J。

最后一项仅比较接触账与保存的应力功；本后处理没有独立重新积分本构。
保存的末态 Ws=0.3049012140693836 J，Es_diag=0.05382773794561942 J，
P=0.2554214501459417 J，独立相减 Ws-Es_diag-P=-0.004347974022177514 J。
step2–40共39步材料警告均仍存在。接触功恒等式通过不意味着材料能量闭合，
未把该差值加到摩擦热、塑性功或其他储能项。

## ID、有限值和本次复核边界

全部history字段有限；27个粒子ID唯一，接触ID19–27唯一；所有41个accepted
存档ID/标签一致，restart_allowed=false。静态总质量0.0027000000000000006kg。
接受后的最小N（不含初始）147.02708804393544N；最大接触根残差
6.828054982765814e-5N；最大ep=0.004181979456166672。

NOT_RUN/不授予认证：材料能量修复、正常开离/再接触、热、FKS、旋转工具、
真实AFSD沉积、8增量敏感性与空间收敛。0.2s属于制造准静态小块路径，
不是论文AFSD沉积时长；本次完整路径不取消旧材料/空间/闭循环FAIL。
