# 整体非线性停滞的独立只读审查（2026-10-03）

范围：阅读实际源码和既有32次Picard日志，独立处理最后保存的候选数组。
没有调用UpdatedTrial/contact/material来计算新物理候选，没有运行新的真实时间步，
没有修改任何已冻结物理、容限、材料、配置或原始失败证据。
本文件不是新求解器的测试通过证书；后续实际运行以主代理报告为准。

## 已确认的调用和确定性边界

`UpdatedTrial.evaluate`确实从同一接受历史克隆PySPH粒子，在候选中点重新求NNPS、
修正矩阵及平均速度梯度；从原S/ep更新材料，再求中点应力内力、事件接触、热和五场map。
没有在迭代间累积ep/热。已有纯性/重排/换标签短测试覆盖该接口；本次只读审查没有重跑这些测试。

- `candidate.velocity_m_s`只影响域检查，不进入实际力学map；其根是返回端点速度。
- 位置根满足`x1=x0+dt*vbar`，但迭代中位置和vbar独立，位置仍影响几何和接触。
- 因而五场历史容易有相关列。裸、无尺度、无秩检查的Anderson正规方程不合适。
- FKS对`max(rate,rate_floor)`的对数在`rate_floor=1/s`处连续，但导数有折角。
  塑性active切换、接触active/滑移及有限feature也使整体映射分段光滑，不具备全局光滑保证。
- NNPS集合跨核支撑可能改变，但现有核在q=2的函数/一阶导数为零；
  不能把所有邻居表变化直接认定为物理跳跃。有限face丢失非零储能仍由旧守卫拒绝。

源码检查没有找到可证实的同史污染、隐藏裁剪或一个明显实现错误能够解释全部平台。
这不等于证明完整映射连续或正确：未知局部接触误差、有限feature及原本构局限仍存在。

## 新的具体证据：塑性点紧贴速率折角

可重跑命令（只读，约1秒外进程；输出必须是未存在文件）：

```powershell
.envs/pysph/python.exe -B scripts/assess_nonlinear_continuation.py --output outputs/nonlinear_continuation_audit_20261003/assessment.json
```

本次结果：`PASS_READ_ONLY_DIAGNOSTIC_NOT_PHYSICAL_ADVANCE`，内部0.165617秒；
`evaluated_new_physical_candidates=0`。独立数组梯度与原保存梯度最大差为0。
输入、脚本SHA和详细值保存在上述JSON，全部输入前后哈希相同。

最后完整候选的速度主导残差仍为PID268的z方向`5.0427445671e-4m/s`。
它没有直接接触、ep=0，是PID269的真实NNPS邻点。
PID269无直接接触，已塑化，ep=`3.2768391479e-5`，rate=`0.9913406664/s`，
恰在FKS的1/s低率正则折角下方。

本次仅冻结保存的几何、rho、B/g和最终ep，计算
`L(c+alpha*(map-c))=L(c)+alpha*L(map-c)`，再用独立FKS直接公式计算Y，
不调用被测函数，不将此线性方向诊断当作完整trial或新接受态：

| alpha | PID269 rate，1/s | 固定ep下Y，Pa |
|---:|---:|---:|
| 0 | 0.9913406664 | 1,517,301.2707 |
| 0.10 | 1.0058190600 | 1,541,432.0391 |
| 0.25 | 1.0275867121 | 1,630,477.7441 |
| 0.50 | 1.0639905202 | 1,775,263.4277 |
| 原Picard 0.75 | 1.1005364354 | 1,915,714.7934 |
| 1 | 1.1372107573 | 2,052,047.1084 |

该点T*=0.76512457，`b1*T*^b2=2.74098736`，Y在floor正侧的导数约
`4.15890361e6 Pa·s`，负侧为0。不到`2.6e-4m/s`的候选均速缺陷，经邻域梯度
即可跨折角；原0.75方向的固定ep屈服变化约26.3%。这不是对原材料的调参建议。

限制：原运行只存最后完整candidate，32行历史只存残差范数，不能证明每一轮的有符号
两周期或全Jacobian谱。本诊断也没有更新方向上的几何/rho/ep/接触，因此只支持机制，
没有证明根的存在、唯一性或新整体求解必能收敛。

## 根因分级

| 级别 | 结论/证据 |
|---|---|
| CONFIRMED | 固定50ns、0.75 Picard在32次预算内失败；最后残差远超过冻结根精度，不是最后一个ulp问题 |
| CONFIRMED | 最大v残差点无直接face；邻点269塑化且位于已声明FKS低率折角附近 |
| SUPPORTED | 同步材料率依赖与邻域内力反馈能把Picard候选跨越强切线折角；整体非线性策略比继续修局部摩擦面更针对证据 |
| POSSIBLE | 两周期、接触远程耦合、多个塑性分支或局部切线失稳；没有全谱证据 |
| RULED_OUT（本次停止原因） | 达32次前没有负摩擦异常、NaN、温域/密度/穿透/条件数拒绝或超时 |
| BLOCKED/NOT_RUN | 任意新整体根、真实多步/输运、全局Jacobian、真实时间步收敛、过程尺度成本 |

## 最小整体策略建议，不改本构

1. 保留原纯trial及所有原始残差精度，新增有限历史、尺度化、秩保护的全局Anderson或
   类secant求根。x/rho/T/v量纲必须缩放，物理单位和验收门不变。
2. 记录每次真实evaluate、残差、候选、是否用于历史、秩/条件数；只接受原始map残差
   达标的实际候选。不得把混合系数、拟合误差或小步长当收敛。
3. 秩不足历史显式丢弃/重启，禁止静默把病态空间Mi改用pinv。
   总evaluate/墙钟预算包含失败线搜索；预算耗尽保存REJECTED，不能偷续32次。
4. 对非光滑/域拒绝，丢弃越域候选并明确记录；必要时回到上一个已评估候选做有界
   残差驱动回退。候选回退不是物理时间缩放，不修改dt/模量/质量。
5. 明确记录PID268/269的有符号残差、269的rate/active/ep及active contact集合，
   使下一次失败不再只有范数而无法辨认分支周期。

独立新求根基准宜包括：负斜率线性map、尺度相差多数量级的五场、秩亏历史、域拒绝、
以及连续折角制造map `g(z)=1.1-9*max(z-1,0)`，其唯一根是`z=1.01`。
此max只定义数学回归，不进入物理FKS或耗散；参照根可直接代数推导，不能用被测算法生成。

完整有限差分Newton若按342点五场共3762未知逐列取导数，按上一轮完整候选约3.81秒估计，
一张Jacobian就约3.98小时，远超600秒，不建议当前实现。
低维矩阵无关Newton/Krylov可以作为后备，但接触/材料折角、有限差分尺度与预条件尚需
独立证明，不能默认比有界历史法便宜或稳健。不得启动未登记成本的新真实物理试算。

连续推进最终还需实际接受根、可重启同变体状态、连续多步及真实材料输运度量；
单个50ns根成功只能解锁这些工作，不能直接授予75%或宣布AFSD沉积已经完成。

## 新增连续driver/提交接口的只读审查

对象：`continuation_state.py`、`run_nonlinear_continuation.py`、`stage_anderson.py`，
并跟读实际旧acceptance、checkpoint、thermal ledger和transport observer。
本节没有执行新物理trial或测试；测试结论需看主代理保存的日志。

实际提交顺序为：同史求根 → 原完整acceptance → 独立克隆commit → 按ID输运观测
→ source/input复核 → 新状态累计账 → checkpoint/VTK写出 → 严格NPZ重载与全属性hash
→ acceptance/transport/tool_pose/accepted_trial落盘 → 再次source/input复核
→ accepted_history发布 → 以重载PA而非旧候选开始下一阶段。
在成功路径上，没有在求根中累积塑性或把未验收candidate提供给下一阶段。

### 审查发现并已反馈的具体问题

1. **真实累计错误，源码已修：** 初版runner将boundary ledger每个数值相加，
   误累计normal potential、最大穿透和当前contact数量。现`accumulate`将这些字段
   放入latest状态，穿透另取全史最大，原始摩擦最小值另取全史最小；只累计J/冲量/事件次数。
   thermal ledger经逐键阅读确实均为J增量，不含热行或时间步界，允许累计。
2. **失败evaluate计数，源码已修：** 初版异常candidate不走journal，导致总evaluate漏记1。
   现Anderson异常路径调用`on_evaluation(row,c,None)`，runner保留异常candidate并计数。
   该错误原本会立即停止而不会暗中多跑，但成本/失败证据不应漏记。
3. **累计范围必须明确：** PA的plastic_work/return从真实step62继承，原IntegratedSolver
   实际按J累加，并非per-volume scratch；新continuation ledger从0起。
   应明确ledger是新变体continuation区间账、逐粒子累计包含parent，二者不得直接混作同起点。
4. **发布事务建议：** accepted_history应先将`history+[entry]`原子写出成功后才更新
   内存history；不能先append再写，在磁盘异常时误报accepted_steps。
5. **恢复授权尚有边界：** 底层restore核完整文件、source/config、全属性/ID，但不核公开
   accepted_history发布记录。它用于发布前内部往返是合理的；不是独立公开的安全续算入口。
   当前CLI没有`--restart`，不得根据一个孤立accepted.npz文件就宣称已验证外部恢复。
   全源同变体重新求解的restart等价仍NOT_RUN。

### 预算与证据边界

源码真实限制最多4阶段、每阶段32次、整个进程72次，并对每步按remaining取较小上限；
因此不是4×32无限追加。原物理dt50ns不变。新probe包括完整evaluate和压缩写盘，
同source/config才能授权`1.5*耗时*72+30<=450s`的一次进程；内540s，外部仍必须570s保护。
成本probe不能因为其candidate完整就计为一个物理步。

每candidate已补268/269有符号速度/T缺陷、rate减floor、塑性增量，后续可判断实际分支，
不需要再凭单个残差范数猜测。没有改原FKS floor或热/摩擦模型。

累计边界normal_potential来自冻结中点contact ledger，而非接受端点；解释其latest字段时
仍须保留这个时间层。接受端点U和normal temporal residual另在acceptance记录。
独立完整轨迹/输运核验还须在真实运行后完成，不能由代码审查自动标PASS。

## 唯一真实continuation完成后的独立读回

新增` scripts/assess_continuation_result.py `（仅反序列化与数组代数，不调用被测求解/
材料/接触/热/acceptance/observer生成参考）；源输入与输出全部前后核哈希。
只把正式accepted_history中的snapshot计为接受态，孤立文件或REJECTED均不计物理进展。

首次执行直接PASS，未修测试追绿或重跑真实物理：外60秒上限，实测外2.9444321秒、
内部2.499117秒。命令：

```powershell
.envs/pysph/python.exe -B scripts/assess_continuation_result.py --output outputs/nonlinear_continuation_20261003/independent_readback.json
```

详细结果在上述JSON，日志为`logs/nonlinear_continuation_20261003/readback_initial.log`
及同前缀stderr。读取43个登记source、97个产物及全部4正式阶段，没有未发布孤立snapshot。
全属性NPZ、真实字段VTK、工具位姿、ID、粒子质量、x/v/rho/T/S/ep和逐步塑性/return累计对应。
检查全部23,411个实际接触子区间的冲量、原始带符号Q、事件均速和位移map、接触功，
独立塑性热/摩擦热/成对与壁导热、mcΔT、空间功和冻结原始根残差，以及固定z=1.5mm穿越算子。

| 项目 | 实际只读结果 |
|---|---:|
| 已接受连续阶段 | 4，50ns/阶段 |
| 接受终时 | 850.393632910904ns |
| 本新变体实际物理增量 | 200ns |
| 全局evaluate | 61（16/15/15/15） |
| driver内部真实墙钟 | 236.4724705s；外进程由主报告单独记录 |
| 四阶段最大归一根残差 | 0.864020 / 0.151771 / 0.581994 / 0.441446，均<=1 |
| 自由棒相对本轮源态最大净位移 | 3.6174180704e-8m，即36.174nm |
| 向下/向上跨肩面事件质量 | 0kg / 0kg |
| 总真实raw摩擦Q | 2.8814781321e-8J |
| 总新塑性功 | 6.4845985868e-9J |
| 最大接触子区间动量恒等式残差 | 1.3234889801e-23N·s |
| 最大绝对mcΔT总收支残差 | 1.3264935784e-15J |

**塑性范围必须说清：** 仅第1阶段PID269有新塑性增量；在实际step62源checkpoint中
PID269的ep为0，故它是本轮新增塑性点，不是该初态里已塑化点的持续流动。
第2–4阶段均无新塑性增量，整个窗口没有“初态已塑化粒子持续增加ep”的证据。
不能把保存着非零ep的后续粒子绘图称为持续塑性输运。

接触Q中约1e-26J的负舍入项原样保留，落在独立重算的浮点证书内，没有abs或裁剪。
累计状态量修复、history先成功写后发布、ledger_scope以及timeout finally退出均已在
本次真实运行前落盘的源码中确认；此次读回又验证了真实累计账。

结论：`PASS_ACCEPTED_CONTINUATION_READBACK_NOT_DEPOSITION`。
原Picard停滞已被同物理整体求根路线跨过，并获得4个相连接受步；但过程尺度、持续材料
穿越/沉积、自由能认证、同终时敏感性、重新求解的重启等价仍未验证。
不能把36nm短窗位移和0穿越质量解释成所要求的宏观AFSD输运完成。

## 成本前瞻中初始材料上下组的独立纠错

主代理发现其首版只读成本报告误取`min(z0-layer)`为待供入距离；本审查随后独立从
真实step62 NPZ读取component_id/z0/m及原配置，不改solver或reader源码，确认：

- 方棒总150点、质量`5.0625e-5kg`。
- 初始肩面上侧/孔内组`z0>=1.5mm`为75点、质量`2.53125e-5kg`。
- 初始肩面下侧组为75点、同质量；没有中心恰在肩面上的点。
- 方棒原始z层为0.25、0.75、1.25、1.75、2.25、2.75mm。
- 上侧到肩面的最小正距离为0.25mm，下侧的最小正距离同样0.25mm；
  全棒`min(abs(z0-layer))=0.25mm`，而冻结最大总位移只有0.05mm。
- 负的`min(z0-layer)=-1.25mm`描述预置肩下材料的位置，不是进料时间或成本的长度依据。

故“最下方棒中心在肩上0.25mm”是不准确的全棒描述，应明确限定为
“初始肩上/孔内待穿出的那组中心”。原已hash冻结合同的措辞保留并用本节和主纠错文档
更正，不回写删除原证据。首版负时间/负成本无效，不能用于可行性判断。

上侧0.25mm除以2.12mm/s为0.1179245283s，仅是规定执行器速度的长度时间标尺，
不是接触驱动自由材料的真实穿越预测。现有0.05mm原始位置位移guard小于所有中心到
该面的距离，因此两侧任一初始组都不能在此冻结安全域内产生中心cross。
已有肩下75点属于预置初态，不能把其质量当作本轮新沉积质量。

这项前瞻统计错误不改变真实四阶段求解器、独立observer或本次读回结果：
step62实际肩面两侧仍75/75；各相邻接受态按固定面before/after计算的上下cross均为0。
