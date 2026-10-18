# 方法类验收：接触认证与本构/边界决策（2026-09-30）

## 结论

本轮完成的是**固定范围的局部接触方法认证**，不是又修过一个粒子就跑AFSD。
432个预定义凸锥子问题均通过独立证书；439的闭式锥参考和完整冻结kick都通过。
固定域边界探针仍不支持节点一致性；本构闭循环非零功仍在。
APPROXIMATE/BLOCKED；真实2us、时空敏感性、完整三维耦合新基准和生产均NOT_RUN。

PySPH仍为1cd3dea9cfaac7d9f1bda2a2491ce4b42a31dbe5；复用原设施，无新依赖、
GPU/并行/性能工程或第二后端。本轮不改material.py、mechanics.py、solver.py、
geometry.py、thermal.py及任何物理配置。默认legacy仍保留，显式event_aware_v1
使用改进数值求解；完整源码hash变化禁止直接恢复旧源码checkpoint。

## 1. 先验合同与基线

docs/CONTACT_CLASS_CONTRACT.md为验收依据；初稿副本logs/contact_class_20260930/contract_pretest.md。
修改前360项旧冻结PASS；日志freeze_before.log。旧代码和错误唯一面力验收快照在source_before/。
全部本轮日志根为logs/contact_class_20260930，产物另在outputs/contact_class_20260930。

新测试首先暴露两项测试参考实现问题：恒零目标的极小分母溢出、把理论零切向
迁移率的舍入噪声归一为O(1)。baseline/和baseline_valid_oracle/均保留。
修正独立参考尺度后，**baseline_final/**仍使用完全未改的被测源码：

- 固定432个锥问题：372通过、60不收敛。
- 原439锥：FAIL。
- 共同法向速度从0变为1e14m/s的数学不变性测试：旧算法把x末速-2/3算成0，
  摩擦热从7/18J变成-0.5J，漏掉反向滑移。此极端boost仅用于坐标不变性测试，
  不是AFSD运行速度或参数。x相对滑移仍精确可表示。
- 原454和12项事件测试PASS，说明它们单独通过不能覆盖整类问题。

## 2. 统一数学问题、独立参考和最小源码范围

未知每面二维切向力f_a，三维合力sum(T_a f_a)。求解原库仑圆盘上的凸问题：

    phi(f)=b:f + c/2 * (sum T_a f_a)^T M (sum T_a f_a)
    ||f_a|| <= mu Fn_a,  M=P/m
    H_ab = c T_a^T M T_b

H的秩最多3，而未知量为2N；多面分力不一定唯一。独立gap证书为
sum(g_a dot f_a + limit_a*||g_a||)，其中g=b+Hf。由圆盘支持函数可得其为
可行解相对最优值的上界；同时检查独立投影驻点残差、可行性、已知目标值及末速响应。
不以“迭代两次变化很小”代替验收，不由求解器输出生成reference。

制造参考先选f*和非负边界乘子lambda，再设b=-Hf*-lambda*f*，解析保证KKT。
覆盖1/2/3面、0/1/2/3自由度、内点/混合/滑动、近乎平行与正交面、
旋转/排序、容量比及1e-6/1/1e6数值尺度。432是子案例数，不与unittest数量相加。
线性/零容量、非法输入、不同墙速、共同boost及事件反转另测。

改动文件/接口：

| 位置 | 改动与限制 |
|---|---|
| contact_cone.py::solve | 新增局部缩放投影KKT/Newton求解，最多200步；原圆盘/迁移率/目标不变。显式秩揭示最小二乘仅求Newton修正，证书检验未正则化原问题，不是空间M_i伪逆 |
| contact_kick.py::cone_solve | 接口保持返回三维面力；旧算法原样保留为cone_solve_projected_legacy；不添加物理选解正则项 |
| contact_kick.py::kick | 在一个kick内选平移速度参考，避免绝对实验室速度扩大停滑阈值；输出仍LAB。动能差改等价点积表达以减少平方相减消去 |
| 同文件::interval | 原始Q=-Jt·u_mid不变；记录舍入预算、速度尺度和参考帧，不裁负值 |
| test_contact454_precision.py | 两面静力分配不再强认(1,1.5)N唯一；检验圆盘/切向/合力。原文件完整备份，并新增(0.75,1.75)N也合法的独立证书回归 |
| test_contact439_open_risk.py | 独立唯一解原容限不变；成功/新证书失败都可落盘 |
| test_contact_class.py / run_contact_class.py | 固定矩阵、独立判据、每命令60s上限、每次新目录；不启动真实AFSD |

求解器内部证书目标2e-14，未放宽原请求2e-12；439/454唯一参考的原力容限未变。
一般类证书与可行性容限在运行前定义，不能拿不同指标直接当旧力增量的宽松替代。
实际通过只覆盖这组有限矩阵，不是对任意输入的全局收敛证明。

## 3. 两轮实现与原始负热的处理

第一轮已通过432矩阵、439和事件测试，但旧曲线停止测试驱动/容量=.5失败：
Q=-3.91653943e-19J，区间速度约1e-16m/s，旧瞬时舍入预算缩至6.11e-31J。
这是参考帧消去后遗漏初始速度累积舍入尺度；细节在round1_diagnostic.log。

第二轮只将不随共同boost变化的已有特征速度加入256eps*|Jt|的舍入预算，
没有改变原始热值、停滑阈值、物理摩擦力或加残差热。新回归独立重算该预算和
冲量中点功；实际最小原始热仍为-3.91653943e-19J，两种参考系均保留该负值。
相对地，旧boost测试的-0.5J不再出现，正确热为7/18J；没有用预算掩盖O(1)错误。
两轮日志implementation_round1/、implementation_round2/均保留。

最终类级日志：final_class_with_holdout/；完整回归：final_suite_with_holdout/。
后者116项：**111 PASS、5 FAIL、0 ERROR，exit1**。

| 层级 | 实测状态 |
|---|---|
| 固定432锥矩阵及9项类级测试 | PASS |
| 439独立闭式锥 | PASS，力误差4.33579e-16N |
| 439完整冻结局部kick | PASS完成/冲量/功恒等式；没有声称存在完整轨迹的解析参考 |
| 原454精度/旋转/停止/热5项、原事件12项 | PASS；1D解析反转仍为t=1/3、u=-2/3、Q=7/18 |
| 旧raw旋转、raw异密度功、legacy标量反转 | 历史路径3 FAIL，继续保留 |
| 当前材料Es/自由能一致性、节点边界patch | 当前路径2 FAIL |
| 闭循环诊断脚本 | PASS仅表示确认了非零残功，绝非材料模型合格 |
| 完整三维耦合新基准 / 真实2us / CFL和空间序列 / 生产 | NOT_RUN，未过后续门禁 |

完整回归包含原有单步/短链路smoke测试，不是本轮新AFSD真实时序试跑。
不能据PASS数推导论文完成度。源代码设施/三维自由度/标签/ID回归仍通过。

## 4. 本构方法决策

当前constant-G Cauchy/Jaumann更新在指定闭合路径上，功的独立值为
(1-Jc)G(A:A)e²V0=0.146111111111J；16/32/64/128子步均趋向此非零值。
当前闭循环Es峰值14.611111111J；此特定路径残功比例为1-Jc=1%，
不能把它推广成任意AFSD轨迹的误差上界。

决定：保留现有论文近似路径、不替换本构、不补偿热；其有限应变自由能和完整能量
验收仍BLOCKED。原论文方法复现与能量一致方法扩展分开声明，禁止为全绿换模型。
FINITE_STRAIN_VARIANT_DESIGN.md中的自由能候选继续DESIGN_ONLY/NOT_RUN。
下一个材料决策必须说明论文模型适用范围/误差证据，或先完成独立材料点自由能
变体并明确方法偏离；未决FKS和热学来源不能因此自动补齐。

## 5. 边界方法决策与固定物理域证据

scripts/probe_method_boundaries.py保持[-1,1]^3立方体、常应力、体积8、
相同物理牵引sigma*n；仅改粒距，h/dx=1.3。外力由独立梯形面面积积分生成，
不是取负的被测内力。仿射v=x和三次v=(x³,y³,z³)的解析表面功均为48。
输出：outputs/contact_class_20260930/boundary_probe/report.json。

| 每轴粒子数 | 节点相对残差 | 仿射内部功 | 三次内部功 | 三次相对误差 |
|---:|---:|---:|---:|---:|
| 3 | 0.32615719 | 48 | 48 | 0（采样混叠，x³=x于-1/0/1） |
| 5 | 0.57263126 | 48 | 40.06580603 | 0.16529571 |
| 7 | 0.52918142 | 48 | 42.76922812 | 0.10897441 |

结论：仿射弱功检验仍通过，但独立牵引的节点相容性没有被验证；不得用粗档混叠
零误差宣称准确或据三点声称收敛阶。5到7的弱功误差下降是有限观察，不是完整一致性证明。
这不直接证明corrected/伴随算子全部错误，也不能把当前节点FAIL改成PASS。

决定：继续保留corrected/伴随内力成对路径、cond(M_i)>1e6失败和原patch反例。
下一阶段工作单元是**离散边界功与牵引施加的配套设计**：先给离散弱形式和独立
常应力/多项式基准，再实现最小边界接口。不得只令外力=-内力让patch通过，
不得修改点云边界/物理牵引来追求结果。域内力mechanics.py、边界加载/测试接口
是候选改动位置，尚未授权为另一套空间方法的静默替换。

## 6. 下一关及可重跑命令

在项目根、沿用.envs/pysph/python.exe；输出目录必须新建：

```powershell
& .envs/pysph/python.exe -B scripts/run_contact_class.py --output logs/contact_class_recheck_01
& .envs/pysph/python.exe -B scripts/probe_method_boundaries.py --output outputs/boundary_recheck_01
& .envs/pysph/python.exe -B scripts/run_coupling_review.py --output logs/method_recheck_01
```

第一条当前exit0；第二条exit0仅为诊断完成，报告显式NOT_METHOD_CERTIFIED；第三条exit1保留5FAIL。
如sandbox账户需要git只读权限，沿用CONTACT454_REPAIR_REPORT.md的进程级safe.directory设置，
不写系统配置。本轮未安装依赖，隔离环境没有SciPy，局部方法仅使用已有NumPy。

下一步不再追新的粒子ID：先解决边界弱形式/牵引相容门禁，并完成本构适用性决策，
再定义非AFSD三维变形—刚壁摩擦—导热的耦合验收。那一新基准当前NOT_RUN。
其通过后才讨论固定AFSD2us与同终时敏感性；每次600s上限不变。
旧失败与全部三个旧freeze不覆盖；新METHOD_CLASS_FREEZE.json只证明证据完整，
绝不授予生产或论文复现标签。
