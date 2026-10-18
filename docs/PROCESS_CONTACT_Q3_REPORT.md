# Q3：刚性平压头压入—脱离结果（2026-10-01）

## 本轮新增的物理能力

独立三维可变形小块现在由单边接触受压，压头卸载后可以脱离；不是Q2的双向端面
位移约束。J2两档增量都实际出现零反力释放和残余塑性形状，脱离后材料不随壁抬升。
这是ENGINEERING_APPROXIMATE、固定参考小应变、等温制造基准，不是AFSD沉积。
PySPH commit仍1cd3dea9cfaac7d9f1bda2a2491ce4b42a31dbe5；src/旧配置测试/旧失败未改。

先验合同：PROCESS_CONTACT_Q3_CONTRACT.md；独立推导：PROCESS_CONTACT_Q3_REFERENCES.md。
新实现scripts/probe_process_contact.py继承Q2的材料/梯度/切线设施；新reader
inspect_process_contact.py核对保存数据。配置configs/benchmarks/process_contact_q3.json。
命令/预算/完整日志：logs/process_contact_q3_20261001/COMMANDS.md；
实际数据：outputs/process_contact_q3_20261001。各完整算例仅运行一次，无参数重试。

## 调用与边界

    Q2初始化/PySPH粒子/三维NNPS -> reset permanent rollers only
    advance(requested_time)
      -> solve_candidate: tentative active contact set
        -> _equilibrate: Q2 trial/material.update -> internal force
                        Q2 tangent -> reduced GLOBAL Newton response
        -> remove tensile contacts / add penetrated contacts / re-equilibrate
      -> if loaded contacts would open: locate first zero-force release on branch
      -> commit: recompute material/contact/KKT/work/metadata -> accepted state
      -> save event frame -> solve remainder of requested physical interval

27点、三轴自由度；初始54自由DOF，9个面节点接触时仅其法向受约束，另有45自由DOF。
永久支承只有X=0/Y=0/Z=0各自法向roller，不将全块或整个cap刚化。
壁初始x=.01m，规定0→-40um→0的位移，各阶段.1s。材料无参考轨迹输入。
条件为g=d-u_x>=0、lambda>=0、lambda*g=0，压头对材料力=-lambda*e_x。
lambda由真实内部力恢复；不把拉力裁成零并保留拖拽位移。
实际解的是活跃约束后的整体切线子块，含粒子间耦合，不用单点P/m；
理想塑性全自由切线可能奇异，不能把它等同于受控接触子块必然奇异。

脱离事件由有界试探定位，不把独立解析时间写入求解器。
事件子步同样走材料/平衡/提交并输出；剩余区间重新求无接触平衡。
两个永久支承横向分量仍存在，所以这里不是完全无支承自由物体释放。

## 实际结果

| 算例 | 每相请求增量/实际接受步 | 终点 s | 完整进程墙钟 s | 释放事件 s | 终轴向材料面位移 um |
| --- | ---: | ---: | ---: | ---: | ---: |
| elastic_n3_i8 | 8/16 | 0.2 | 0.958933 | 无内部事件 | 0 |
| j2_n3_i8 | 8/17 | 0.2 | 2.157425 | 0.132142857142753 | -27.142857142899 |
| j2_n3_i16 | 16/33 | 0.2 | 2.076379 | 0.132142857104918 | -27.142857158033 |

三组UNILATERAL_PATH_COMPLETE，exit0；三组读回PASS_COMPLETE_READBACK，exit0。
elastic仅回到零位移/零反力，未有严格正gap的释放后段；J2两组终gap约27.14um，
全部接触解除、接触力为0。实际初接触时刻为0，本轮没有任意再碰撞验证。
终点都是本基准0.2s，不是AFSD工艺已推进0.2s。

| 算例 | 累计壁功 J | 累计材料接触功 J | 塑性功 J | Wstress-Es-P J |
| --- | ---: | ---: | ---: | ---: |
| elastic_n3_i8 | -3.643e-17 | -3.123e-17 | 0 | 3.469e-18 |
| j2_n3_i8 | 0.242142857144 | 0.242142857144 | 0.244285714286 | -0.002142857143 |
| j2_n3_i16 | 0.244017857132 | 0.244017857144 | 0.244285714286 | -0.000267857143 |

J2两组最大单区间壁功—材料接触功差分别3.949e-14/1.138e-11J，均实际检查。
保留返回算法诊断0.068035714286/0.035223214286J，不作为热或总能量缺口。
支承不移动，支承功实测0。未求温度，不将P或账差自动转热。
本基准回读可从原始反力/位移独立重算各区间功，而不仅信报告中的PASS。

### 独立数值解释

制造J2为E70GPa、nu0、Y90MPa、rho2700kg/m3，非AA6061生产标定。
独立连续体参考：峰反力-9000N，终ep=.002714285714286；释放时间
.132142857142857s、残余轴向面位移-27.142857142857um，累计塑性功.244285714286J。
实际两组体积平均ep为.002714285714284607/.002714285714285692，释放后ep/P不再增长。

独立分段解析应力在8/16节点的加载梯形功为.300/.301875J；卸载段释放后不做功，
所得净功恰为.242142857143/.244017857143J，解释上表未补热的屈服拐点积分差。
同终点将时间增量减半，功差绝对值由.002142857降为.000267857J；
不需要为了这个采样误差重写材料或继续细分到全绿。

事件二分括区间宽度两组均约7.276e-13s，小于事前1e-12s。
这不等于相对解析事件时间的误差一定小于1e-12s：细增量组约偏3.79e-11s，
因为事件分支还受整体平衡停止余量影响。原始终应力最大值分别.00159/.47524Pa，
没有裁为零；细组终平衡relative2.587e-9仍在事前1e-8内。
弹性终某节点原始lambda=-9.70e-14N，保留舍入符号，不称数学上逐点精确非负。

## 回归、最小修复和历史

核心13项PASS：q3_tests_final.log，2.435s，完整进程2.831484s，exit0。
包含独立整体2DOF凝聚刚度1.5N/m（非对角2N/m）、约束后子块条件数、54自由度、
试探/切线/活跃集不污染、过期/不平衡候选拒绝、提交重核、弹性路径、
J2解析释放/残余形状/功、真实标签变更不影响力学，以及配置来源一致性。

原q3_tests_initial.log为12项中11PASS/1FAIL：RELEASE_LIMIT缺事件metadata时
在改变accepted后才抛异常。只将事件kind/bracket检查和拷贝移到所有提交之前，
修复后q3_tests_fixed.log12PASS保留；没有为过关改力学、材料、物性或容限。
随后发现标签测试使用index%3恰好等于原z%3标签；改为71+index%5并先断言确实不同，
增强测试有效性，不把原无效标签覆盖当证据。追加配置漂移检查后为最终13项。

还加Q3声明参数必须与实际Q2使用值一致、production必须false的快速拒绝，
以及REJECTED清单的candidate_time/dt/base_step；这些均在完整算例前实施。
未调用被测接触函数生成解析reference；完整模型2D/3D的比较不在本轮伪造。

读回回归6/6 PASS、0.254s、完整进程.571578s、exit0（readback_tests.log），
本轮最终新测试13+6=19项PASS；验证原始记录、截断账本、
篡改字段/哈希后物理定义、ID/标签、BLOCKED不晋升完成等，不把读回当新的求解器验证。
复用Q2的23项本轮重新执行PASS（q2_compatibility.log，0.774s，exit0），与新测试分开计数。
旧主116项本轮NOT_RUN，旧111PASS/5FAIL和全部失败数据保留，不被上述PASS覆盖。

## 性能、范围与下一步

三组advance耗时.571164/1.775138/1.701000s。全程单线程、外timeout60s；未启动长任务。
selected_array_bytes235764，不含临时切线/求解缓冲和RSS；峰值进程内存NOT_RUN。
这不是显式CFL步速或宏观AFSD成本估计，不可由此宣传完整AFSD的加速倍数。

本轮Q3空间敏感性NOT_RUN（沿用n3）；Q2 n5约22.11%反力偏差仍是离散风险。
新接触的候选只有初始平面节点、固定法向和小变形；一般碰撞、角部多面切换、
有限变形/本构自由能、旋转摩擦、FKS接触、热、开放方孔和供料尚未接入此准静态变体。
旧Jaumann闭循环反例不被这条固定主轴路径消除。constructor异常仍可能无REJECTED文件，
只在形成candidate后承诺可得候选的时间层/数组；checkpoint不支持生产重启。

下一小单元Q4：两个相交刚壁的法向接触和脱离，先独立耦合参考，再三维角部，
每个壁独立记录力/间隙/功、共享整体材料切线和accepted状态；不并行加入摩擦和热。
预计扩展新Q3 companion的接触几何/法向映射接口和专用测试，src生产默认仍不改。
若双面法向耦合可行，再进入短切向运动/摩擦及实际工具几何；当前这些均NOT_RUN。
已足够进入上述下一机制开发，不授权三工况生产、长跑、GPU或完整复现标签。
真实AFSD仍是原2us/191步、无持续沉积；后续绘图只能读实际冻结字段，不能改物理。
