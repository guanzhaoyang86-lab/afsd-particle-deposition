# E3-2结果：局部数值尺度修复通过，真实2us未重跑（2026-10-01）

## 结论与范围

IMPLEMENTED / LOCAL_REGRESSION_PASS / FULL_2US_NOT_RUN / ENGINEERING_APPROXIMATE。
只有src/afsd_backend/pysph/contact_cone.py的solve实现被修改。
PySPH仍固定1cd3dea9cfaac7d9f1bda2a2491ce4b42a31dbe5；上游未改。
几何/执行器/接触容量/摩擦定律/事件层/材料/空间/热学/物理配置未改。
全局/各面容限、200次Newton上限和32次回溯上限未增加。
原E3真实运行仍是FAILED、21步/.220296us、ep=0、无沉积。
旧REJECTED不能恢复；本轮局部成功不能改写旧真实运行状态。

## 实际调用及机制

SolidSolver中点内力 -> domain_kick -> kick.instantaneous/midpoint ->
cone_solve -> contact_cone.solve -> 同一凸圆盘QP -> kick.interval记冲量/原始热。
前后调用顺序不变，方法不是全局隐式，也不改变原冻结几何近似。

原QP在共同力/速度尺度fs、vs下，H=(c*fs/vs) A^T M A，q=b/vs，
约束每个二维x_a落在半径caps_a/fs的圆盘；保留全部非对角接触耦合。
476快面控制vs，慢面q/H只有约1e-14；原统一投影步长的
x-Proj(x-step*g)损失相对精度，Newton秩检测也难区分弱方向。
旧候选与80位参考力差.45866N，而真实局部区间仅1.32e-22s。
不能因为末速度影响小就忽略力方向，也不应把它诊断成宏观物理爆炸。

CONFIRMED：统一尺度对慢面不够分辨；新增独立64组矩阵原算法46组失败，
其中有些全局证书很小仍有显著唯一解力误差。新测试没有调用被测函数生成答案。

### 最终最小改动

每个圆盘使用正标量：

    rate_a = max(||H||_2, c*||M||_2*fs/vs, ||q_a||_2)
    step_a = 1/rate_a
    F_a = x_a - Proj_disk(x_a-step_a*g_a)
    J = I - D_projection*(I-diag(step)*H)

任意正标量逐盘步长有相同KKT零点；H、q、caps、M均未改变。
完整H范数保证step_a<=1/||H||；完整mobility尺度防止旋转后
理论零切向行的舍入噪声被放大。这只是求解预条件，不是物理正则化。
仍同时核对投影残差与原全局gap，内部2e-14、调用方2e-12不变。
外部solve接口和三维返回力不变；异常上下文增加disk_projected_steps。

### 初版失败没有删掉

第一版用局部H块行与||q_a||/r_a缩放：新64组和冻结476通过，
但旧432组15组失败。因此没有进行真实试跑。
原因是共同力单位fs已经归一化所有未知量，再除小圆盘半径会恶化尺度，
仅切向局部H还会放大受约束近零行的舍入噪声。
源快照contact_cone_first_candidate.py.txt、contact_class.log及class_artifacts都保留。
只读去除/r的诊断仍有5组失败，见common_force_scale_diagnostic.log。
最终采用完整H/mobility数值尺度，未逐个粒子设分支或调整物理参数。
这是同一预条件方案的回归修正，不是连续重跑AFSD寻找下一个粒子。

## 实际测试

| 入口/证据（logs/contact476_20261001/） | 结果 | 含义 |
|---|---|---|
| freeze_before.log | PASS767 | 修改前旧证据完整，非科学验收 |
| before_scaling.log / before_artifacts | 2FAIL+1ERROR；64组46失败 | 未改源码时独立反例；ERROR为预期的旧求解异常，原样保存 |
| after_scaling.log | 初版3PASS | 不能单独批准初版 |
| contact_class.log | 初版8PASS/1FAIL；432组15失败 | 被后续回归否定的初版证据 |
| final_scaling.log / final_local_artifacts | 3PASS；64/64 | 最终代码的参考力、跨尺度和完整冻结kick |
| final_contact476_original.log | 1PASS | 原失败测试文件没有改动；历史FAIL仍保留 |
| final_contact_class.log / final_class_artifacts | 9PASS；432/432 | 原多面/约束/非唯一/变换/反转/功/静动墙等 |
| gates_driver.log / gates/summary.json | 116项111PASS/5FAIL/0ERROR，exit1 | 没有新增失败；不把总退出码改成绿 |
| 新真实2us、CFL/空间细化、三工况、GPU | NOT_RUN | 本轮未运行，不继承局部PASS |

最终完整回归的5FAIL仍是legacy反转1、边界patch1及原scientific的
raw客观性/raw功/当前材料能量3；闭循环诊断PASS仅确认旧.146111J非零功，
不是认可该材料的普遍能量一致性。E2应力FAIL未改。

### 冻结476定量结果

最终力与上一轮80位独立参考之差3.789257351616192e-16N，
冻结局部2*c*M*sum(delta_f)末速度差1.3923246060282154e-31m/s。
参考数值和阈值在CONTACT476_SCALE_CONTRACT.md及原JSON先于本次修改冻结。

整个冻结kick现在覆盖1.0490073656537613e-8s，258个接受子区间，
5个内部拒绝细分。两面原始摩擦热约9.86817e-14/2.25771e-11J，
动能-冲量功残差6.18003194e-26J。没有裁热、速度归零覆盖或残差热。
发生近零滑移事件后仍按原可行性/反转逻辑处理。
完整kick的冲量/功/时间检查是实现一致性，不是独立精确轨迹答案；
1D反转的解析轨迹/热参考由原event_contact/接触类测试另外覆盖。

## 为什么没有再跑整体

gates/artifacts/preflight_timing.json：657粒子，初始dt1.0490609943676047e-8s，
2us约191步，首步实测3.5764233s，线性启动估计683.0968503s。
按本轮先验“估计>600s不启动”，新真实算例NOT_RUN。
首步含初始化/缓存因素，不是稳态步速；不能据单次数据断言修复造成固定性能退化，
也不能重测直到挑出一个低于600s的数字。本轮没有增加上限或偷偷缩短目标终时。

## 下一工作单元与止损

本地476阻塞已解除，停止继续针对这个粒子打补丁。
下一步是集成预算与实现核验：先预先定义一个<=600s的固定步数/短窗口检查，
或独立测量初始化与稳定步耗时后重新评估原2us预算；不直接宣称已跨过真实第22步。
若采用较短窗口，必须显式命名实现检查，不能拿不同终时声称收敛，
也不能改几何/质量/模量/时间尺度或把它算作沉积完成。
只在下一轮事先确定的预算和范围内启动一次；不自动授权GPU/并行/长跑/生产。

保留风险：一般Newton有限步收敛没有普遍保证、事件几何仍冻结；
旧求解循环最后一次更新后未重判接受的边界问题，本轮未扩范围修改。
材料自由能、空间边界、宏观计算成本、无限供料/真实接合仍未解决。
交接与最终哈希见HANDOFF.md及CONTACT476_SCALE_FREEZE.json。
