# 1961 DOF 方向求解成本及有界改进（2026-10-05）

本次只读成本审计、小代数资格、新变体、真实 A/B 及独立保存数据复核已经完成。
不把求解成本改进当作接触起始时间分辨率失败的修复。

## 实际源码和已保存候选

`scripts/secant_equilibrium.py::_full_matrix` 在初始/刷新切线时建立完整自由矩阵 A，
同时计算 H=inv(A)；`_broyden` 后续同时维护 A 和 H，并检查条件估计、分母及两探针漂移。
冻结 `scripts/visibility_equilibrium.py` 的每次方向却仍执行 `np.linalg.solve(A,rhs)`，
因此已有 H 未用于减少逐方向的重新稠密分解。
首个割线更新可把原稀疏结构填满，直接改用 sparse solve 并无已验证的稀疏性依据。

实际文件 `outputs/expanded_hot_budget_prefix_20261005/window_001/fine/step_014/broyden/iterations.json`：

- 918 粒子 / 1961 自由位移 DOF；18 次方向、38 次完整 map evaluate、3 次切线构造。
- 图 fingerprint 全段不变；两次刷新原因是 `HEAVY_DAMPING`，发生在迭代 5、11。
- 初始 A 非零项 788873，一次秩一更新后为完整 3845521 项；一份 float64 A 或 H 约 29.34 MiB。
- 完整候选耗时 15.6624274 s；日志只保存 19 个接受迭代的完整 map 时长，合计 3.5439231 s。
  另外 19 次 evaluate 的独立时长未完整保存，不能据此精确算出线性分解占比。
- 最大条件估计约 45056.99；最大已记录线性余项 3.3621e−16 N。
  最终力余项 4.3761e−6 N，原门槛 1.01e−5 N。

该候选所属时间筛选失败，未提交；其求根合格不改变失败结论。

## 新变体的实现边界

新增 `verified_inverse_direction.py` 与 `visibility_inverse_equilibrium.py`；旧文件未改。
先 H rhs，再用当前权威 A 实算余项；只在原 1e−9 max(1,‖rhs‖∞) 门槛内返回。
一次校正只有确实降低实算余项才保留；否则或仍不达标，使用原 dense solve，并再次检查。
A/H/rhs 不写入、不裁剪、不正则化；旧条件估计、图切换刷新、非线性残差、Armijo、
24/100/两刷新、所有 fresh/材料/热/接触/几何/能量条件不变。
截止时间失败不会被当作可再花一次刷新预算的普通线性失败。

这是利用已维护逆矩阵的搜索方法，并非新的物理解，也不是把 A 的精度换成 H 的断言。
逐 rhs 校验很重要：两个固定探针可能看不到某个逆误差方向，不能代替当前 rhs 的真实余项。
条件数仍是原 ‖A‖∞‖H‖∞ 估计；若 H 陈旧，它不等同精确条件数证明。
图/free 排序/矩阵版本由同一求解器作用域绑定，本轮没有跨时间步缓存。

官方方法参照是已有因子或逆更新的重复应用，而非直接换一个默认求解器：
[SciPy lu_factor](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.lu_factor.html)
将带主元 LU 分解与
[lu_solve](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.lu_solve.html)
应用分开；当前 A 每次割线更新会变化，不能不更新因子就套用前次 LU。
[SciPy broyden1](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.broyden1.html)
说明 inverse Broyden 更新；
[PETSc SNESQN](https://petsc.org/release/manualpages/SNES/SNESQN/)
提供多种拟 Newton 逆近似与重启策略。
本变体没有安装或调用这些外部后端，仍使用原 PySPH 设施与 NumPy；上述默认容差不是项目资格。

## 现有验证及真实 A/B 边界

`outputs/verified_inverse_direction_20261005/algebra_v1/receipt.json`：15 项 PASS、exit 0。
包括 24 次原 Broyden 代数更新、故意陈旧 H、两探针漏检的误差、一次有益/无益校正、
病态及奇异拒绝、输入不修改、期限、静态源码差分，以及 3 个制造残差/图切换例子。
没有真实 AFSD evaluate 或 NNPS。192×192 小矩阵的中位方向时间
0.0001537 s（dense+余项）与 0.0000653 s（H 路由含条件和余项）仅是代数测速，
不包含初始逆构造，不代表 1961 DOF 或过程加速。

`INVERSE_DIRECTION_AB_CONTRACT_20261005.md` 固定的同真父一次 A/B 为：
62.5 μs / 14 步，dt 7.8125 μs，各根独立重载，先同一真实初始切线方向比较，再完整原/新根。
原全部物理资格不变；算法比较另用逐点应力 0.1%、位置 1 nm、ep 1e−7 等严格预登记筛选。
完整粒子/材料/接触/热/图字段实际差要记录，不能先称机器位等价。
已测最大完整根 28.0122514 s，A/B 预测 117.0367542 s；内层 540 s / 外部 600 s。
`outputs/verified_inverse_direction_20261005/ab_preflight_v1/receipt.json` 证明语法和真父预检 PASS，
2.044 s，无物理/NNPS、无过程时间推进；它不证明 A/B 已通过。

## 真实 A/B 与独立读回结果

真实输出 `outputs/inverse_direction_ab_20261005/v1`：exit 0，两个完整候选均通过原资格及事前数值筛选。
`scripts/readback_inverse_direction_ab.py` 随后只读取保存矩阵、数组、JSON 和哈希，
`outputs/inverse_direction_ab_readback_20261005/v1/summary.json` 为 PASS / exit 0 / 2.507 s，
没有调用物理求解、tangent、NNPS 或几何重证。

| 量 | 原 dense 方向 | verified H 方向 |
|---|---:|---:|
| 完整候选成本 | 13.3346941 s | 12.0790631 s |
| 搜索方向 / 完整 evaluate | 17 / 29 | 17 / 29 |
| 切线构造 / refresh | 2 / 1 | 2 / 1 |
| 同初始实际 A d − rhs 最大余项 | 5.82867e−16 N | 1.33227e−15 N |
| 最终 fresh 力余项 | 9.800232055875e−6 N | 9.800232055910e−6 N |

两根均无 Newton 回退；新根 17 个方向全部直接 H 应用，无校正或 dense 回退。
独立读回同一保存的 1961×1961 A/H/rhs，两条方向都通过 1e−9 N 线性门槛；
方向最大差 1.86116e−22 m。候选位置最大差 5.42101e−20 m、逐点应力相对差
5.33062e−15、ep 最大绝对差 1.28749e−19，温度/密度/接触历史差为 0。
图、4 个接触布尔/整数分支、strict crossing IDs 一致；每支完整 14 条旧记录和第 15 条候选能量通过。
所有输出保持 NOT_ACCEPTED / NOT_RESTARTABLE，正式过程时间没有因此增加。

单次完整候选减少 1.2556310 s，即 9.4163%，原/新耗时比 1.10395。
这次改动只减少方向的重复分解；29 次完整 map、2 次切线及初始 inverse、Broyden 矩阵更新、
fresh/scalar 检查、接触/路径资格与输出仍执行。新方向自身也需要逐次扫描矩阵、估计条件及实际余项，
并非只有一次 H 向量乘法：17 个 direction 总计 0.6035656 s，其中 H 应用与余项部分合计 0.0474204 s。
独立初始系统的 dense+余项为 0.1257764 s，完整 verified H 路由为 0.0418039 s；
该一次测量不能当作其余 16 次方向的精确成本。

旧方向没有逐方向 wall 计时，29 次 map 也仅有 18 个迭代端点的独立时长；
所以不能准确分摊全 map、线性代数、切线、I/O 的份额，不能声称已证明某一项占主导。
这一个固定顺序 A/B 未测运行波动，也没有检验跨步成本、不同接触图或未来时间细化。
后续每轮 Broyden A/H 未保存，独立 reader 对后续线性守卫是日志检查；
只有共同初始系统能够独立重算 A d − rhs。完整候选字段和全部能量历史已从保存数组独立复核。
