# 边界第一步：只定位闭合缺口，不修改求解器

## 预注册范围与容限（2026-09-30）

本轮只新增独立诊断与回归，不改src/、configs/、已有测试、材料、接触或空间算子。
不跑AFSD或新的耦合试算，不增加GPU/依赖。旧patch失败保留。

固定物理立方体[-1,1]^3 m、rho=1kg/m³、常应力
sigma=[[2,.4,-.3],[.4,3,.2],[-.3,.2,1]] Pa；测试n=5、9，h/dx=1.3。
这是制造解测试域，不是修改AFSD几何/物性。张量积梯形体积/面积权重；
每个面的面积4m²、体积8m³。外牵引独立由sigma*n计算，不取负内部力。

虚速度（不作为粒子演化驱动），ell=1m、v*=1m/s、X=x/ell等：

- affine：v=v*(X,Y,Z)，解析功率∫sigma:grad(v)dV=48W；
- bubble：v=v*(X(1-X²)(1-Y²)(1-Z²),0,0)。六面速度均为0，
  故∫边界(sigma*n)·v dA=0W；恒应力下体积积分亦为0W。
  解析梯度用多项式直接求导，另用固定3点张量Gauss精确积分交叉核验。

将缺口分解为：

    W_sph - W_surface = (W_sph-W_analytic_gradient_quadrature)
                      + (W_analytic_gradient_quadrature-W_exact)
                      + (W_exact-W_surface)

独立记录每项W值、体内/面/棱/角残差、最强体内节点，以及
W_sph=-sum(v·F_internal)的离散伴随恒等式。不把诊断恒等式PASS称作patch通过。

预设数值门槛：体积/面面积1e-12绝对；解析Gauss参考1e-12W；
仿射/伴随/分解代数恒等式1e-10W；解析梯度独立中央差分1e-8(1/s)。
非零残差若超过1e-8W或对应力尺度1e-8才算明显离散缺口，非浮点误差。
只保留旧节点patch验收1e-10，不放宽。n=5/9仅定位，不拟合收敛阶。

若零边界速度的bubble仍有缺口，且有非零体内节点力：仅重分配到几何边界节点
的牵引无法使此固定体积梯度/力算子满足任意节点虚速度的离散分部积分。
这不排除边界迹重构/配套体积积分或ghost等方法；本轮不实现任何一种。

计划产物：probe_boundary_closure.py、独立回归、原始JSON/NPZ、旧patch失败日志。
每个入口60s上限。本构自由能变体/边界修复/新耦合基准/AFSD均NOT_RUN。

## 实测结果与本轮结论

诊断完成，4项参考/恒等式测试PASS；原节点patch重新运行仍FAIL（0.32615719），
未删除或改容限。全部src/和configs/未改。没有重跑完整116项；其本轮状态NOT_RUN，
上一轮111PASS/5FAIL只能作为历史结果。修改前METHOD_CLASS_FREEZE的521项校验PASS。

| 每轴粒子数 | 解析/边界bubble功率W | SPH功率W | 梯度重构项W | 体积求积项W | 表面求积项W |
|---:|---:|---:|---:|---:|---:|
| 5 | 0 | 1.41221751 | 2.19346751 | -0.78125000 | 0 |
| 9 | 0 | 0.61350901 | 0.82884104 | -0.21533203 | 0 |

解析梯度的独立Gauss积分为-4.16e-17W，符合精确值0。
离散体积功与负内力功的差仅约1e-16W，仿射参考48W通过。
体内节点外牵引严格为0，但内力残差范数分别为2.42482896N、1.18057376N。
节点/面/棱/角分组及最强体内节点坐标见JSON；NPZ保留完整数值数组。
两档bubble误差下降只是有限观察，不声称收敛阶，也没有继续增加测试网格。

CONFIRMED：在此固定域/常应力测试中，当前体积梯度及求积与几何边界节点牵引
不满足一般节点虚速度下的离散分部积分。缺口存在于体内，不是表面求积误差；
现有梯度/内力仍保持所测离散功共轭。不能据此把整个corrected算子判成全错。

数学上的限制很直接：若v在所有几何边界节点为0，任何只加载这些节点的
外力都有sum(F_ext dot v)=0；而本次W_internal非零。因此只调这些节点的
牵引权重，不可能让原体积算子对此场闭合。此结论不涵盖具有独立迹重构的
边界离散；若边界插值依赖体内自由度，必须重新定义其一致性和虚功，而非
把伪体力或F_ext=-F_int当成物理牵引。

下一轮仅建议：在固定5³点云上设计一个体积—边界离散分部积分的配套原型，
先用本轮零迹bubble与常应力基准验证；尚不接入AFSD，不替换现有空间算子。
若需要修改体积梯度，必须同时保留伴随内力、仿射/客观性/合力矩及cond门禁，
显式说明方法偏离。本轮原型实现NOT_RUN，本构问题也没有因这次定位而解除。

## 重跑与证据路径

项目根执行（目录须使用新名称）：

```powershell
& .envs/pysph/python.exe -B scripts/probe_boundary_closure.py --output outputs/boundary_step1_recheck_01
& .envs/pysph/python.exe -B tests/review/test_boundary_closure_step1.py
& .envs/pysph/python.exe -B tests/review/test_spatial_probes.py --boundary-red
```

本轮各命令实际通过subprocess.run(timeout=60)调用；探针约0.066s，单入口不足1s。
前两条exit0只表示诊断/参考验证完成，第三条exit1保留方法FAIL。
日志在logs/boundary_step1_20260930/；JSON/NPZ在outputs/boundary_step1_20260930/。
新增脚本/test与本文件，无生产源码改动；HANDOFF更新。
新BOUNDARY_STEP1_FREEZE.json为完整性快照，不是patch或论文验收。
