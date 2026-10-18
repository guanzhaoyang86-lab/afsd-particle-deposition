# 统一工具可见邻域：开源源码核对（2026-10-04）

本轮查到了可直接借鉴的**实际邻域过滤结构**：Peridigm 在构造最终邻域前应用几何排除标记；MOOSE 在排除跨裂纹键后对称建立邻域，并同时重建支撑统计。它们支持在当前 PySPH 候选搜索后增加统一可见图的路线，不能证明 AFSD 移动工具、热塑性、接触和图切换的联合问题已经解决。

本轮只读检索、下载少量源码作证据并检查本地固定 PySPH；没有安装后端、编译或运行外部测试，也没有修改冻结物理源或调用物理求解。下载的是审查材料，不是加入运行时的依赖。11 个源码/许可证文件的固定提交、URL、Git blob（外部文件）及 SHA256 见 `outputs/open_source_visible_graph_review_20261004/source_review_receipt.json`；相邻目录保存原文件。

## 已核对的仓库与版本

| 官方仓库 | 本轮固定提交 | 已读取许可证 | 直接检查的核心路径 |
|---|---|---|---|
| [Peridigm](https://github.com/peridigm/peridigm) | `f25c4928e850b2aa0012bc5d3b8b1ec5c78c8903` | [BSD 3 条款](https://github.com/peridigm/peridigm/blob/f25c4928e850b2aa0012bc5d3b8b1ec5c78c8903/LICENSE.md) | `src/io/pdneigh/BondFilter.cxx`、`NeighborhoodList.cxx`、`unit_test/utFinitePlaneFilter.cxx` |
| [MOOSE](https://github.com/idaholab/moose) | `e6a7269a5c8e884ea82a7b64ac687043d219022c` | [LGPL 2.1](https://github.com/idaholab/moose/blob/e6a7269a5c8e884ea82a7b64ac687043d219022c/LICENSE)，目标文件头同样声明 | `modules/peridynamics/src/mesh/PeridynamicsMesh.C` |
| [SPHinXsys](https://github.com/Xiangyu-Hu/SPHinXsys) | `846029ffef0906913685c4ba0a0911a16d98396d` | [Apache 2.0](https://github.com/Xiangyu-Hu/SPHinXsys/blob/846029ffef0906913685c4ba0a0911a16d98396d/LICENSE) | `src/shared/particle_neighborhood/neighborhood.cpp`、`src/shared/body_relations/inner_body_relation.cpp` |
| [PySPH](https://github.com/pypr/pysph) | `1cd3dea9cfaac7d9f1bda2a2491ce4b42a31dbe5`，本地 HEAD 核对一致 | [BSD 3 条款](https://github.com/pypr/pysph/blob/1cd3dea9cfaac7d9f1bda2a2491ce4b42a31dbe5/LICENSE.txt)，目录内例外另计 | `pysph/base/linked_list_nnps.pyx` |

固定提交的具体文件经 GitHub contents API 取得原始内容，不把网页摘要当成源码。这里只记录许可证事实；本轮没有向项目物理模块复制这些代码。若以后复制，应连同相应许可证和来源说明处理。

## 哪些代码真正有帮助

**Peridigm：可借鉴过滤阶段与证据布局，几何判定不能原样搬来。** `NeighborhoodList.cxx` 703–715 行先做半径搜索、排序、清空排除标记，再依次调用全部 `filterBonds`；732–737 行仅将未排除节点装入最终邻域。这是实际使用路径，不只是一个孤立几何函数。`FinitePlaneFilter` 在 `BondFilter.cxx` 123 行起按候选粒子原索引写排除标记；同文件还实现圆盘、三角形、四边形过滤。其 `utFinitePlaneFilter.cxx` 164 行起用八节点网格逐粒子检查预期排除掩码，可借鉴为“完整掩码、ID、交换端点”几何夹具。[构图调用链](https://github.com/peridigm/peridigm/blob/f25c4928e850b2aa0012bc5d3b8b1ec5c78c8903/src/io/pdneigh/NeighborhoodList.cxx#L703)、[过滤实现](https://github.com/peridigm/peridigm/blob/f25c4928e850b2aa0012bc5d3b8b1ec5c78c8903/src/io/pdneigh/BondFilter.cxx#L63)、[实际单测源码](https://github.com/peridigm/peridigm/blob/f25c4928e850b2aa0012bc5d3b8b1ec5c78c8903/src/io/pdneigh/unit_test/utFinitePlaneFilter.cxx#L164)。

边界很明确：该有限平面代码把平行或共面线段返回为无交，使用固定 `1e-14` 分母判据，并判断表面相交；它不计算穿过三维工具实体的连续内部长度，也不处理本项目圆柱减开放方孔的运动实体。因此采用它的“候选→排除掩码→最终图”结构；本项目仍应使用已独立核对的解析实体区间几何，不借平面容差放行真实遮挡。[上述几何分支](https://github.com/peridigm/peridigm/blob/f25c4928e850b2aa0012bc5d3b8b1ec5c78c8903/src/io/pdneigh/BondFilter.cxx#L63)。

**MOOSE：筛后邻域和支撑统计必须一起建立。** `createNodeHorizBasedData` 246–313 行先按 horizon 选候选，再检查预设裂纹；没有相交时同时加入 i→j 和 j→i，并同步累加各自 `_pdnode_horizon_vol` 与 `_pdnode_weight_normalizer`。本项目可借鉴“一个对称筛选结果支配所有依赖量”，不能照搬其 horizon 权函数。`checkCrackIntersectBond` 566 行起是二维裂纹中心线和宽度两侧线段检查，底层还带共享端点/角度特例；它不是可靠的本项目三维实体夹取器。该路径也是网格建立时的裂纹几何，不能据此声称已经解决持续移动工具下的构型更新。[固定源码](https://github.com/idaholab/moose/blob/e6a7269a5c8e884ea82a7b64ac687043d219022c/modules/peridynamics/src/mesh/PeridynamicsMesh.C#L246)。

**SPHinXsys：有清楚的邻域构造扩展点，但本次所读路径没有工具遮挡筛选。** `NeighborBuilderInner` 85–98 行和 adaptive 版本主要按核截止距离选邻居；`NeighborBuilderSelfContact` 123–142 行比较当前距离和初始距离，用途是区分自接触；`InnerRelation::updateConfiguration` 负责重建邻域。`Neighborhood::removeANeighbor` 会一起移动 j/W/dW/r/e 数组，说明删除邻居时缓存字段也必须保持一致；这个删除函数本身不重算本构/梯度，也不证明后处理删边有效。本次只检查列明文件，不能把它扩大为对整个仓库不存在可见性方法的断言。[邻域代码](https://github.com/Xiangyu-Hu/SPHinXsys/blob/846029ffef0906913685c4ba0a0911a16d98396d/src/shared/particle_neighborhood/neighborhood.cpp#L16)、[构型更新](https://github.com/Xiangyu-Hu/SPHinXsys/blob/846029ffef0906913685c4ba0a0911a16d98396d/src/shared/body_relations/inner_body_relation.cpp#L14)。

**PySPH：继续提供现有三维候选搜索即可。** 固定版本 `LinkedListNNPS.find_nearest_neighbors` 180–189 行计算三维平方距离，当其落在任一粒子的支撑半径内便加入邻居，接口没有第三方工具实体。由此可明确新筛选应插在本项目候选边生成与算子构造之间；不需要为了这项功能切换整个后端。[固定源码](https://github.com/pypr/pysph/blob/1cd3dea9cfaac7d9f1bda2a2491ce4b42a31dbe5/pysph/base/linked_list_nnps.pyx#L180)。

## Jacquemin 2021：找到作者仓库地址，尚未取得可审查代码

[作者论文](https://orca.cardiff.ac.uk/id/eprint/141777/1/nme.6703.pdf)第 19 页 Data Availability 明确给出 [tjacquemin/collocation-visibility](https://gitlab.com/tjacquemin/collocation-visibility)。它不是前轮未找到的 GitHub 仓库。本文研究 2D/3D GFD 线弹性凹域的可见支撑选择，包含角度阈值，能提供凹角/孔洞基准思路。

本轮 GitLab 网页读取未取得内容，直接 API 查询也未在本次有界研究内完成，已停止；所以**没有核实该仓库提交、许可证、具体函数或测试**，不将其列为已经取得的开源代码。无法访问不等于仓库已删除或私有。论文中的角度容差、包括示例采用的 5°，也不能转化为 AFSD 工具穿墙边的放行条件。

## 本项目可实施的边界

以下是结合实际源代码与本项目守恒链得到的工程设计判断，不是任何外部库已验证的 AFSD 功能。

1. 在独立变体中，由每次实际候选的 `x/pose` 和 raw NNPS 无向边产生一次可见图；保存 raw 索引、PID、active/blocked 掩码、严格实体区间、最长内部段、容差政策及图指纹。标签和 AS/RS 不参与该函数。明确区分新的离散模型阈值与旧输出拒绝阈值，做其敏感性，不能隐藏内缩。
2. 密度/修正矩阵、本构梯度、配对内力、粒子间导热及搜索切线全部读取该图。重算修正矩阵及秩/条件数；不能只删输出 i/j，不能缓存旧矩阵。导热仍采用同图上的原非负对称径向系数，不把机械修正梯度替换进去。工具反力、摩擦及壁热独立保留。
3. 本轮没有找到可直接移植并保证本问题收敛的“图切换 Broyden”实现。建议将 active 图指纹变更作为搜索矩阵刷新触发，预算内重复切换/失秩就拒绝。每次试点、线搜索和 fresh 接受仍重算完整残差；固定图近似切线不冒充跨切换的完整导数。这个建议须用切换两侧残差与守恒测试证明，不能靠仓库名称背书。
4. 首先验证无阻挡 map 等价、对称性/顺序/旋转协变、穿壁隔绝力热、过滤后仿射再现与虚功/净内力/力矩、等温与热守恒、支撑退化拒绝、固定图切线以及接触和图同时切换。通过后才能从明确合格父态登记有界新过程；旧失败候选只能用于差异诊断，不可改判。端点图可见不自动证明整步扫掠图可见。

最直接的复用结果是 **Peridigm 的过滤构图顺序 + MOOSE 的筛后双向图/支撑同步思想 + 当前 PySPH 候选搜索**。几何、热塑性、图切换和接受链仍需本项目独立实现与验证；本报告不声称已经修复原 39.875 ms 遮挡失败、获得肩外沉积或完成论文验证。
