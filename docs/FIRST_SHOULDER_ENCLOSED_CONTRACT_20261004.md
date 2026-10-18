# 同一首肩候选的独立几何补证与事件筛选（运行前）

旧first_shoulder_event/approach运行接受15步至21.875ms、累计151步。
step016完整力学/能量候选触及SHOULDER，但旧连续罚壳审计256次查询后为UNRESOLVED，
因此没有发布22ms候选。旧状态FAILED_WITH_ACCEPTED_PREFIX及所有源/输出保持不变。

新的独立目录first_shoulder_enclosed从同一21.875ms接受父态重新求解，父SHA为
6d5284d82600fe89dac405d2ca26545b448397ec230d08fa43bb3a12c836699a。
不可重启的旧22ms候选SHA为904714c47961a5325e259a69b61468d9d6213e289c914ab00ab1d8144526fd94，
只作补证输入，绝不作为接受历史或新求根初态。

## 受控复用与来源

新驱动显式复用未修改的run_first_shoulder_event.py worker及接受/候选/局部对照函数。
configure仅在独立worker进程设置新BASE、新preflight、audit_contact_shell替代函数、
MAX_APPROACH_STEPS=1；另为新目录manifest增加controlled_adapter元数据的JSON写出包装。
原物理trial.evaluate、材料、摩擦、导热、实体扫掠、接受门槛、Broyden及Newton回退
均原样复用。此适配及所有替换项在每份新manifest中明确记录，旧源文件不修改。

运行前独立核验旧15条接受记录的hash/直接父链、151条累计能量前缀、材料历史、
源、参数及时间，核旧step016不可重启标志与原UNRESOLVED报告；新几何测试receipt、
同一已存父/候选的只读补证、独立曲率review必须绑定上述两个输入SHA并通过。
旧driver与9项fixture及receipt、新源/合同/证书/测试/receipt/review均冻结。

## 仅扩展保守几何证明

contact_shell_enclosure先原样执行旧审计，保留original_audit。只有原query/depth
预算耗尽的剩余平面patch区间才可补证；FOUND、实体错误、非平面剩余或无法证明仍拒绝。
固定同侧无限平面距离的端点弦减去M h²/8及舍入余量，是有限patch距离的下界，其中
M=2θ|dx-do|+θ²max(|x0-o0|,|x1-o1|)，导数针对全局步内fraction，h为子区间宽度。
原有限patch查询仍最多256次；额外解析平面路径求值单独最多256次。不增加真实接触力
或材料求值，不把密集抽样当证明。只有全部剩余区间证明后才为CERTIFIED。
该证书只覆盖声明插值的两端不活动patch对，仍不代替活动接触时间精度或真实轨迹验证。

## 有界首肩筛选

approach只重新求1×0.125ms，从21.875至22ms。若出现首次SHOULDER，数值候选及事件记录
保留但不发布；未出现也不自动扩展。macro从同一151步父态独立重解同一时间步，仅记
数值合格候选。halves从同一父态求2×0.0625ms。以21.875ms父作为共同origin，保留151条
旧能量记录，不重置历史。

各候选标量/批量完整接触和能量对照；全局位置5um、T.5K、S10%、ep20%原门槛，
以及事件粒子和一跳邻域逐点同门槛（S地板1MPa、ep地板1e-5）全部保留。
非SHOULDER新feature、并列/实体/吸引法向/负耗散/热和累计能量失败仍拒绝。
macro不发布，两half都先未提交；只有half2全局和局部对照均通过才依直接父链原子发布。

预测使用已测最大完整根成本（包括旧失败阶段剩余耗时），乘尝试步数、1.5余量及30s预留，
超过450s不启动；内部540s/外部600s硬停和每根剩余预算检查保持。单根24/100的Broyden
失败后仅允许一次同父同初值原Newton回退。通过也仅是首次肩接触工程筛选，不宣称沉积。
