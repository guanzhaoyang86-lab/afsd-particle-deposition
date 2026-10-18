# Predictor + Onset 独立制造参考（执行前冻结，2026-10-02）

范围：三自由度 SPD 线性弹簧的可行初猜、接触起滑定位和原子提交。
不是 AFSD、AA6061、热学或三维空间离散验证，不产生材料参数校准。
测试复用既有 LinearFixture 的实际 ConeBlock.state、cone_rows 和纯 predictor；
驱动必须是未修改的 advance_onset，参考答案由以下解析消元得到，不能调用被测函数生成。

## 物理定义与独立推导

SI：位移 m、力 N、刚度 N/m、功 J。自由度 (x,s,q) 的
K=[[2,1,-1],[1,3,-1],[-1,-1,2]] N/m。
x=-1 m 为固定法向位移，q 是待求自由变量，s 是接触材料切向位移。
刚壁位移 w(t)=2.5t m，路径时间 t∈[0,1] s；绝不直接规定释放后材料 s。
初始 q=-1/2 m、s=0、N=3/2 N、tau=-1/2 N、Es=3/4 J。
mu=10 是使制造路径具有明确有限起滑点的输入，不是 AFSD 摩擦系数。

自由平衡第三行给 q=(s-1)/2。由 K 独立消元：
N=(3-s)/2，tau=(5s-1)/2。黏着时 s=w，tau=mu*N 的首个正向交点
s*=31/15 m，q*=8/15 m，t*=62/75 s。
之后正向滑动保持 s=s*、q=q*、N=7/15 N、tau=14/3 N。
终点壁 w=5/2 m，累计相对滑移 s-w=-13/30 m。

独立功积分：

- 黏着段 Q=0；W_stick=∫[0,31/15] (5s-1)/2 ds=155/36 J。
- 滑动段材料位移为零；Q_slide=(14/3)*(13/30)=91/45 J。
- W_wall=155/36+91/45=1139/180 J。
- W_material=Delta Es=155/36 J，Es_final=91/18 J。
- W_wall-W_material-Q=0，固定法向与固定支承功均为零。

原目标初猜只替换 s=5/2 而保留 q=-1/2，会得到 N=-1 N，旧 guard 必须失败。
预测修正给 q=3/4、N=1/4、tau=23/4；这只是超摩擦圆的强制黏着分支根，
不能直接提交。全释放目标根虽满足末端摩擦圆，若从初始态直接提交则跨过
未分割的黏着/滑动事件，也必须拒绝。
对壁 w=4 m 的强制黏着分支，自由平衡 q=3/2、N=-1/2，不得裁成零或接受。

## 执行前冻结的容限和停止条件

沿用 Onset 独立回归的原事件/功门槛：force 9e-5 N、slip 1e-11 m、
work 1e-6 J、negative dissipation 1e-10 J、event time 1e-12 s、
最多40次二分/12事件、mode probe 1e-4。路径时间归一化到1s，
不增加二分次数；旧 predictor + cone 共用14修正/8次线搜索折半不变。
良态制造线性根沿用原 LinearFixture 的 relative root 1e-12。
独立解析输出核对采用绝对1e-10 m/N/J；事件时间1e-12 s。
此处解析误差容限不替换真实小块容限，也不对真实材料/AFSD授予认证。

提交前重新计算状态，核对旧 cone 可行性/自然映射残差、自由平衡、gap、
起点黏着/滑动边界和方向、signed midpoint Q、功恒等式和事件 bracket。
所有检查完成后才更新状态。禁止 abs(Q)、max(Q,0)、能量残差热、伪逆。
测试另用独立指定的相反方向功例要求 vector_contact_work 返回 -1 J 原值。

计划执行：单次独立 unittest，外 subprocess timeout=60 s，单线程。
预期小于2s；失败原样保留，不通过调整物理或容限追绿。
旧 negative-normal guard 保留 expectedFailure RED；其余未运行项目 NOT_RUN。
