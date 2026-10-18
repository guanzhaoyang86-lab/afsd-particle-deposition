# 读取兼容修正及预算说明

本轮原合同的“最多两轮新接口/读取器修复”自设预算需一次明确例外，不能假装未发生。
测试前发现逐面汇总可被交换，已修并保留11PASS+1FAIL；随后实际读取分别遇到：

1. assessment.json：相对路径未正规化，入口失败，未处理真实step。
2. assessment_final.json：把上游PySPH文件不在项目freeze索引误报source drift，未处理真实step。

只读证据显示：上游integrator.py的当前SHA256与已冻结manifest完全相同；
原verify_scientific_freeze以git commit/受跟踪源码洁净度核上游，而项目files_sha256不重复包含它。
修正为：先验证run manifest本身在freeze中的hash；再逐源核当前hash，
项目文件还必须直接匹配freeze索引，上游文件允许由已冻结manifest中的hash传递核对。
没有跳过hash、放宽物理容限或改动旧源码；首次/第二次失败JSON及日志保留。

这是总计第3次本轮新读取器修正，超过原自设2次预算；允许的唯一额外动作仅为此
有证据的登记结构兼容修正及其回归，随后最后一次完整只读调用。
若仍失败，停止新读取器修复并把真实惯性评估标BLOCKED；不新增物理步或延长任何60s预算。
该说明在最后一次完整读取前登记。无需新增系统/外部授权，不影响冻结物理。
