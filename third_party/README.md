# 第三方后端

`xsets/X-Sets.cpp` 和 `uber/Uber.cpp` 源自 James S. Plank 的 X-Set / Uber-CSHR 工具；作者、来源和完整 BSD 条款保留在各源文件开头。分发二进制时也须提供相应版权和许可声明。仓库此前已经增加了数学线路输出。

- **X-Sets**：迁移中只移动路径，源码字节不变。SHA-256 为 `2a05bf43ff8909b5dabea0990f522003fb445f018261bfcfbf3a2a5ffce9afde`。`ringlab` 使用它运行 `MW_SQ` 和 `UBER_XSET`。
- **Uber**：保留独立可执行文件 `build/uber`。修复 `Schedule::Do_Uber()` 删除完成目标后又读取 `t->Min` 的问题：删除前保存完成标志，再使用标志。原实现由三行矩阵 `tests/uber-smoke.txt` 稳定触发 ASan 错误；`tests/uber_asan.sh` 检查修复。历史 JSON 线路不因此重写。

原始后端假定输入为格式正确的二进制矩阵。推荐从 `ringlab analyze` 或 `ringlab profile` 调用；适配层校验输入，去除免费输出别名，恢复输出顺序，逐门核验并剔除死门。

`UBER_XSET` 是 X-Sets 中的一种方法，与独立的 `uber` 程序不同。修改后端或改用 `--backend` 时，必须分别记录来源，不能把历史源码哈希当作任意可执行文件的认证。
