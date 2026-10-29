# 第三方后端

`xsets/X-Sets.cpp` 和 `uber/Uber.cpp` 源自 James S. Plank 的 X-Set / Uber-CSHR 工具；作者、来源和完整 BSD 条款保留在各源文件开头。分发二进制时也须提供相应版权和许可声明。仓库此前已经增加了数学线路输出。

- **X-Sets**：保留作者、来源和完整 BSD 条款；当前仓库 SHA-256 为 `8fe109bc7e02b65a756c5ab1495f28a0cbc9e8fd3d0063120c3a9c772f45ef45`。`ringlab` 使用它运行 `MW_SQ` 和 `UBER_XSET`。
- **X-Sets 本地补丁记录**：在不改变候选选择规则或线路语义的前提下，增加可选的 `--quiet` 输出模式；`ringlab` 后端启用该模式以避免逐步进度输出和 `fflush`。`Bitset::SetKey()`/`XOR_With()` 改为复用预分配字符串，`Diff()`/`Ones()` 改用字节级 `popcount`，并把两个内层临时 `Bitset` 从堆分配改为栈对象。补丁只影响运行时间和诊断输出，不重写历史 JSON。
- **候选评估并行化**：`MW_SQ` 和 `MW_SMALLEST_SUM` 的独立前瞻评估可通过 `XSETS_THREADS=N` 并行执行；未设置时使用硬件并发数，设置为 `1` 可复现串行路径。worker 只写入按候选位置分配的结果槽，主线程仍按原有候选顺序评分和处理平局，因此线路、计算历史和 XOR 数保持一致。线程数只影响运行时间和内存占用。
- **固定基准**：在 `xordepth/data/matrices/vandermonde_8to17_direct.txt`、`-O3`、同一进程参数 `0 3 0` 下，优化前后均逐门核验相同线路计数：`MW_SQ` 为 219，`UBER_XSET` 为 212。macOS 本机墙钟时间从约 27.9 s 降至 8.1 s、从约 4.5 s 降至 1.2 s；时间用于确认运行优化，不构成跨机器性能声明。
- **Uber**：保留独立可执行文件 `build/uber`。修复 `Schedule::Do_Uber()` 删除完成目标后又读取 `t->Min` 的问题：删除前保存完成标志，再使用标志。原实现由三行矩阵 `tests/uber-smoke.txt` 稳定触发 ASan 错误；`tests/uber_asan.sh` 检查修复。历史 JSON 线路不因此重写。

原始后端假定输入为格式正确的二进制矩阵。推荐从 `ringlab analyze` 或 `ringlab profile` 调用；适配层校验输入，去除免费输出别名，恢复输出顺序，逐门核验并剔除死门。

`UBER_XSET` 是 X-Sets 中的一种方法，与独立的 `uber` 程序不同。修改后端或改用 `--backend` 时，必须分别记录来源，不能把历史源码哈希当作任意可执行文件的认证。
