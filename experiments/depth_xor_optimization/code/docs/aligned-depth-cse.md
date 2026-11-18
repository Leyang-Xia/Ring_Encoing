# 仓库原流程复现与 OUR_DCCSE 接入

## 范围与来源

基准版本为 `Leyang-Xia/Ring_Encoing` 的 `4b6a4e83b00c2cd66cfd9798e41959ca0361a25a`，实验开始时也是远端 main。此次严格复现的是 `xordepth/docs/count-depth-frontier.md` 的四组正式矩阵和审计，以及 `depth-bounded-search.md` 的两组深度 6 搜索；另运行根 README 的 RS direct UBER_XSET profile。不是重新执行全部家族、环参数及专利研究搜索。

- 输入逐字节使用 `xordepth/data/matrices/` 原文件，没有生成随机矩阵。
- `frontier` 使用原 `--share 2,3,0`，Vandermonde/RS 使用原 `--also` 证据文件；完整输入到完整输出、不计接线与免费别名的成本口径不变。
- 新算法直接进入共享核心 `src/depth_cse.cpp`；根 `ringlab analyze --method OUR_DCCSE` 和 `xordepth frontier --include-ours` 调用同一实现。
- 不加 `--include-ours` 时，四组输出文件与未修改版本逐字节一致。原 60 条候选保留，新增 8 条候选与原候选共用 `candidate()`、Pareto 标记和 `verify_result()`。验证额外检查有声明时的深度上限。
- 原有证据 JSON、第三方后端和矩阵均未改写。

## 参数

主比较新增候选为每矩阵 `D0` 和 `D0+1` 两档，`restarts=1, seed=0, lookahead=false, refine=true`。这与原 frontier 单次 seed0 的调用方式对齐，不将新增多重启最好结果伪装为同一基准。可选 `--ours-lookahead` 可启用前瞻，本次主比较没有启用谱方法。

两条补充 `ringlab analyze` 使用 `depth=6, seeds=8, seed=1`，匹配 Beam 的深度、种子标签和重启次数。二者的搜索空间/每次工作量不同，不能据此声称等计算预算。Beam 完整保留 `beam=64, branch=16, representations=8, restarts=8, seed=1, timeout=180`。

最小深度表的原基准是仓库原 frontier 中所有已验证候选在该深度约束下的最好值，包括导入、组合和重平衡的完整线路；不是对该深度下所有可能算法的最优性证明。历史保存候选经历的搜索预算不统一，因此比较是已有可实现上界的更新，而不是同运行时间排名。

## 运行

在指定 commit 上应用补丁后：

```bash
make -j3 test
(cd xordepth && make -j3 test)
bash docs/run-aligned-benchmark.sh /tmp/aligned-results
```

脚本默认也运行官方 180 秒 Beam 搜索及 60 秒 profile。若只重跑快速 frontier 部分，可显式设置 `RUN_BOUNDED=0 RUN_PROFILE=0`；此次交付中的两组 Beam 和 profile 已真实执行，没有用归档结果代替执行。

单独分析自己的矩阵：

```bash
./build/ringlab analyze matrix.txt --method OUR_DCCSE --depth 6 --seeds 8 --seed 1 --output /tmp/ours.json
./build/ringlab verify --input /tmp/ours.json --samples 0
```

## 实现边界

OUR_DCCSE 是已讨论的深度可行公共子表达式提取，带相同 GF(2) mask 的等价信号复用、按到达深度合并剩余项，以及受限局部门重代入。不是完整 BP/DCLO、SAT 或全局最优算法；门级深度不包含扇出和布线时延。实现当前支持 0..60 的深度上限，系数 mask 使用原动态位向量，不受 64 位限制。

共享正确性与功能验证采用原 `Circuit`；深度由核心和原 `xd::depth_profile` 两个路径交叉检查。新测试覆盖 200 输入 mask、零/单位/重复输出、固定种子重现、不可行深度拒绝和两种评分/细化开关。根 Make 和 xordepth Make 的全部测试通过；CMake 源列表/测试目标同步更新，当前运行环境缺少 cmake/ctest，未执行 CMake 构建。

本次没有推送或修改远端仓库。完整补丁和原生 CLI 结果可独立审阅。

补充检查：新增核心单元用 AddressSanitizer/UndefinedBehaviorSanitizer 运行通过；当前 ptrace 环境不支持 LeakSanitizer，因此该次检查使用 `ASAN_OPTIONS=detect_leaks=0`，不声称完成泄漏检测。
