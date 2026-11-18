# KP4 300×200：深度上限 8–13

输入是 `xordepth/data/matrices/generated/kp4_300x200.txt`，200 输入、300 输出；SHA-256 `40761ee587f563f08f27983afb43e3e5cbea6f5d992ab7989e8840c0240f1ca5`。行权重 20–108，均值 82.62；扇入为 2 时，行权重给出深度下界 7，不说明深度 7 能否在 8452 门内达到。代码版本 `85fe1a6232894d6e8f7c123ff87e862b0a2a1287`，使用 Release 构建。

## 主扫参

六档并行运行，每档单线程、墙钟上限 240 秒。JSON 保留完整 target masks 和固定输入到固定输出的活门电路；逐档运行 `xordepth verify`。

| 深度上限 | 活 XOR | 实际深度 | 最佳来源 | Beam 状态 | rewrite 尝试 |
|---:|---:|---:|---|---:|---:|
| 8 | 8452 | 8 | CSE baseline | 0 | 400 |
| 9 | 8452 | 8 | CSE baseline | 0 | 400 |
| 10 | 8452 | 8 | CSE baseline | 0 | 400 |
| 11 | 8452 | 8 | CSE baseline | 0 | 400 |
| 12 | 8452 | 8 | CSE baseline | 0 | 400 |
| 13 | 8452 | 8 | CSE baseline | 0 | 400 |

六档输出电路完全相同。300 个输出在深度 5/6/7/8 分别有 9/10/139/142 个，最大扇出 15。`verified` 表示返回的 incumbent 已验证，不表示搜索完成：主扫参耗尽墙钟预算，只发布了 CSE、重结合及 rewrite 候选，没有发布 greedy/Beam 候选。不能据此推断放宽深度无益。

## 对照与优化方向

- `beam-depth-{8,13}.json`：同时关闭 greedy 和 rewrite，每端 120 秒、最多 256 状态。两端都展开 153 状态、搜索到第 39 轮，仍只返回 8452/8 的 CSE。短时采样中，展开主要耗在 `SearchState::apply` 调用的 `rebuild_pair_index()`，每个子状态都重扫所有保留表示中的 pair。优先增量更新变化的表示并复用 pair 的 mask，而不是先加大 Beam。
- `plain-depth-{8,13}.json`：关闭 MW_SQ 式重评分与 rewrite，8 个 greedy restarts，180 秒。没有发布完整 greedy 或 Beam 候选。`greedy1-depth-8.json` 把 greedy 和 Beam 分别缩至 1 次 restart 和 1 状态，221.77 秒完成并发布了 greedy 候选，但最好仍为 CSE。普通 greedy 的短时采样也看到了候选生成、排序的开销；可先按改变的输出维护 pair 支持度，再考虑重评分优化。
- `cse-control-depth-{8,13}.json`：不做 greedy/rewrite，Beam 仅 1 状态，无墙钟限制，复现 8452/8。`cse-frontier.json` 的链式 CSE 为 8498/27，rebalance 后为 8498/8；depth-aware `share=2/3/4` 分别为 8452/9562/10255 门，实际深度均为 8。提高共享阈值减少了深度 8 的输出数，但代价是门数上升，最大深度未变。
- `ringlab-cse-{8,32}seeds.json`：独立的多种子 CSE 计数基线，两档最优均为 8486/30，audit 测得 rebalance 后为 8486/8，仍不如 8452/8。
- 当前优先把预算分给可交付的候选：先保存可验证 baseline，为 greedy 和 Beam 分别留预算，最后才 rewrite；结果增加分阶段耗时及停止原因。质量方向是对 Beam 保留的部分状态尝试可行收尾，或小规模试 residual/bridge 以改变共享选择，而不只是把同一共享电路重新平衡。上述是后续实验方向，不是已获得的门数收益。

## MW_SQ 基线

用户提供的 `MW_SQ 0 0 0` 门数为 8241，比 8452 少 211 门（以 8241 为基准，xordepth 多 2.56%）。`mw-sq-000.json` 是 `ringlab analyze` 直接路径的独立复跑：240 秒超时，无电路，不能据此确认 8241 或它的深度。更高参数可能改变可枚举的表示和贪心走向，但不保证门数单调降低。已有原始 schedule 时，优先跑 `xordepth rebalance RAW.txt --matrix "$matrix" --output MW.json`；仅当实测深度满足该档上限，才能将它作为 `--incumbent MW.json` 参与同约束比较。

另外用原矩阵行序直接运行 `XSETS_THREADS=4 build/xsets 0 0 0 MW_SQ`，限时 600 秒。超时时只有 399 条进度行，没有 `Schedule:`，所以没有原始或 rebalance 深度。记录见 `mw-sq-000-raw-timeout.json`；未启动更高参数。

## 复现

从仓库根目录运行：

```sh
matrix=xordepth/data/matrices/generated/kp4_300x200.txt
out=$(mktemp -d)
for depth in 8 9 10 11 12 13; do
  (xordepth/build/xordepth bounded "$matrix" --depth "$depth" \
    --beam 4 --branch 4 --representations 4 --restarts 2 \
    --greedy-restarts 2 --greedy-mwsq --rewrite-trials 400 \
    --max-states 256 --threads 1 --seed 1 --timeout 240 \
    --output "$out/depth-$depth.json" &&
  xordepth/build/xordepth verify "$out/depth-$depth.json") &
done
wait
```

有墙钟预算的运行，状态数随机器和并发负载变化。`make test` 和 `make -C xordepth test` 均通过。
