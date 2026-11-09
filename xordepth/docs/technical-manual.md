# xordepth 技术手册

深度硬约束 XOR 电路合成的启发式搜索——功能、流程、组件与数学原理

本手册面向两类读者：需要**使用** xordepth 的工程师（第 1–4、7 章），以及希望**理解其数学机制**的研究与调参人员（第 5–6 章）。全部机制描述均以当前源码实现为准（`src/xd_search.cpp`、`src/xd_greedy.cpp`、`src/xd_depth.cpp`、`src/xd_rewrite.cpp`），并在关键处标注对应的函数名，便于对照阅读。

---

## 目录

1. [功能定位与适用场景](#1-功能定位与适用场景)
2. [符号约定](#2-符号约定)
3. [整体执行流程](#3-整体执行流程)
4. [核心组件及其作用](#4-核心组件及其作用)
5. [启发式搜索的数学原理](#5-启发式搜索的数学原理)
6. [完整示例推导](#6-完整示例推导)
7. [使用参考](#7-使用参考)
8. [正确性边界与已知限制](#8-正确性边界与已知限制)

---

## 1. 功能定位与适用场景

### 1.1 问题形式化

xordepth 解决的是**深度受限的 GF(2) 线性变换电路合成**问题：

> 给定一个 GF(2) 上的 m×n 矩阵 M（m 行输出、n 列输入），综合一个由**二输入 XOR 门**组成的有向无环网络 N，使得对任意 x ∈ F2^n 有 N(x) = Mx；并且网络中**所有信号的最大拓扑深度不超过给定上限 D**（硬约束）。在此前提下**最小化活的 XOR 门数**。

矩阵以文本文件输入：每行一个输出、每列一个输入，取值 0/1。输出 `cost` 只统计活的 XOR 门（补零、接线、输出别名不计门），`depth` 按二输入门的拓扑深度计算（输入与常量零深度为 0）。

问题的本质困难：门数最小化（即使不约束深度）在一般情形是 NP难的，深度约束进一步把可行空间切成非凸区域。因此 xordepth 是一个**启发式上界产生器**，不宣称全局最优。

### 1.2 适用场景

| 场景 | 例子 |
|---|---|
| 纠删编码 / RS 类编码线性层 | Vandermonde 8→17、systematic RS [7,4]、Cauchy 7→4、GF(2⁸) 展开矩阵 |
| 密码学线性扩散层 | 任意 GF(2) 线性映射（MDS 分支数矩阵等） |
| 通用 GF(2) 线性映射 | 随机/结构化稠密、稀疏、分块、模加矩阵（测试套件 10×10 至 100×100） |
| 时序关键路径受限实现 | 深度是硬约束而非事后指标——搜索中每一步都保证可满足 |

**不适用 / 不建模**：LUT 合并、多输入 XOR 门、布线与扇出时延、非 GF(2) 运算。`search_failed` 与 `timeout` 不构成"不存在深度 ≤ D 的线路"的证明。

### 1.3 方法谱系中的位置

xordepth 的候选生成借鉴 X-Sets / Uber-XSet（Plank, DSN'12）与 MW_SQ 的"共享 pair 优先"思想，但以**深度硬约束**为第一原则：

- X-Sets 系方法把"输出行的 pair 出现频率"作为共享节点的选择依据，门数表现好，但**不保证深度**（事后靠重结合补救）；
- xordepth 把每个候选的接受条件改为"**替换后所有受影响输出仍可在 D 层内完成**"（Kraft 预算，见 §5.2），从而搜索的每一步都落在可行域内；
- 最终通过统一候选池把 CSE 基线、greedy、Beam、重写、外部 incumbent 的结果按 (cost, depth) 字典序合并，形成 anytime 式的逼近。

---

## 2. 符号约定

| 符号 | 含义 |
|---|---|
| n, m | 输入数、输出数（目标行数） |
| D | 深度硬上限（`--depth`） |
| M | m×n 的 GF(2) 目标矩阵 |
| r | 一个输出行（目标行），对应位掩码 mask(r) |
| w(r) | 行权重：mask(r) 中 1 的个数 |
| v, a, b | 已物化节点（输入或中间门），`mask(v)` 是它计算的 GF(2) 向量 |
| d(v) | 节点深度：`d(输入) = 0`；`d(a⊕b) = 1 + max(d(a), d(b))` |
| R | 一个**表示**（representation）：当前用若干节点之和表达某目标的集合，`R = {t₁,…,t_k}`，要求 mask(R) = 目标掩码 |
| \|R\| | 表示中的项数（terms 数） |
| Φ(R) | 表示的完成预算占用：`Φ(R) = Σ_{t∈R} 2^{d(t)}`（§5.2） |
| edge(a,b) | 边掩码 `mask(a) ⊕ mask(b)`，即门 a⊕b 计算出的值 |
| supp(e) | 支持数：边掩码 e 在全部表示的 pair 中出现的总次数 |
| LB | 门数下界（分支定界用，§5.6.2） |
| ⟨·⟩ | 排序用的字典序元组 |

深度合并规则全文通用：**合并 a、b 得到的新节点深度为 1 + max(d(a), d(b))**。

---

## 3. 整体执行流程

`xordepth bounded MATRIX.txt --depth D …` 调用 `synthesize_bounded()`，总体是一条**从廉价到昂贵、逐级发布候选**的流水线。所有阶段的结果进入同一个统一候选池（§4.5），随时可能成为当前最优。

```text
输入矩阵 M, 深度上限 D, 预算参数
        │
        ▼
[预检] 行权重深度下界 LB_d = max_r ⌈log₂ w(r)⌉
        若 D < LB_d → 直接返回 infeasible_by_lower_bound（数学上不可行）
        │
        ▼
[阶段 0] 外部 incumbent（--incumbent FILE，可重复）
        逐条验证功能与深度后零成本入池；重复/错值/超深跳过并计数
        │
        ▼
[阶段 1] 廉价基线（在昂贵搜索之前，保证墙钟截止也不会丢基线）
        CSE 基线（synthesize_schedule，share 阈值 2，depth-aware）
        重结合基线（rebalanced：单消费者锥按到达时间 Huffman 重排，门数不变深度可能降）
        │
        ▼
[阶段 2] 可完成 greedy portfolio（synthesize_completable，--greedy-restarts N）
        N 次尝试，普通/带 lookahead 交替；每次成功的尝试经 observer 立即入池
        可选 MW_SQ 式并列候选重评分（--greedy-mwsq）
        │
        ▼
[阶段 3] 深度约束 Beam 搜索（XOR-Set 状态空间，§5.5–§5.6）
        restarts 轮；每轮受 round_limit 限制；全局受 max_states 限制
        完成状态导出电路入池（"beam" 来源）
        │
        ▼
[重写环] 池中出现新的最优时（允许重写且预算未耗尽）
        rewrite_depth_gates：尝试删除多消费者共享门（§5.8）
        结果再入池（"rewrite" 来源，不递归触发重写）
        │
        ▼
[输出] 池内最优线路按 (cost, depth) 字典序胜出
        逐门符号验证 + 深度复检 → status = verified
        预算耗尽而池为空 → search_failed / timeout
```

### 3.1 Anytime 语义

候选池的设计使整个流水线具备 **anytime 性质**：

- 墙钟（`--timeout`）、状态数（`--max-states`）、轮次上界任何一处触发截止，都**不会丢弃已经入池的已验证线路**（代码中以异常 `SearchDeadlineExceeded` 贯穿，捕获点保留 incumbent）；
- 便宜来源（incumbent、CSE 基线）先发布，昂贵来源（greedy、Beam）后发布，所以加时间总是单调不劣；
- `merge` 子命令支持把多次运行的结果离线合并（零搜索时间），实现跨预算、跨参数的最优汇总。

### 3.2 有限终止（无时间限制也可停）

`--no-timeout` 只关闭墙钟，不关闭结构性预算：

- 每 restart 的轮次上界：`round_limit = m · (2^(D+1) − 1)`（每目标至多物化 2^(D+1)−1 个树节点，§5.2 的对偶）；
- 跨 restart 的全局状态数上界：`--max-states`；
- greedy 主循环每轮**严格减少表达式总项数**（§5.3），迭代次数先天有限；
- 重写预算 `--rewrite-trials` 有限且跨发布共享。

因此无墙钟运行也必然终止。预算耗尽**不是**不可行证明。

### 3.3 结果状态机

| status | 含义 | cost 字段 |
|---|---|---|
| `verified` | 找到完整线路，逐门功能、门数、深度全部复验通过 | 有效 |
| `infeasible_by_lower_bound` | D < 行权重下界，数学上不存在解 | null |
| `search_failed` | 启发式把搜索空间耗尽仍未找到 | null |
| `timeout` | 墙钟/状态预算耗尽且池为空 | null |

---

## 4. 核心组件及其作用

### 4.1 数据底座：BitMask 与 Circuit

- `BitMask`（`include/xd/core.hpp`）：精确多字位掩码（≤128 位内联小缓冲，零堆分配），是所有 GF(2) 向量、"这个节点算什么值"的唯一载体。等值即等价：两个节点 mask 相同 ⇔ 计算同一变换。
- `Circuit`：`ops`（`[node, left, right]` 三元组列表）、`outputs`、`inputs`、`zero_node`。`cost()` 统计活门，`verify(targets)` 做符号功能验证，`pruned()` 剔除不可达门。

### 4.2 表示层：Expression 与 ParetoFrontier

搜索维护的不是"一个输出一个表达式"，而是**每个目标一个 Pareto 前沿**（`ParetoFrontier`），前沿中的每个**表示** R 是该目标的一种当前分解。表示的坐标是：

```text
( |R|, arrival_depth(R) )     arrival_depth = combine_depth({d(t) : t∈R})
```

`combine_depth` 是 Huffman 式最早就绪合并的深度（§5.2.4）：把 R 中所有项两两合并所需的最小最终深度。

Pareto 支配关系（`dominates()`）：R₁ 支配 R₂ ⟺ `|R₁| ≤ |R₂|` 且 `arrival(R₁) ≤ arrival(R₂)` 且至少一项严格小。前沿保留**项数多但深度低**的备选——它们是后续深度受限合并的"缓冲库存"。每目标的前沿容量由 `--representations` 限制，溢出时按 `future_support` 降序保留。

### 4.3 搜索状态：SearchState

一个 Beam 状态是完整的搜索快照：

| 成员 | 作用 |
|---|---|
| `nodes_` | 已物化节点表（mask、深度、两个生产者） |
| `node_by_mask_` | mask → 节点索引，等值复用查询 |
| `frontiers_` | 每目标的 Pareto 前沿 |
| `pair_index_` | **边掩码 → 出现列表**：所有表示中每对项 (tᵢ, tⲻ) 记录 `(target, rep, tᵢ, tⲻ)`，按 `edge(tᵢ,tⲻ)` 聚合。这是候选生成与 support 统计的索引 |
| `residual_pair_support_` | 跨表示 residual 候选的支持计数（`--residual-pairs` 开启时） |

初始状态：节点 = n 个输入，每个目标一个表示 R = {该行所有置位输入}，若 `arrival_depth > D` 直接标 invalid（深度下界预检的另一种体现）。

### 4.4 转移算子

三种状态转移（每个都返回子状态、失败则 invalid）：

- `apply(pair)`：物化门 `a⊕b`（若 mask 已存在则复用旧节点）。对所有包含 {a,b} 的表示执行 `R' = (R \ {a,b}) ∪ {a⊕b}`（若 a⊕b 已在 R 中则相消：`R \ {a,b,a⊕b}`），重算 Pareto 前沿与 pair 索引。
- `apply_toggle(a, b, target, rep)`：residual 模式。对单个表示执行对称差替换 `R' = normalize(R ⊖ {a, b, a⊕b})`——三项中在 R 里的被移除、不在的被加入。暴露跨表示共享节点。
- `apply_bridge(a, b, ta, tb)`：跨两个目标/表示同步执行 toggle，用于两个操作数不在同一表示中出现、但组合可行的门（MW_SQ trace 回放证实普通 apply 表达不了的门，bridge 可以表达）。

每个转移都伴随深度检查与前沿重建；`valid_ = false` 的子状态被丢弃并计入 `pruned_depth`。

### 4.5 统一候选池：CandidatePool

所有来源的**唯一入池口**是 `publish_candidate()`：

```text
offer(c)：
  1. c.verify(targets)     —— 符号功能验证（错值拒绝）
  2. depth_profile(c) ≤ D  —— 硬深度复检（超深拒绝）
  3. metric = (cost, depth)，字典序与当前最优比较
     更优 → 替换 best，记录 best_source
  4. 新最优且允许重写且预算有余 → rewrite_depth_gates 派生新线路再入池
```

来源标签（`best_source` 字段）：`incumbent` / `cse_baseline` / `rebalanced_baseline` / `greedy` / `beam` / `rewrite`。`load_incumbent` 额外做 JSON 级去重。

### 4.6 后处理

- `rewrite_depth_gates`（§5.8）：尝试删除"被多个消费者使用的非输出共享门"，一次性重写全部消费者，只在严格降门数时接受。
- `rebalanced`：单消费者锥按到达时间重排（Huffman），门数严格不变、深度可能下降。
- `merge` 命令：离线把多份结果 JSON 合并为逐门验证后的最优。

---

## 5. 启发式搜索的数学原理

本章是手册的核心。按"可行性数学 → 选择函数 → 评估函数 → 剪枝 → 逼近机制"的顺序展开。

### 5.1 深度下界：行权重的对数下界

`depth_lower_bound(targets)`（`xd_depth.cpp`）：

```text
LB_d = max over rows r of  ⌈log₂ w(r)⌉
```

**推导**：一个深度为 d 的二输入 XOR 树，其输出至多覆盖 2^d 个不同的输入（每个门至多把覆盖集翻倍；归纳：叶子覆盖 1 = 2^0 个）。要覆盖 w(r) 个输入，必须 2^d ≥ w(r)，即 d ≥ ⌈log₂ w(r)⌉。该界对每个输出独立成立（深度是网络级的 max），故对行取最大。

**用途**：这是唯一的数学不可行判据——`D < LB_d` 时直接返回 `infeasible_by_lower_bound`，无需搜索。注意它忽略共享与扇出，是很松的界：Vandermonde 8→17 的 LB_d = 5，而实际已知最优深度解恰好是 5；但对多数矩阵真实最优深度会大于 LB_d。

### 5.2 完成可行性的 Kraft 预算（核心不变量）

这是整个算法最重要的数学构件。**问题**：表示 R = {t₁,…,t_k}（各项深度 dᵢ 已知）还能不能在 D 层内合并成单个输出？

#### 5.2.1 合并规则与势函数

合并规则 `d(a⊕b) = 1 + max(d(a), d(b))`。定义**势函数**：

```text
Φ(R) = Σ_{i=1..k} 2^{d(tᵢ)}
```

对 R 中任意两项 a, b（不妨 d(a) ≤ d(b)）做一次合并：

```text
ΔΦ = 2^{1+max(d(a),d(b))} − 2^{d(a)} − 2^{d(b)}
   = 2^{d(b)+1} − 2^{d(a)} − 2^{d(b)}
   = 2^{d(b)} − 2^{d(a)}
   ≥ 0                                    （当且仅当 d(a) = d(b) 时取等）
```

**Φ 在任何合并序列下单调不减**。终止状态是单项 {t}，深度 d_final ≤ D，势为 2^{d_final} ≤ 2^D。于是得到必要条件：

```text
Φ(R) ≤ 2^D            —— "Kraft 预算"（budget = 2^D）
```

（代码：`row_weight()`；`feasible_update()` 做的是它的增量算术。）

#### 5.2.2 Kraft 对偶解释（为什么恰好充分）

把合并过程倒过来看：最终输出是根，位于深度 ≤ D；项 tᵢ 从深度 dᵢ 处**接入**这棵树，从根向下走到 tᵢ 还剩 `ℓᵢ = D − dᵢ` 层。把 tᵢ 视作这棵"完成树"中深度为 ℓᵢ 的叶子，则整棵树是一棵高度 D 的二叉前缀树。**二叉前缀树存在 ⟺ Kraft 不等式成立**：

```text
Σᵢ 2^{−ℓᵢ} ≤ 1   ⟺   Σᵢ 2^{−(D−dᵢ)} ≤ 1   ⟺   Σᵢ 2^{dᵢ} ≤ 2^D
```

即 Φ(R) ≤ 2^D 不只是必要条件，**恰好刻画可完成性**（在 1+max 合并规则下）。直观：深度 dᵢ 的项在高度 D 的树里要占据 2^{dᵢ} 个叶子槽位（例：D=3 的树有 8 个叶子槽；一个 d=1 的项自身已消耗掉一半的槽位——因为它下方还挂着自己的生产子树）。

**特例自检**：k 个深度全为 0 的项，Φ = k ≤ 2^D ⟺ k ≤ 2^D，退化为"高度 D 的二叉树至多 2^D 个叶子"这一显然事实。

#### 5.2.3 增量式判定与预算增量

搜索中每次替换都在 O(1) 内维护 Φ（`feasible_update()`）：

```text
R' = (R \ {a,b}) ∪ {v}，v = a⊕b，d(v) = 1+max(d(a),d(b))
Φ(R') = Φ(R) − 2^{d(a)} − 2^{d(b)} + 2^{d(v)}
```

- 若 v 已在 R 中（相消）：`Φ(R') = Φ(R) − 2^{d(a)} − 2^{d(b)} − 2^{d(v)}`（三项同时消失）；
- 若 v = 常量零（a = b 的情形）：不加项。

配套的**预算增量**（`candidate_delta()`）：

```text
Δ = 2^{d(v)} − 2^{d(a)} − 2^{d(b)}
```

Δ < 0 表示这次替换**释放**了预算（合并等深项），Δ > 0 表示消耗预算（把浅项挂到深项下）。greedy 的并列排序里 Δ 越小越优先（§5.3）：等门收益下优先选择为将来的合并**保留最多深度余地**的替换。

D ≥ 63 时 2^D 溢出 uint64，代码切换到等价的显式最早就绪合并检查（`combine_depth` 试算），语义完全一致。

#### 5.2.4 combine_depth：最早就绪合并的最优性

`combine_depth(depths)` 反复取**深度最小的两项**合并，返回最终深度。它正是"固定项集合下最小化 1+max 合并深度"的最优算法——与 Huffman 编码同构（把 2^{d} 视作频率即可由 5.2.1 的 ΔΦ ≥ 0 论证：任意时刻合并非最小对都会留下更大的 Φ 终值，而 Φ 终值 = 2^{d_final} 单调对应最终深度）。它出现在三处：

- 表示的 `arrival_depth`（这个表示若现在就收尾，最少需要多少层）；
- greedy 收尾的 `earliest_ready`（把每个表达式实际合并成输出）；
- 深度可行性的通用检查。

#### 5.2.5 Kraft 预算如何引导搜索

预算把"深度约束"从**事后过滤器**变成**事前导航**：

- 任何使某行 Φ 超过 2^D 的候选根本不会生成（greedy 中 `feasible_update` 拒绝，Beam 中转移后前沿为空 → invalid）——搜索树从一开始就只在可行域内伸展；
- 预算余量 `2^D − Φ(R)` 是表示的"呼吸空间"，`candidate_delta` 让等收益候选按呼吸空间排序；
- Pareto 前沿保留的"项数多但深度低"的表示，正是**预算占用更低**的备选（浅项便宜：2^{d} 随 d 指数增长）。

### 5.3 Greedy 的选择函数（synthesize_completable）

每轮迭代，对每个表达式 R 的每一对项 (a,b) 做预算可行过滤（§5.2.3），把幸存的出现按 (a,b) 分组聚合为候选。候选的关键量：

```text
gain = |uses| − (v 为新门 ? 1 : 0)
```

`|uses|` 是该 pair 出现在多少个表达式中（每次替换使总项数减 1），新门成本 1。**gain > 0 才成为候选**——这保证每轮总项数严格下降（`require(after < before)`），主循环有限。

**排序字典序**（`(-gain, delta, depth, left, right)`）：

1. 门增益最大（净省门最多）；
2. 预算增量最小（释放/最少消耗深度余地）；
3. 新节点深度最浅；
4. 节点编号序（决定论平局打破）。

**lookahead 变体**（交替 restart）：对 gain ≥ best−1 的前 8 个候选，试算替换后每个受影响行中"下一步仍可行的后续 pair"数量 `future`，取最大值 height，评分：

```text
score = gain + max(0, height − 1)
```

这是对"这次替换能激活多少下一步共享机会"的一步前瞻。

**并列候选的随机化与 MW_SQ 重评分**：对 rank 元组 `(score, gain, −delta, −depth)` 完全并列的候选集，默认用种子随机选一（`--seed` 固定时完全可复现）。开启 `--greedy-mwsq` 时先做一步**结构重评分**（§5.4）。

### 5.4 MW_SQ 式并列重评分（mwsq_support_summary）

**动机**（源自 Plank DSN'12 的 MW2/MW_SQ）：贪心的"当前最高频 pair"是短视的——两个候选当前收益相同时，应选择**替换后剩余共享结构更集中**的那个，为后续轮次保留高质量的合并机会。

**度量**：把每个表达式视为节点 id 的行，统计所有行内无序对 (u,v) 的多重集。定义：

```text
S(rows) = ( max_support, ties_at_max )
max_support  = 出现次数最多的 pair 的出现次数
ties_at_max  = 达到该最大值的 pair 个数
```

**评估流程**：对每个并列候选 c，把 c 的替换虚拟应用到所有行（新节点用虚拟 id，不物化），对结果行集计算 S。选择字典序**最大**的 S：

1. 先比 `max_support`：替换后剩下的最热 pair 越热越好——热 pair 意味着下一轮仍有高 gain 的候选，共享机会没有被这次替换打散，而是**被集中**到少数几个 pair 上；
2. 再比 `ties_at_max`：热度相同时，并列数越**多**越好——剩余的共享结构分布在更多等价选项上，随机选择的鲁棒性更强（不会把希望押在唯一 pair 上）。

这是 MW2"most frequent pair first"准则在**深度预算约束**下的翻版：MW_SQ 原始算法用它做主选择，xordepth 只在 gain/delta/depth 全部并列时用它打破平局，因此关闭该开关时结果与旧版逐位一致。

**实现**：`mwsq_support_summary()` 用平铺 64 位键 + 排序游程统计，与 map 计数逐位等价（单元测试锁定金标准）。

### 5.5 Beam 的评估函数

两级评估：**候选级**（这一步走哪）与**状态级**（哪个状态更有前途）。

#### 5.5.1 候选级：candidates() 的两段排序

对 pair 索引扫描聚合出的每个候选 (a, b)：

```text
support(a,b)   = 边掩码 edge(a,b) 在全部表示中的出现次数
arrival_depth  = 1 + max(d(a), d(b))
progress_first = 受影响表示中"替换后剩余项数 ≤ 1"的个数（临近完成的表示数）
progress_second= 受影响表示的剩余项数总和
```

其中单个表示的剩余项数 = `|R| − (v 相消 ? 3 : 1)`（两项并一项，或三项相消）。这是 Uber-XSet 式的**完成进度**估计：`progress_first` 大意味着这次替换让许多表示逼近单叶状态（Beam 状态的出口是每个目标存在 ≤1 项的表示，`complete()`）。

**两段排序**（`candidates()`，稳定排序叠加）：

1. 第一段（决定论基序）：`support` 降序 → `arrival_depth` 升序 → `(left, right)` 升序；
2. 第二段（UBER 进度覆盖）：`support` 降序 → `progress_first` 降序 → `progress_second` 升序。

`--seed ≠ 0` 或 restart > 0 时先做确定性洗牌（`SeedRandom`），使 restart 探索不同路径。

#### 5.5.2 状态级：cost_priority

```text
cost_priority(S) = |live(S)| + max over targets of max(0, min_R |R| − 1)

live(S) = 所有前沿所有表示的项展开成生产子树后的并集节点数（不含输入）
```

即"**已物化的结构规模 + 最薄表示还差的合并数**"。它是**排序启发式**，故意包含互斥的 Pareto 备选（估计状态的物质化程度），**不可**用作下界或最优性证明（头文件注释明确区分了这一点）。

#### 5.5.3 选择性前瞻（--lookahead-width K / --depth L）

对状态的前 K 个候选各试走至多 L 步（每步取试走状态自己的首候选），用试走终点的 `cost_priority` 作为 bound 重排。这是把"状态级评估"从 1 步展宽到 L 步的受控前瞻，代价是 O(K·L) 次状态转移。

### 5.6 剪枝机制

Beam 主循环里，一个子状态要活着进入下一轮，必须闯过五道闸门：

#### 5.6.1 深度可行性剪枝（结构性）

`candidate_allowed`：`1 + max(d(a), d(b)) ≤ D`，加上转移内部的 Kraft 预算（前沿为空的表示组合直接 invalid）。计入 `pruned_depth`。这一层保证**搜索空间整体在可行域内**——不被剪掉的状态全都还能在 D 层内完成。

#### 5.6.2 分支定界（对 incumbent）

```text
若 pool 已有最优且 child.best_cost_lower_bound() ≥ pool.best.cost() → 剪掉
```

`best_cost_lower_bound()`（与 `cost_priority` 的关键区别——**保守**）：

```text
LB(S) = max over targets of  min over R∈该目标前沿 of  |closure(R)|

closure(R) = R 中所有非输入项展开成生产子树后的**互斥**闭包节点数
```

**为什么这是安全下界**（源码注释的论证）：每个表示的生产者是消0的（XOR 网络无 cancellation），单表示的闭包是从 n 个单位输入出发、逐步 pair 替换构造出的支撑集——任何计算同一目标掩码的 XOR DAG（即便允许 cancellation）都至少需要 w(target) − 1 个门；表示的互斥闭包数 ≤ w − 1 的部分即已支付门数。**对每个目标取"最便宜表示"的闭包（min）**、再对目标取 max：因为互斥的 Pareto 备选不能同时计入成本，取 min 保证不因备选而高估。

随着池内 incumbent 改善，这道闸门越来越紧——**分支定界是算法"逐步逼近"的主要收敛机制**。

#### 5.6.3 转移表（transposition）去重

子状态计算规范化签名 `signature()`：节点按 `(mask, depth)` 排序重编号后的（节点序列 + 每目标前沿的项重写）字符串。若同签名状态已以 ≤ 当前下界的值出现过 → 剪掉。等价搜索路径只算一次。

#### 5.6.4 指纹多样性分桶

所有幸存子状态按 `fingerprint()`（每目标前沿的 min 项数/min 深度/表示数/深度谱 + pair 支持谱 + 闭包规模等结构摘要）分桶；**每桶只保留 cost_priority 最小的 1 个**，然后轮转填满 beam_width 个。防止 Beam 被同构状态灌满（多样性），同时保证桶内最优。

#### 5.6.5 宽度截断

每状态至多展开 branch_width 个候选（截断发生在 5.5.1 排序之后，即**先排序后截断**——启发式排序的全部意义就在这里）。

### 5.7 逼近机制总述：各部件如何协同收敛

把上述机制串起来，算法逼近最终解的路径是：

1. **可行域内行走**（Kraft 预算）：任何时刻任何状态都满足"所有表示可完成"，所以只要走到 `complete()`（每目标有 ≤1 项的表示），导出的一定是深度合法的线路——不存在"找到解但深度超限"的回溯浪费；
2. **共享优先**（support + gain）：XOR 门数最优性的来源是**共享**——一个新门被越多的输出行复用，摊销越充分。support 降序保证优先物化"被最多表示需要的值"；
3. **完成压力**（progress_first）：驱动 Beam 状态朝出口（全目标单叶化）移动，避免在中间结构上徘徊；
4. **预算敏感**（delta / arrival_depth）：并列时选择为深度留余地的替换，避免把浅项过早挂到深链上耗尽预算；
5. **定界收紧**（5.6.2）：incumbent 每改善一次，剪枝阈值抬高一次，搜索在收缩的窗口内集中；
6. **多源冗余**（候选池）：CSE / greedy / Beam / rewrite / incumbent 是同一问题的五个不同偏置的求解器，池取字典序最优——即使某个来源在某个矩阵上失效（如纯 Beam 在 12 矩阵样本上仅 4/12 verified），整体仍然可靠；
7. **随机重启**（seed/restarts）：并列候选的种子随机化 + 多 restart 提供路径多样性，MW_SQ 重评分则把随机选择替换为结构性选择（实测在 12 样本上 3 降门 1 降深 0 回退）。

**多线程并行**（`--threads N`，N>1）改变的是调度而非结果：greedy portfolio 的各 restart 相互独立（各自 NodeStore 与种子），Beam 的每层拆为三段——Phase 1 顺序发布 complete 状态并快照每个待扩展父状态的 incumbent 阈值（复现串行路径的"发布与剪枝交错"语义）、Phase 2 各父状态独立并行扩展（候选洗牌按状态种子，天然线程安全）、Phase 3 按父顺序回放转移表查重并合并。因此任意线程数的结果与串行逐门一致（墙钟超时运行除外，它们本就依赖机器速度）。实测 sample-12 无时限 4 线程 1.66×（11.35s → 6.83s）、全量 57 矩阵 1.52×（64.73s → 42.59s）。修复 bounded 管线 greedy 阶段漏传线程数的缺陷、并进一步并行化 greedy 内部 MW_SQ 并列重评分（按行分块生成 key 并行排序、固定块序归并、按 tied 顺序串行归约，重启并行与评分并行共享线程预算）后，6 个最大矩阵合计 4 线程 1.98×（47.95s → 24.16s），全量 57 矩阵 1.68×（64.73s → 38.63s）；所有线程数的输出仍与改动前串行逐门一致，剩余串行份额主要是 Beam 的 Phase 1/3 与 CSE 基线。

### 5.8 事后重写的调度松弛（rewrite_depth_gates）

给定已验证线路 C，重写尝试**删除一个多消费者的非输出门** v 并把它的所有消费者一次性改接：

1. **反向时限表 `latest`**：从输出（latest = D）出发反向扫 ops，每个门给其父节点设 `latest(parent) = latest(child) − 1`（取更紧者）。这是 ALAP（as-late-as-possible）调度的松弛：门 u 的结果最晚要在 `latest(u)` 层就绪。它精确刻画"这个门还有多少深度余量"；
2. **受害者选择**：非输出门按消费者数 `uses` **升序**排序——消费者越少的门越容易被改接绕开；
3. **rebuild_without(v)**：跳过 v 重放所有其他门。重放某门 u 时优先找 mask 相同且就绪时间 ≤ `latest(u)` 的等值信号（**等值别名替换**）；找不到则用**补对搜索**：在已就绪信号中找 c 使得 `mask(c) ⊕ mask(u)` 也有现成信号，且新深度 `1 + max(ready(c), ready(补))` 最小、不超过 `latest(u)`；
4. **接受准则**：重写结果逐门验证功能、复验深度，且 `cost` **严格变小**才替换。最多 4 个 pass、总预算 `--rewrite-trials` 共享。

数学上这是一个受深度时限约束的**局部搜索 / 消解步骤**：它修的是贪心与 Beam 的系统性偏差——"过早物化的共享门在后期成了冗余"。

### 5.9 Rebalance（单消费者锥重排）

`rebalanced()` 对每个**单消费者、非输出**的锥按到达时间做 Huffman 重排。锥内门数严格不变（\|frontier\| − 1 前后不变），但把深链上的早熟合并推迟，深度可能下降。它是对 X-Sets 系输出的标准深度补救，在 xordepth 中只作为基线来源之一。

---

## 6. 完整示例推导

一个手工可验的小例子，贯穿 greedy 与 Kraft 预算。（本节结果已用当前二进制实测核对：同参数运行输出 `verified / cost 5 / depth 2 / best_source "greedy"`，电路 ops 与下文推导逐门一致。）

矩阵（3 输出、4 输入），D = 2：

```text
y0 = x0 ⊕ x1 ⊕ x2
y1 = x1 ⊕ x2 ⊕ x3
y2 = x0 ⊕ x1 ⊕ x2 ⊕ x3
```

### 6.1 预检

行权重：w(y0) = w(y1) = 3，w(y2) = 4。深度下界 = max(⌈log₂3⌉, ⌈log₂4⌉) = **2** ≤ D = 2 ✓，可以搜索。

初始预算（全部项深度 0，budget = 2² = 4）：

| 行 | 初始表示 R | Φ(R) = Σ2^{d} | 可行？ |
|---|---|---|---|
| y0 | {x0, x1, x2} | 1+1+1 = 3 ≤ 4 | ✓ |
| y1 | {x1, x2, x3} | 3 ≤ 4 | ✓ |
| y2 | {x0, x1, x2, x3} | 4 ≤ 4 | ✓ |

### 6.2 Greedy 轮次演算

**第 1 轮**。枚举所有 (pair, row) 出现：

| pair | 出现行 | gain = \|uses\| − 新门1 | 备注 |
|---|---|---|---|
| (x1, x2) | y0, y1, y2 | **3 − 1 = 2** | 最高 |
| (x0, x1) | y0, y2 | 2 − 1 = 1 | |
| (x0, x2) | y0, y2 | 1 | |
| (x1, x3) | y1, y2 | 1 | |
| (x2, x3) | y1, y2 | 1 | |
| (x0, x3) | y2 | 1 − 1 = 0 | ≤ 0 淘汰 |

选 (x1, x2) → 物化 n5 = x1⊕x2，深度 d = 1。替换后（含预算演算，Δ = 2¹ − 2⁰ − 2⁰ = 0，等深合并不改变 Φ）：

| 行 | R | Φ | 演算 |
|---|---|---|---|
| y0 | {x0, n5} | 1 + 2 = 3 ✓ | 3 − 1 − 1 + 2 = 3 |
| y1 | {n5, x3} | 2 + 1 = 3 ✓ | |
| y2 | {x0, n5, x3} | 1 + 2 + 1 = 4 ✓ | 预算刚好用满 |

**第 2 轮**。逐对做增量预算检查（关键步骤逐个演算）：

- **(x0, n5)** → 新节点深度 1+max(0,1) = 2，权 2² = 4。
  - 行 y0：`Φ' = 3 − 2⁰ − 2¹ + 2² = 3 − 1 − 2 + 4 = 4 ≤ 4` ✓ 可行（y0 完成成单叶 {n05}）；
  - 行 y2：`Φ' = 4 − 1 − 2 + 4 = 5 > 4` ✗ **预算拒绝**——若替换，y2 剩 {n05, x3}，Φ = 4+1 = 5，从此 y2 永远无法在 2 层内完成（Φ 单调不减）。
  - ⇒ uses 只剩 {y0}，gain = 1 − 1 = **0 ≤ 0，淘汰**。
- **(n5, x3)** → 同理：y1 可行、y2 被拒（`4 − 2 − 1 + 4 = 5`），gain = 0 淘汰。
- **(x0, x3)** → 深度 1，y2 预算 `4 − 1 − 1 + 2 = 4` ✓ 可行！但 uses = {y2} 单行，gain = 1 − 1 = 0 淘汰。

所有候选 gain ≤ 0 → **pair 交换阶段终止**（这正是贪心的近视性：第 1 轮的全局最优共享 n5 把 y2 的预算锁死在 4，任何进一步的省门替换都无处落脚）。

**收尾（earliest_ready，Huffman 合并）**。Huffman 合并不受 gain 限制，只受预算约束：

- y0 = 合并(x0 d0, n5 d1)：先排序取最早就绪两项 → 新门 g0 = x0⊕n5，深度 2 ✓（y0 = g0）；
- y1 = 合并(n5 d1, x3 d0) → g1 = n5⊕x3，深度 2 ✓；
- y2 = 合并{x0 d0, n5 d1, x3 d0}：最早就绪的是两个 d0 项 → g2 = x0⊕x3（深度 1）；再合并 g2 d1 与 n5 d1 → g3 = g2⊕n5，深度 2 ✓（等深合并，Φ 走 2¹+2¹ → 2²，始终 ≤ 4）。

**结果**：门 = {n5, g0, g1, g2, g3} 共 **5 门、深度 2**，逐门验证通过。可情形枚举验证该实例不存在 4 门解：y2 的深度 2 四叶树必须恰好用满 3 个门（两对 + 根），y0、y1 至少各需 1 个额外门，而 5 门构造恰好达到此界。另注意贪心给出的 5 门与另一条 Beam 路径（n01 = x0⊕x1、n23 = x2⊕x3、y0 = n01⊕x2、y1 = n23⊕x1、y2 = n01⊕n23，同为 5 门）互补——这正是多来源候选池的意义。

### 6.3 Beam 视角（同一实例的状态与剪枝）

初始状态：节点 = x0..x3；前沿如 6.1；pair 索引的 support：edge(x1,x2) = 3（y0/y1/y2 各一次），edge(x0,x1) = edge(x0,x2) = edge(x1,x3) = edge(x2,x3) = 2，edge(x0,x3) = 1。

- 候选排序首位：(x1, x2)，support = 3；其 progress：受影响表示剩余项数 y0: 3−1=2，y1: 2，y2: 3 → progress = (0, 7)。
- 对候选 (x1,x2) 应用 `apply` 得子状态 S₁（即 6.2 第一轮后的状态）。其 `best_cost_lower_bound`：三个目标的 closure 都是 {n5}，每目标 min = 1，max = 1（弱但安全的界）；`cost_priority` = live(1) + max(0, min\|R\|−1) = 1 + (y2 的 3−1) = 3。
- 假如此时池内 incumbent 是 6.2 的 5 门解：任何 `LB ≥ 5` 的子状态将被分支定界剪除。S₁ 的 LB = 1 ≤ 5，继续展开。
- 转移表：不同节点物化顺序但 (mask, depth) 规范化后相同的状态会被 `signature()` 去重。
- 完成判定：当某状态每个目标都有 ≤ 1 项的表示（本例：物化 n01, n23 后 y2 = {n01, n23} 再物化根门即达）→ `complete()` → 导出电路入池。

---

## 7. 使用参考

### 7.1 命令一览

```bash
cd xordepth && make -j4          # 生成 build/xordepth

xordepth bounded MATRIX.txt --depth D [选项] --output OUT.json   # 主命令
xordepth rebalance XSETS_OUTPUT.txt --matrix M.txt --output OUT.json  # X-Sets 调度导入 + 重结合
xordepth verify OUT.json                                        # 复验功能/门数/深度
xordepth merge BASE.json CAND.json... --output OUT.json         # 离线合并最优
xordepth audit FILE.json...                                     # 证据审计
```

`rebalance` 子命令解析 Plank X-Sets（UBER_XSET/MW_SQ）输出末尾的 `Schedule:` 段（`sK = xJ` 别名、`sK = sA + sB` 门、`bT = sK` 目标绑定），逐门校验其确实计算 `--matrix` 给出的矩阵，然后调用 `rebalanced()` 重结合并输出前后 `(cost, depth, max_output_depth)` 与完整线路 JSON。任何无法识别的行、先引用后定义、重复绑定或 trailer 门数不符都会报错拒绝。

### 7.2 bounded 选项与默认值

| 选项 | 默认 | 作用 | 数学对应 |
|---|---:|---|---|
| `--depth D` | 必填 | 深度硬上限 | budget = 2^D（§5.2） |
| `--beam N` | 64 | Beam 宽度（多样性桶轮转后总量） | §5.6.4 |
| `--branch N` | 16 | 每状态展开候选数（排序后截断） | §5.6.5 |
| `--representations N` | 8 | 每目标 Pareto 前沿容量 | §4.2 |
| `--restarts N` | 8 | Beam 重启轮数（洗牌种子逐轮变化） | §5.7-7 |
| `--greedy-restarts N` | 2 | greedy portfolio 次数（普通/lookahead 交替，0=纯 Beam） | §5.3 |
| `--greedy-mwsq` | 关 | 并列候选 MW_SQ 重评分（关闭=逐位旧版行为） | §5.4 |
| `--lookahead-width/depth K L` | 0/0 | 选择性前瞻 | §5.5.3 |
| `--residual-pairs/targets N M` | 0/0 | residual toggle/bridge 候选 | §4.4 |
| `--rewrite-trials N` | 0 | 重写预算（共享） | §5.8 |
| `--seed N` | 0 | 并列随机化的种子（固定则完全可复现） | |
| `--threads N` | 1 | greedy restarts 与 Beam 层内扩展的工作线程数；>1 时结果仍与串行逐位一致（发布/合并按原顺序回放） | §5.7 |
| `--timeout SEC` | 180 | 墙钟（anytime 截止不丢 incumbent） | §3.1 |
| `--no-timeout` | — | 关墙钟；有限性由 max-states/round-limit 保证 | §3.2 |
| `--max-states N` | 10000 | 全局状态预算（跨 restart） | |
| `--incumbent FILE` | — | 外部已验证线路零成本入池（可重复） | §4.5 |

### 7.3 推荐配置

```bash
# 协议对齐的可复现实验（12 矩阵 sample 协议）
--depth D --beam 4 --branch 4 --representations 4 --restarts 2 \
--greedy-restarts 2 --greedy-mwsq --rewrite-trials 400 --max-states 256 \
--seed 1 --no-timeout

# 长预算探索
--beam 64 --branch 16 --representations 8 --restarts 8 --timeout 180

# 跨运行结果汇总（零搜索时间）
xordepth merge run1.json run2.json run3.json --output best.json
```

### 7.4 结果 JSON 关键字段

`status`（见 §3.3）、`cost`/`depth`、`best_source`（哪个来源胜出）、`explored_states`/`pruned_depth`/`pruned_bound`（搜索与剪枝统计）、`toggle_candidates`/`bridge_candidates`/`diversity_buckets`（residual 与多样性统计）、`published_candidates`/`incumbents_used`/`rewrite_trials_used`（池活动统计）、`depth_lower_bound`/`round_limit`/`state_limit`（预算参数）、完整 `circuit` 与 `target_masks`（供 verify/merge/audit 复验）。

---

## 8. 正确性边界与已知限制

- **启发式上界**：`verified` 只保证功能正确、门数与深度如实计量，不保证门数最优；
- **不可行判据唯一**：只有 `infeasible_by_lower_bound` 是数学证明；`search_failed`/`timeout` 只是该预算下没找到；
- 深度按二输入 XOR 拓扑深度计；不建模 LUT、布线、扇出时延、多输入门；同一状态内不复制重复 mask 节点；
- `--greedy-mwsq` 关闭时全部路径与历史版本逐位一致（回归测试锁定）；开启时只改变并列候选的选择；
- 结果的可复现性：固定 `--seed` 与全部预算参数即逐位可复现（含 `--no-timeout`；`--timeout` 受机器速度影响）；
- 大矩阵（如 100×100）上无预算不设限运行仍很贵，靠 `--max-states`/`--timeout` 控制；预算耗尽时已有 incumbent 不丢失。

### 对照数据（摘要）

57 矩阵全套件、与 PR #6（V1/V2 greedy）及各阶段自身的对比数据与复现脚本存于 `data/results/`（详见 [pr6-sample-comparison.md](pr6-sample-comparison.md)、[matrix-suite.md](matrix-suite.md)、[runtime-optimization-2026-11-09.md](runtime-optimization-2026-11-09.md)）。要点：统一池 + MW_SQ + 离线 merge 后，对 PR 最优版本 52 个可比矩阵 **26 优 26 平 0 劣（净 −140 门）**；同预算（5s）对比 **4 优 48 平 0 劣**；无时限模式经 2026-11-09 的分配优化提速 2.18×，全部 57 矩阵结果与优化前逐位一致。
