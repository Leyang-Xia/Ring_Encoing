# C++ 工具使用说明

构建方法见 [README](../README.md)。以下命令均在仓库根目录执行；默认 `make` 生成 `build/ringlab`、`build/xsets` 和 `build/uber`。也可使用 CMake。

## 1. 核验保存的证据

```bash
./build/ringlab verify --exhaustive
./build/ringlab verify --group original --samples 1000
./build/ringlab verify --group family --output /tmp/family-verification.json
./build/ringlab verify --group practical
```

`--group` 可选 `all`、`small`、`original`、`family`、`practical`。默认每例执行 1,000 组随机输入；`--samples 0` 只略过随机执行，矩阵、线路和 MDS 检查仍会进行。`--exhaustive` 仅影响四位域小例，额外枚举全部 65,536 个输入并确认最小距离。

- 从矩阵定义重新生成原域变换、环主变换和投影，检查二进制矩阵恒等式。
- 精确解析 JSON 大整数，对每个 XOR 做符号执行，检查依赖顺序、输出及所有计入的门确实存活。
- 重新拼接主计算与投影，复核实际共享与死门消除后的计数；检查保存的分解域基线。
- 枚举全部所需的方形子式，确认系统化编码的 MDS 条件。
- 重跑小例和原始 8→17 例子的 100 个 CSE 种子，核对历史直方图及选定线路。

所有检查使用运行时错误处理，`-DNDEBUG` 不会禁用它们。符号执行与矩阵恒等式覆盖全部输入；随机执行是额外检查。

`--data-dir DIR` 指定证据目录，默认指向编译时的仓库路径。拷贝二进制到别处时，把 `xsets` 放在同目录，并传入证据目录；也可用 `--backend PATH` 指定后端。

## 2. 重新优化一个阶段

```bash
./build/ringlab profile original --path projection --method MW_SQ \
  --params 0 3 0 --timeout 180 --output /tmp/projection.json
./build/ringlab profile systematic-rs-7-4 --path direct --method UBER_XSET \
  --params 0 3 0 --timeout 60 --output /tmp/direct.json
./build/ringlab verify --input /tmp/direct.json
```

例名包括：

| 分组 | 名称 |
|---|---|
| 基础实例 | `small`、`original` |
| 未缩放 Vandermonde | `vandermonde-11d-k8`、`vandermonde-11d-k12`、`vandermonde-1d7-k8` |
| 缩放扩展 | `scaled-eight`、`scaled-ten` |
| 系统化 RS | `systematic-rs-7-4`、`systematic-rs-7-4-ring12`、`systematic-rs-7-4-powers`、`systematic-rs-17-14-preliminary` |
| Cauchy | `cauchy-7-4-0` 至 `cauchy-7-4-3` |

`--path direct` 计算固定原域输出；`main` 只计算环中间量；`projection` 从去重后的环输出信号返回原域，必要的校验混合也计入其中。进入环使用多项式代表补零，成本为 0。若以后改用其他输入映射，其实际成本必须另行计算。

`--params` 顺序为 `thresh nstart give-up-when`。它们控制后端候选枚举策略，三个数同时变大不保证搜索更强。`--timeout` 是单次后端的正数墙钟秒数，不限制整个实验的累计 CPU 时间；搜索中的 CSE 和矩阵核验也不属于该计时。

输出状态：

| 状态 | 含义 |
|---|---|
| `verified` | 已完成并逐门核验，保存 `circuit` 与实际 `cost` |
| `timeout` | 到时终止进程组，`cost: null`，没有线路结论 |
| `cancelled` | 收到 SIGINT/SIGTERM，结束并回收后端；profile 保存取消状态并非零退出 |
| `error` | 后端执行失败或输出未通过核验，保存原因，程序非零退出 |

超时结果不能作为 XOR 下界。自定义 `--backend` 不由仓库源码哈希认证；记录中的哈希标识随项目提供的 X-Sets 源码，另有实际可执行文件路径。

## 3. 搜索与比较

```bash
./build/ringlab search --family vandermonde --limit 3 --timeout 60 \
  --output /tmp/search.json
./build/ringlab verify --input /tmp/search.json --samples 1000
./build/ringlab compare --timeout 20 --output /tmp/small-comparison.json
./build/ringlab verify --input /tmp/small-comparison.json
```

`search` 细化家族证据中已有的候选定义：按已保存的初筛百分比排序，主计算使用种子 0，投影使用 10 个种子，直接路径使用 20 个种子；再分别对投影和直接路径运行两种 X-Sets 方法。它不代表重新生成整个数学参数空间。默认严格检查第一行全为 1、取点非零互异及其余行的幂关系；`--family all` 才纳入缩放扩展。

每完成一个候选就写出完整结果，保留负例和超时。取消某次后端时，保存已完成候选和当前取消测量；该文件是中断记录，不能当作完整候选证据。

`compare` 对四位域旧例重跑四套参数、两种方法和直接/环两条路径，共 16 次；完整环路线另接回原域投影并核验。所有输出要求显式 `--output`，避免长搜索结束后只剩终端摘要。

## 4. 自定义矩阵与商环

二进制矩阵文件每行一个输出、每列一个输入，支持连续 `01` 或以空白分隔的比特。第 0 列对应输入 `x0`；拒绝非 01 字符和不等长行。

```bash
./build/ringlab analyze matrix.txt --method CSE --seeds 100 --output /tmp/cse.json
./build/ringlab analyze matrix.txt --method MW_SQ --params 0 3 0 \
  --timeout 60 --output /tmp/mw.json
./build/ringlab verify --input /tmp/cse.json
```

零、单位和重复输出只涉及免费别名。后端只接收去重后的非平凡目标；返回后恢复原输出顺序，并独立校验线路。

自定义域/环定义使用 JSON。`w,P,n,M,k,m,A,lifted` 分别指定域位宽、不可约多项式、环维数、环模数、数据数、校验数、域矩阵和提升矩阵；多项式用其二进制系数对应的整数表示。要求 `P` 整除 `M`、`degree(M)=n>w` 且 `lifted mod P=A`。可选 `Q` 代替 `M`（此时 `M=P Q`）。实现支持域位宽 1 至 16、环维数小于 63；线路掩码不受 64 位限制。这是实现范围，不是理论上的限制。

[8→10 定义示例](../data/examples/original-8-to-10.json)沿用早期程序的模数，提供一个小型未缩放 Vandermonde，便于跑通：

```bash
./build/ringlab profile --input data/examples/original-8-to-10.json \
  --path main --method UBER_XSET --output /tmp/ring10-main.json
./build/ringlab profile --input data/examples/original-8-to-10.json \
  --path projection --method MW_SQ --output /tmp/ring10-projection.json
```

该定义用于展示接口，不声称有收益。其他定义形式包括家族中的单项式幂参数，以及带 `kind: systematic_rs` / `cauchy` 和 `definition` 的实际编码构造。系统化 RS 必须从评价点重新构造 `V * inverse(V_data)`；如果给出后处理矩阵，则核验其与主矩阵的乘积。实现细节见 [矩阵构造](../src/matrices.cpp)。

## 5. 开发核验与后端来源

```bash
make test
make BUILD=build-sanitize CXXFLAGS='-O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer' \
  LDFLAGS='-fsanitize=address,undefined' -j3 test
bash tests/uber_asan.sh
```

默认检查包括核心测试、全量证据、分解基线、超过 64 位的掩码、投影复用、未决状态、后端失败、超时和取消清理。独立 Uber 测试重现曾经的释放后读取；其历史实现持有进程生命周期的分配，因此该专用测试关闭退出时泄漏检测，只检查无效访问和未定义行为。

X-Sets 源码保持原样；其输入约定与版权见 [第三方说明](../third_party/README.md)。`UBER_XSET` 与独立 `uber` 程序不同。默认工具不会运行历史演示中缺失实现或硬编码外部依赖的路径。

迁移保留了历史 JSON 原始字节，并用迁移前程序的报告核对了新实现。构建和运行无需 Python；后续新增实现遵循 [C++ 开发约定](../AGENTS.md)。
