# Mini Robust Control（迷你鲁棒控制）

**从零开始、零依赖的 C 语言实现**，涵盖大学级鲁棒控制理论。每个模块对应 MIT、Stanford、Caltech 等顶尖大学的课程，将教科书中的公式转化为可运行的 C 代码，实现理论与实践的桥接——从小增益定理、H∞ 综合到 μ 分析、无源性控制和 Kharitonov 多项式框架。

## 子模块

| 子模块 | 主题 | 参考课程 |
|--------|------|----------|
| [mini-gap-metric-robustness](mini-gap-metric-robustness/) | 归一化互质分解、间隙度量、回路成形鲁棒稳定性 | MIT 6.241J, Caltech CDS 110 |
| [mini-hinf-synthesis](mini-hinf-synthesis/) | 有界实引理、H∞ 范数计算、Riccati 方程求解 | MIT 16.30, Stanford AA278 |
| [mini-kharitonov-theorem](mini-kharitonov-theorem/) | Kharitonov 多项式、区间多项式 Hurwitz 稳定性、参数鲁棒性验证 | MIT 6.241J, UT Austin ME 381R |
| [mini-mu-synthesis](mini-mu-synthesis/) | 结构奇异值 μ、D-K 迭代、混合 μ 综合、鲁棒性能 | Stanford AA278, Caltech CDS 110 |
| [mini-parametric-uncertainty](mini-parametric-uncertainty/) | 多面体不确定性、LMI 稳定性条件、基于 Kharitonov 的鲁棒分析 | MIT 16.30, MIT 16.31 |
| [mini-passivity-based-control](mini-passivity-based-control/) | 无源性定理、能量整形、端口哈密顿系统、IDA-PBC | MIT 6.241J, Imperial College EE4-48 |
| [mini-small-gain-theorem](mini-small-gain-theorem/) | 小增益定理、输入-输出稳定性、互联鲁棒性 | MIT 6.241J, Caltech CDS 110 |
| [mini-structured-singular-value](mini-structured-singular-value/) | 结构奇异值（μ）分析、上/下界、LFT 不确定性建模 | Stanford AA278, MIT 16.30 |

## 设计理念

- **零外部依赖** — 纯 C（C99/C11），仅使用 `libc` 和 `libm`
- **模块自包含** — 每个目录自带 `Makefile`、`include/`、`src/`、`examples/`、`demos/`、`tests/`
- **理论到代码的映射** — 每个模块包含 `docs/` 目录，内有课程对齐说明及教材引用（Zhou–Doyle–Glover、Skogestad–Postlethwaite、Green–Limebeer）
- **实用演示程序** — 飞机俯仰控制、直流电机鲁棒调参、四旋翼稳定性分析、磁悬浮控制、功率变换器调节

## 构建方式

每个模块相互独立。进入模块目录后运行：

```bash
cd mini-gap-metric-robustness
make all    # 构建全部
make test   # 运行测试
```

需要 **GCC** 和 **GNU Make**。

## 项目结构

```
mini-robust-control/
├── mini-gap-metric-robustness/    # 互质分解、间隙度量、回路成形鲁棒性
├── mini-hinf-synthesis/           # H∞ 范数、有界实引理、Riccati 求解器
├── mini-kharitonov-theorem/       # Kharitonov 定理、区间多项式 Hurwitz 稳定性
├── mini-mu-synthesis/             # μ 综合、D-K 迭代、鲁棒性能
├── mini-parametric-uncertainty/   # 多面体不确定性、LMI 条件、参数稳定性
├── mini-passivity-based-control/  # 无源性定理、能量整形、IDA-PBC
├── mini-small-gain-theorem/       # 小增益定理、互联鲁棒性
└── mini-structured-singular-value/ # 结构奇异值 μ、LFT 不确定性
```

## 许可证

MIT
