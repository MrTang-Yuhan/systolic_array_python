# 脉动阵列（Systolic Array）

## 1. 介绍

本目录提供一个 **二维脉动阵列** 的 Python 仿真，用于计算矩阵乘法。

- 输入矩阵：`A` 的形状是 `(m, k)`，`B` 的形状是 `(k, n)`
- 输出矩阵：`C = A @ B`，形状是 `(m, n)`

阵列由 `m × n` 个处理单元（PE）组成。每个 PE 在每个周期执行一次乘加（MAC）：

- `acc[i, j] = acc[i, j] + in_a * in_b`
- `in_a` 向右传递
- `in_b` 向下传递

![systolic_array](./img/systolic.png)

---

## 2. 关键规律（与代码一致）

对于 `A(m×k)` @ `B(k×n)`：

- **输入倾斜（skew）规则**
  - `A[i, p]` 在周期 `t = i + p` 注入第 `i` 行左侧
  - `B[p, j]` 在周期 `t = j + p` 注入第 `j` 列上侧

- **脉动阵列计算矩阵总耗费周期数**
  - `total_cycles = m + n + k - 2`

---

## 3. 代码说明（[`systolic_array.py`](./systolic_array.py)）

- **`PE (脉动阵列cell)`**
  - `compute(in_a, in_b)`：执行乘加，保存下一拍输出
  - `clock_tick()`：更新 `out_a` 和 `out_b`

- **`generate_skewed_inputs(A, B)`**
  - 根据 `A`、`B` 生成倾斜输入序列 `A_inputs` 与 `B_inputs`
  - 返回 `A_inputs`、`B_inputs`、`total_cycles`

- **`simulate_systolic_array(A, B)`**
  - 构建 `m × n` 的 PE 网格
  - 每个周期执行“先 compute、后 clock_tick”
  - 返回 `C`、`history`、`A_inputs`、`B_inputs`

- **`print_cycle_history(history)`**
  - 打印每个周期的 PE 累加值

- **`plot_systolic_history(history, m, n)`**
  - 可视化每个周期的数据输入、PE 状态和数据流向

---

## 4. 使用方式

```powershell
python .\systolic_array.py
```

运行后会：

1. 生成并打印 `A_inputs` / `B_inputs`
2. 打印每个周期的累加状态
3. 输出 `C` 与 `A @ B` 对比结果
4. 绘制每个周期的数据流动图

