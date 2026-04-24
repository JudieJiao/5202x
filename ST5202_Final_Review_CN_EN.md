# ST5202 Applied Regression Analysis — 期末复习资料 (Final Review)
## 覆盖范围 Scope: Lecture 4 (Multiple Comparisons) → Lecture 9 (General Linear Models)

> **考试形式 Exam Format**: Pen-and-paper, open-book, calculator allowed, 不允许 R 运行。
> **题型 Question Format**: 读 R 输出 (R output reading) + 理解 (understanding) + 手算 (calculation)。
> **核心技能 Core Skill**: 能从 `lm()` / `glm()` 的 Coefficients / Anova / Var-Cov matrix 中提取数字,代入正确公式,得出结论。

---

## 0. 总览:六大知识模块的勾稽关系 (How the modules connect)

```
              ┌───────────────────────────────────────────┐
              │  基础模型 y = Xβ + ε   ε ~ N(0, σ²I)       │
              └───────────────────────────────────────────┘
                                │
   ┌─────────────┬──────────────┼─────────────┬─────────────┐
   ▼             ▼              ▼             ▼             ▼
 L4 多重比较   L5 一般线性   L6 模型选择   L7 回归诊断    L8–9 补救
 Multiple      General        Model         Regression     Remedies
 Comparisons   Linear Tests   Selection     Diagnostics    & GLM
   │             │              │             │             │
   │  Contrast    F-test        AIC/BIC/Cp    残差/杠杆/    WLS/Ridge/
   │  (L, c=Lβ)   Full vs       Forward/      Cook's D /    LASSO/PCR
   │  Bonf/Tukey  Reduced       Backward/     QQ plot       MCP/SCAD
   │  Scheffé     SSE_R-SSE_F   Stepwise      VIF           Box-Cox
   ▼
 控制 FWER
```

**一条主线 Main thread**: 先拟合 (fit) → 发现问题 (诊断 diagnose) → 补救 (remedy) → 推断 (inference,含多重比较)。

---

## 1. Lecture 4: Multiple Comparisons 多重比较

### 1.1 核心问题 Core problem
当同时做 **m 个** 假设检验时,每个单独用 α 会导致 **overall Type I error rate**(即 **family-wise error rate, FWER**) ≫ α。需要调整 critical value。

### 1.2 四种准则 Four criteria (★必考★)

| 准则 Criterion | 适用范围 When | Critical value (two-sided) |
|---|---|---|
| **Bonferroni** | 任何有限多个预先指定的检验 (any finite pre-specified contrasts / tests) | **$t_{n-p}(\alpha/(2m))$** |
| **Scheffé** | **所有可能的 contrasts** (any contrast, even data-driven) among q group means | **$\sqrt{(q-1)\, F_{q-1,\,n-p}(\alpha)}$** |
| **Tukey (HSD)** | **成对比较 pairwise comparisons** $\mu_i - \mu_j$ among q means (**需要平衡 equal sample sizes**) | **$\dfrac{q_{q,\,n-p}(\alpha)}{\sqrt{2}}$** (q = Studentized range) |
| **Dunnett** | **比较 vs 一个对照组 (vs control)** | Dunnett's 表格,单/双边 |

> **公式口诀 (重要!)**:  
> **`|t*| > c` 拒绝**, 其中 $t^* = \dfrac{\hat{L}}{SE(\hat{L})}$,`c` 由上表给出。

### 1.3 Contrast (对比) 的构造
一个 contrast 是 $L = \sum_i c_i \mu_i$,其中 **$\sum c_i = 0$**。
- Pairwise: $\mu_i - \mu_j$ → $(c_1,\dots)=(1,-1,0,\dots)$
- 一对多: $\mu_2 - (\mu_1+\mu_3)/2$ → $(-1/2, 1, -1/2)$

**同步置信区间 (Simultaneous CI)**:  
$$\hat{L} \pm c \cdot SE(\hat{L})$$
其中 `c` 按选用的准则取 critical value。

### 1.4 如何挑选准则 (past exam 高频题 ★)
**一次题目通常问: "Bonferroni / Scheffé / Tukey 哪个最好?"**
逻辑判断:
1. **检验的是 pairwise** 且 **m 较大** → 比 Bonferroni 和 Tukey 的 c 值,**Tukey 通常更紧 (smaller critical value) → better power**。
2. **检验的是一般 contrasts** → Scheffé 保守但通吃;若 m 很小,Bonferroni 可能更紧。
3. **Tukey 要求 equal n**, 不满足就不能用。
4. **criterion 越小 → CI 越短 → 越优**。

### 1.5 单侧调整 (one-sided adjustment)
若题目说 "one-sided",Bonferroni 用 $t_{n-p}(\alpha/m)$(不再除 2)。
Scheffé 不易单侧化;Tukey 的 Studentized range 本质是双侧的,单侧时乘 $1/\sqrt{2}$ 要谨慎使用(考题常让你只在 Bonferroni 上调整)。

### 1.6 R 代码要点 R code
```r
pairwise.t.test(y, group, p.adjust.method = "bonferroni")
TukeyHSD(aov(y ~ group))                 # Tukey HSD
# Scheffé: 手动, F_{q-1, n-p}(alpha) → 临界值
library(multcomp)                        # Dunnett
summary(glht(fit, linfct = mcp(group = "Dunnett")))
```

### 1.7 典型考题示例 (Past paper Q1)
> 已知 3 brands,`lm(y~brand)` 输出 brand2 estimate=-15.8 (SE=2.498), brand3=5.2 (SE=2.498),n-p = 12。预期 μ₂<μ₁<μ₃。
> - H0₁: β₂ ≥ 0 vs H1₁: β₂<0 (即 μ₂<μ₁)
> - H0₂: β₃ ≤ 0 vs H1₂: β₃>0 (即 μ₃>μ₁)
> **test statistics**: t₁ = -15.8/2.498 = -6.325; t₂ = 5.2/2.498 = 2.082
> **Bonferroni one-sided critical**: $t_{12}(0.05/2) = t_{12}(0.025) = 2.1788$
> **Decision**: |t₁| = 6.325 > 2.1788 → 拒绝 H0₁; |t₂| = 2.082 < 2.1788 → 不拒绝。
> **iia**: Tukey 不能直接做 directional contrast 除非转化为 pairwise;Scheffé 可以;Bonferroni 可以。
> **iib**: 比 c 值 → 选小的。

---

## 2. Lecture 5: General Linear Test & Contrasts & Interaction 一般线性检验、对比、交互

### 2.1 一般线性检验 (General Linear F-test) ★高频公式★
对假设 $H_0: L\beta = 0$ (q 个线性约束):
$$\boxed{F^* = \dfrac{(SSE_R - SSE_F)/q}{SSE_F/(n-p)} \;\sim\; F_{q,\,n-p} \text{ under } H_0}$$

- **SSE_F** = Full model residual sum of squares
- **SSE_R** = Reduced model (imposing $H_0$) residual sum of squares
- **q** = df of numerator = number of linear restrictions (= df reduction)
- **n − p** = df of Full model residual

**等价形式 (Wald form)**:
$$F^* = \dfrac{(L\hat\beta)^\top [L(X^\top X)^{-1}L^\top]^{-1}(L\hat\beta)}{q\,\hat\sigma^2}$$

### 2.2 交互模型 (Interaction model) 考点 (Past paper Q2)
模型 $y = \beta_0 + \beta_2 T_2 + \beta_3 T_3 + \alpha z + \gamma_1 T_2 z + \gamma_2 T_3 z + \epsilon$

**检验交互是否显著**:
- $H_0: \gamma_1 = \gamma_2 = 0$, 即 $L = \begin{pmatrix} 0&0&0&0&1&0 \\ 0&0&0&0&0&1 \end{pmatrix}$
- 用 partial F 或从 var-cov matrix 构造:  
  $F^* = \frac{1}{2\hat\sigma^2} \hat\gamma^\top \Sigma_{\gamma\gamma}^{-1} \hat\gamma$
- **2×2 矩阵求逆提示**: $\begin{pmatrix}a&b\\c&d\end{pmatrix}^{-1}=\frac{1}{ad-bc}\begin{pmatrix}d&-b\\-c&a\end{pmatrix}$

### 2.3 计算 $\text{Var}(\hat L)$ 与 $SE(\hat L)$ ★
已知 var-cov matrix $\Sigma = \hat\sigma^2 (X^\top X)^{-1}$,则对 contrast $L\hat\beta$:
$$\boxed{\text{Var}(L\hat\beta) = L\,\Sigma\,L^\top, \quad SE(L\hat\beta) = \sqrt{L\,\Sigma\,L^\top}}$$

例如 contrast $\mu_2-\mu_1 = \beta_2$ (dummy coding): $Var = \Sigma_{22}$。
Contrast $\mu_2 - (\mu_1+\mu_3)/2 = \beta_2 - \beta_3/2$: $Var = \Sigma_{22} - \Sigma_{23} + \Sigma_{33}/4$。

### 2.4 Dummy variables 虚拟变量
**q levels → q-1 dummies** (conventional / reference coding, level 1 = baseline)。
- Intercept = $\mu_1$ (baseline mean)
- brand2 coef = $\mu_2 - \mu_1$
- brand3 coef = $\mu_3 - \mu_1$

### 2.5 ANCOVA (协方差分析)
模型同时含 factor 和 covariate。"when x is held at its average value" 意味着 $z = x - \bar x = 0$,此时 $\mu_j$ 只由 treatment 部分决定。

### 2.6 R 代码
```r
fit_f <- lm(y ~ T + z + T:z)       # Full
fit_r <- lm(y ~ T + z)             # Reduced (no interaction)
anova(fit_r, fit_f)                # General linear F-test
vcov(fit_f)                        # Variance-covariance matrix
```

---

## 3. Lecture 6: Model Selection 模型选择

### 3.1 准则 (Criteria)
| 准则 | 公式 | 方向 |
|---|---|---|
| **AIC** | $n\ln(SSE/n) + 2p$ | 越小越好 smaller better |
| **BIC** | $n\ln(SSE/n) + p\ln n$ | 越小越好;惩罚更重 |
| **Mallows' Cp** | $SSE_p/\hat\sigma^2_{\text{full}} - (n - 2p)$ | 接近 p 越好 |
| **Adj. R²** | $1 - \frac{(1-R^2)(n-1)}{n-p}$ | 越大越好 |
| **PRESS** | $\sum (y_i - \hat y_{(i)})^2$ | 越小越好 (LOO-CV) |

### 3.2 三种逐步算法 (Stepwise)
- **Forward selection**: 从 null 开始,每步加入 **使 SSE 下降最大** (等价于 partial F 最大 / AIC 最小)的变量。
- **Backward elimination**: 从 full 开始,每步剔除最不显著的。
- **Stepwise (both)**: 每步先考虑加,再考虑删。

### 3.3 读 R stepAIC 输出 ★ (Past paper Q3)
```
Step: AIC=276.73
y ~ x3
        Df     SSR   SSE
+ x2     1   52502   2210      ← SSR 最大 / SSE 最小 → 候选
+ x4     1   10639  44073
<none>                 54712   ← 当前模型的 SSE
+ x1     1    2887  51825
```
- **SSR** = sum of squares **reduction** (加入该变量后 SSE 下降量) = $SSE_{\text{old}} - SSE_{\text{new}}$。
- **候选**: SSR 最大的 x2。
- **AIC 判断是否能加**: 计算加入 x2 后的 AIC,若小于 276.73 则加入。  
  $\Delta \text{AIC} = n\ln(SSE_{\text{new}}/SSE_{\text{old}}) + 2$

### 3.4 两种起点 (Two starting models, Tutorial 9 Q2)
- **null1**: `lm(y~1)`,从零变量开始
- **null2**: `lm(y~x10)`,强制保留 factor region

结果可能不同,用 **AIC** 比较两个最终模型,取 AIC 较小者。

### 3.5 R 代码
```r
library(MASS)
null <- lm(y ~ 1, data=d); full <- lm(y ~ ., data=d)
stepAIC(null, scope=list(upper=full, lower=null), direction="forward")
stepAIC(full, direction="backward")
stepAIC(null, scope=list(upper=full), direction="both")
```

---

## 4. Lecture 7: Regression Diagnostics 回归诊断 ★重点章节★

### 4.1 残差三兄弟 (Three types of residuals) ★必背★
令 $H = X(X^\top X)^{-1}X^\top$, $h_{ii}$ = 第 i 个杠杆值 (leverage)。

| 名称 | 公式 | 用途 |
|---|---|---|
| **Raw residual**   | $e_i = y_i - \hat y_i$ | 基础 |
| **Standardized / internally studentized** | **$r_i = \dfrac{e_i}{\hat\sigma\sqrt{1-h_{ii}}}$** | ~N(0,1),看 homoscedasticity/normality |
| **Studentized deleted / externally studentized** | **$t_i = \dfrac{e_i}{\hat\sigma_{(i)}\sqrt{1-h_{ii}}} = r_i\sqrt{\dfrac{n-p-1}{n-p-r_i^2}}$** | ~ $t_{n-p-1}$, 严格检测 outlier (inconsistency) |

### 4.2 杠杆 Leverage $h_{ii}$
- $h_{ii} \in [0,1]$, $\sum h_{ii} = p$, **平均 = p/n**。
- **高杠杆阈值**: $h_{ii} > 2p/n$ (rule of thumb)。
- 杠杆只看 X 空间,不看 y。

### 4.3 Cook's Distance (影响) Cook's D ★
$$\boxed{D_i = \dfrac{r_i^2}{p}\cdot\dfrac{h_{ii}}{1-h_{ii}}}$$

- 综合了 **inconsistency (r_i²)** 和 **leverage (h_ii)** → 衡量 **influence**。
- **阈值 rule of thumb**: $D_i > 1$ 严重, $D_i > 0.5$ 关注; Tutorial 中常用 **$D_i > 0.1$** 做筛选。

### 4.4 三个诊断 (重要概念区分 ★)
| 概念 | 英文 | 度量 | 解释 |
|---|---|---|---|
| **高杠杆** | high leverage | $h_{ii}$ | X 空间中位置极端 |
| **不一致** | inconsistency (outlier) | $|t_i|$ (studentized deleted) | y 偏离模型预测 |
| **高影响** | high influence | $D_i$ | 删掉后对 $\hat\beta$ 影响大 |

> **记忆**: leverage = **X 异常**; outlier = **y 异常**; influence = **两者结合决定删除后的变化**。

### 4.5 诊断流程 (考试流程题 ★)
1. `plot(fitted, residuals)` → 看 **homogeneity** (应水平散布,漏斗形 → 异方差)
2. `qqnorm(standardized residuals)` → 看 **normality** (偏离直线 → 非正态)
3. `qqnorm(hat)`, `qqnorm(rstudent)`, `qqnorm(cook)` → 找极值点
4. **formal test** for outlier: 创建 dummy u_i, 拟合 `lm(y ~ X + u_i)` , 检验 u_i 的 t 值
5. 若 m 个可疑观测同时检验 → **Bonferroni**: critical = $t_{n-p-m}(\alpha/(2m))$

### 4.6 Past paper Q4 (iv) 解析
三个 dummy u,v,w 同时检验 3 个点。控制 FWER = 0.05:
- Bonferroni one-sided/two-sided critical = $t_{n-p-3}(0.05/(2\cdot 3)) = t_{df}(0.00833)$
- 比较每个 |t value| 和 critical value
- 仅 |t_u| = 3.947 超过 → 只有观测 12 是 significant outlier

### 4.7 R 代码 (Tutorial 9/10 都用)
```r
fit = glm(y ~ X, data=d)
yhat = fit$fitted.values
r    = residuals(fit, type="pearson")
h    = hatvalues(fit, type="diagonal")
infl = influence(fit, do.coef=FALSE)
rsta = rstandard(fit, infl, type="pearson")   # 内学生化
rstu = rstudent(fit, infl, type="pearson")    # 外学生化 (studentized deleted)
cook = cooks.distance(fit, infl, res=infl$pear.res,
                      dispersion=summary(fit)$dispersion, hat=infl$hat)

# 诊断四图
par(mfrow=c(2,2))
plot(yhat, r)                              # homogeneity
qqnorm(rsta)                               # normality
qqnorm(h);    qqnorm(rstu);  qqnorm(cook)  # outliers 3 measures

# 阈值筛查
p=...; n=length(h)
h[h > 2*p/n]                               # high leverage
rstu[abs(rstu) > 2*sd(rstu)]               # inconsistent
cook[cook > 0.1]                           # influential

# Bonferroni 同时检验 m 个 outlier
qt(alpha/(2*m), n-p-m, lower.tail=FALSE)
```

---

## 5. Lecture 8: Remedies — Non-constant Variance & Transformations 补救 — 异方差与变换

### 5.1 Weighted Least Squares (WLS) 加权最小二乘 ★
当 $\text{Var}(\epsilon_i) = \sigma_i^2$ 非常数,取权重 $w_i = 1/\sigma_i^2$,估计:
$$\boxed{\hat\beta_W = (X^\top W X)^{-1} X^\top W y, \qquad W = \text{diag}(w_1,\dots,w_n)}$$

### 5.2 估计权重 (Estimating weights) — 两种方法 (Tutorial 10)
两种方法都基于 OLS 残差 $|e_i|$ 对 $\hat y_i$ 的依赖:

**(A) Log-log regression**:
$$\log|e_i| = \gamma_0 + \gamma_1 \log(\hat y_i) + u_i \Rightarrow \hat s_i = \exp(\hat{\log}|e_i|), \; w_i = 1/\hat s_i^2$$

**(B) Linear regression**:
$$|e_i| = \gamma_0 + \gamma_1 \hat y_i + u_i \Rightarrow \hat s_i = \hat{|e_i|}, \; w_i = 1/\hat s_i^2$$

**R 代码**:
```r
fit = lm(y ~ X)
r = fit$resid; yhat = fit$fitted
# 方法 A
wreg = lm(log(abs(r)) ~ log(yhat)); s = exp(wreg$fitted); w = 1/s^2
fit.w = lm(y ~ X, weights = w)
# 诊断加权后的残差: sqrt(w)*r vs sqrt(w)*yhat 是否平稳
```

### 5.3 Box-Cox 变换
$$y^{(\lambda)} = \begin{cases}(y^\lambda - 1)/\lambda, & \lambda\ne 0\\ \log y, & \lambda = 0\end{cases}$$
选使 log-likelihood 最大的 λ。常见: λ=0 (log), 0.5 (sqrt), 1 (none), -1 (1/y)。

```r
library(MASS)
boxcox(lm(y ~ X))
```

### 5.4 方差稳定化变换 (Variance stabilizing)
| y 的类型 | Var(y) ∝ | 建议变换 |
|---|---|---|
| 计数 Count | μ | $\sqrt y$ |
| 比例 Proportion | μ(1-μ) | $\arcsin\sqrt y$ |
| 右偏 | μ² | log y |
| 方差 ∝ μ⁴ | μ⁴ | 1/y |

---

## 6. Lecture 9: Multicollinearity Remedies & Regularization 多重共线性补救 & 正则化

### 6.1 VIF (Variance Inflation Factor) ★核心公式★
对第 j 个预测:
$$\boxed{\text{VIF}_j = \dfrac{1}{1-R_j^2}}$$

其中 $R_j^2$ 来自 **auxiliary regression**: $x_j$ 对其他所有预测变量的回归 R²。

**等价形式**(从 R 的 `lm` 输出直接算):
$$\text{VIF}_j = \dfrac{\text{Var}(\hat\beta_j)}{\hat\sigma^2 / \left[\sum_i (x_{ji}-\bar x_j)^2\right]} = \dfrac{SE(\hat\beta_j)^2 \cdot \sum_i(x_{ji}-\bar x_j)^2}{\hat\sigma^2}$$

- VIF > 10 (或 5) → 严重 multicollinearity
- **Past paper Q3(iia)** 就让你从 SE、σ̂、$\sum(x_{ji}-\bar x_j)^2$ 手算 VIF。

### 6.2 两个预测时两 VIF 相等 (Past paper Q3 iib) ★
当只有 2 个预测 x₂, x₃: 互相回归的 R² 是对称的:$R_2^2 = R_3^2 = r_{23}^2$(皮尔逊相关系数平方)。故 VIF₂ = VIF₃ = 1/(1 − r²)。

### 6.3 Ridge Regression 岭回归
$$\boxed{\hat\beta^{\text{ridge}} = \arg\min_\beta \left\{\|y-X\beta\|^2 + \lambda\|\beta\|_2^2\right\} = (X^\top X + \lambda I)^{-1}X^\top y}$$

- λ 由 **GCV (Generalized Cross-Validation)** 选,取 GCV 最小值的 λ。
- 所有 coef 向 0 收缩但**不为 0** (不做变量选择)。
- **标准化 X** 后再做。

```r
library(MASS)
r = lm.ridge(y ~ ., data=d, lambda=seq(0,50,0.01))
plot(r$lambda, r$GCV, type="l")
lambda.best = r$lambda[which.min(r$GCV)]
```

### 6.4 LASSO
$$\boxed{\hat\beta^{\text{lasso}} = \arg\min_\beta \left\{\tfrac{1}{2n}\|y-X\beta\|^2 + \lambda\|\beta\|_1\right\}}$$

- $\ell_1$ 惩罚会产生**稀疏解 sparse** → **自动变量选择** (variable selection)。
- 用 **10-fold CV** 选 λ,`lambda.min` (最小 CV error) 或 `lambda.1se` (保守)。

```r
library(glmnet)
cv <- cv.glmnet(X, y, alpha=1, foldid=foldid)   # alpha=1: LASSO, alpha=0: Ridge
plot(cv)
cv$lambda.min                                   # 选中的 λ
coef(cv, s="lambda.min")                        # 系数
pred <- predict(cv, newx=X_test, s="lambda.min")
mean((y_test - pred)^2)                         # Test MSE
```

### 6.5 MCP & SCAD (非凸惩罚)
都是为解决 LASSO bias 而设计,对大 |β| 的惩罚 plateau。
- **MCP**: penalty plateau 在 $|\beta| > \gamma\lambda$ 时完全不惩罚。
- **SCAD**: Smoothly Clipped Absolute Deviation,三段式。

特点: **比 LASSO 更稀疏 / 估计偏差更小**(大系数)。

```r
library(ncvreg)
cv.ncvreg(X, y, penalty="MCP",  fold=foldid)
cv.ncvreg(X, y, penalty="SCAD", fold=foldid)
```

### 6.6 Principal Component Regression (PCR) 主成分回归
1. 对 X 做 PCA,得到主成分 $Z_1,\dots,Z_p$。
2. 用前 k 个 Z 做 OLS: $y = \alpha_0 + \alpha_1 Z_1+\cdots+\alpha_k Z_k + \epsilon$。
3. k 由 CV MSEP 选。

```r
library(pls)
pcr.fit = pcr(y ~ ., data=d, scale=TRUE, validation="CV")
validationplot(pcr.fit, val.type="MSEP")
```

### 6.7 四种方法对比 (★高频考点★)
| 方法 | 变量选择? | 偏差 | 适用 |
|---|---|---|---|
| **OLS** | 否 | 无偏 | p << n, 无多重共线 |
| **Ridge** | 否 | 有 | 多重共线,保留所有变量 |
| **LASSO** | 是 (sparse) | 有 (大 β 时偏) | 高维、想选变量 |
| **MCP/SCAD** | 是 (sparser) | 小 | 高维、大系数不想被缩太多 |
| **PCR** | 否 (降维) | 有 | 严重共线,预测为主 |

### 6.8 Past paper Q4 (v) 思路
- VIF=7.03 → 共线性严重,**$\hat\beta_1$ 对 x1 单独效应的推断不可靠 (inference unreliable)**,但 $\hat y$ 的预测仍可。
- **若目的是预测**: 可用 ridge / LASSO / PCR 降 SE。

---

## 7. Lecture 9(别名): General Linear Models (GLM) 广义线性模型 (若考)

> 注:Ng 教案中 "General Linear Models" 可能指 GLM(Logistic / Poisson)。若课件侧重 ordinary 'general' linear model(即 matrix form $y=X\beta+\epsilon$),以上 L5 内容已覆盖。以下简要备忘。

### 7.1 GLM 三要素 Three components
1. **Random**: $y_i$ ~ exponential family
2. **Systematic (linear predictor)**: $\eta_i = x_i^\top\beta$
3. **Link function**: $g(\mu_i) = \eta_i$

| 模型 | 分布 | link | 场景 |
|---|---|---|---|
| Gaussian (OLS) | Normal | identity | 连续 y |
| **Logistic** | Binomial | logit $\log\frac{p}{1-p}$ | y ∈ {0,1} |
| **Poisson** | Poisson | log | count y |

### 7.2 Deviance 偏差
$$D = 2[\ell(y;y) - \ell(\hat\mu;y)] \quad (\text{饱和模型 − 当前模型})$$
- Residual deviance ~ $\chi^2_{n-p}$ (渐近)
- **模型比较**: $D_R - D_F \sim \chi^2_q$

### 7.3 R 代码
```r
glm(y ~ X, family=binomial, data=d)
glm(y ~ X, family=poisson, data=d)
anova(fit_r, fit_f, test="Chisq")
```

---

## 8. R 输出逐行解读手册 (R output cheat-sheet)

### 8.1 `summary(lm)` 输出
```
Coefficients:
             Estimate   Std.Error   t value    Pr(>|t|)
(Intercept)  95.200     1.766       53.896     1.1e-15 ***
brand2      -15.800     2.498       -6.325     3.8e-05 ***
```
- Estimate = $\hat\beta_j$
- Std.Error = $SE(\hat\beta_j) = \sqrt{\hat\sigma^2 [(X^\top X)^{-1}]_{jj}}$
- **t value = Estimate / SE** → ~ $t_{n-p}$ under $H_0:\beta_j=0$
- p-value → **双侧**!要单侧须除以 2 后与 α 比。

```
Residual standard error: 9.803 on 23 degrees of freedom
```
- $\hat\sigma = 9.803$, **df = n − p = 23**, 故 **n = 23 + p**。

```
F-statistic: 918.3 on 2 and 23 DF, p-value: <2.2e-16
```
- 整体显著性 overall F-test: $H_0: \beta_1=\cdots=\beta_{p-1}=0$
- $F = \frac{MSR}{MSE} = \frac{SSR/(p-1)}{SSE/(n-p)}$

### 8.2 `anova(lm)` 输出
序贯 (Type I) SS: 按输入变量顺序依次 reduce SSE。
**注意**: 顺序改变 → partial SS 改变(非正交)。

### 8.3 Variance-covariance matrix
对 contrast $L\hat\beta$:
$$\text{Var}(L\hat\beta) = L\,\text{vcov}\,L^\top$$
**Past paper Q2 (ii)**: $\mu_2 - \mu_1 = \beta_2$,  $\text{Var}= 1.9539$, SE = $\sqrt{1.9539}=1.398$。

---

## 9. 高频公式速查 ★★★ (Formula Cheat Sheet)

### 核心分布
- $\hat\beta \sim N(\beta, \sigma^2(X^\top X)^{-1})$
- $\dfrac{\hat\beta_j-\beta_j}{SE(\hat\beta_j)} \sim t_{n-p}$
- $\dfrac{(n-p)\hat\sigma^2}{\sigma^2} \sim \chi^2_{n-p}$

### 预测
- **Mean response** at $x_0$: $\text{Var}(\hat y_0) = \hat\sigma^2 x_0^\top(X^\top X)^{-1}x_0$
- **New obs**: $\text{Var} = \hat\sigma^2[1+x_0^\top(X^\top X)^{-1}x_0]$

### 多重比较(必背)
| 类型 | Critical value (two-sided) |
|---|---|
| Single t | $t_{n-p}(\alpha/2)$ |
| Bonferroni (m tests) | $t_{n-p}(\alpha/(2m))$ |
| Scheffé (any contrast, q groups) | $\sqrt{(q-1)F_{q-1,n-p}(\alpha)}$ |
| Tukey pairwise | $q_{q,n-p}(\alpha)/\sqrt 2$ |
| Dunnett | 查表 |

### 一般线性 F 检验
$$F = \dfrac{(SSE_R - SSE_F)/q}{SSE_F/(n-p_F)} \sim F_{q,\,n-p_F}$$

### 诊断度量
- Leverage: $h_{ii}$, threshold $2p/n$
- Standardized residual: $r_i = e_i/[\hat\sigma\sqrt{1-h_{ii}}]$
- Studentized deleted: $t_i = r_i\sqrt{(n-p-1)/(n-p-r_i^2)} \sim t_{n-p-1}$
- Cook: $D_i = \frac{r_i^2}{p}\cdot\frac{h_{ii}}{1-h_{ii}}$

### 共线性
$$\text{VIF}_j = \dfrac{1}{1-R_j^2} = \dfrac{SE(\hat\beta_j)^2 \cdot \sum(x_{ji}-\bar x_j)^2}{\hat\sigma^2}$$

### 正则化目标函数
- Ridge: $\|y-X\beta\|^2 + \lambda\|\beta\|_2^2$
- LASSO: $\tfrac{1}{2n}\|y-X\beta\|^2 + \lambda\|\beta\|_1$

---

## 10. 题型与策略 (Question types & strategies)

### 题型 A: Multiple Comparison 构造检验
**策略**:
1. 写 **Hypotheses** (用 β 或 μ): "H0: μ₂-μ₁ ≥ 0, H1: μ₂-μ₁<0"
2. 从 Coefficients 表找 $\hat L$ 和 $SE(\hat L)$ → 算 $t^*$
3. 查表得 critical value (注意 one-sided 调整和 m 个检验的 α 拆分)
4. **Decision rule**: "拒绝 H0 if |t*|>c";下结论

### 题型 B: General linear F-test (interaction / nested)
1. 写 $H_0:L\beta=0$
2. 用 partial F or var-cov matrix 计算 $F^*$
3. 查 $F_{q, n-p}(\alpha)$
4. 比较、结论

### 题型 C: stepAIC 输出读懂
- 辨认候选 (SSR 最大 / SSE 最小)
- 判断能否加入 (与 AIC=<none> 比较,或算 ΔAIC)

### 题型 D: VIF 手算 + 解释
- 套公式(从 SE, σ̂, $\sum(x-\bar x)^2$)
- 判断是否严重 (>10)
- 讨论对 inference vs prediction 的影响

### 题型 E: 诊断图阅读
- Residual vs fitted: **漏斗形** → 异方差; **弯曲** → 非线性
- QQ plot: 尾部离直线 → 非正态
- Hat QQ: 右尾跳出 → high leverage 个数
- Studentized deleted QQ: 尾部极值 → outlier 个数
- Cook QQ: 极值跳出 → influential 个数

### 题型 F: 正式 outlier test with Bonferroni
- m 个 dummy, 找 critical = $t_{n-p-m}(\alpha/(2m))$
- 比较每个 |t| 与 critical

### 题型 G: WLS / Ridge / LASSO 概念题
- 何时使用、与 OLS 区别、λ 如何选、结果解释

---

## 11. 考前 48 小时 must-do checklist (中英)

- [ ] 背熟 **三类残差** 公式与用途 / Memorize **3 residual** formulas
- [ ] 能默写 **Bonferroni / Scheffé / Tukey** 的 critical value / Critical values
- [ ] 能从 `vcov` 手算任意 contrast 的 SE / Compute SE of any contrast
- [ ] 能读 `stepAIC` 输出并判断候选变量 / Read stepAIC output
- [ ] 能从 `summary(lm)` 推出 **n**, **p**, **σ̂** / Derive n, p, σ̂ from summary
- [ ] 能手算 **VIF** (两种形式) / Compute VIF two ways
- [ ] 清楚 **leverage vs outlier vs influence** 三者区别 / Distinguish 3 concepts
- [ ] 会 **Bonferroni 同时检验 m 个 outlier** / Simultaneous outlier test
- [ ] 熟悉 **WLS 两种权重估计** / 2 weight estimation methods
- [ ] **Ridge / LASSO / MCP / SCAD / PCR** 的适用场景比较 / Compare shrinkage methods
- [ ] **2021 past paper** 完整做一遍 / Complete past paper once

---

## 12. 常见英文术语对照 English Glossary

| 中文 | English |
|---|---|
| 多重比较 | multiple comparisons |
| 族错误率 | family-wise error rate (FWER) |
| 对比 | contrast |
| 同步置信区间 | simultaneous confidence interval |
| Studentized range | Studentized range distribution |
| 一般线性检验 | general linear test |
| 残差平方和 | residual sum of squares (SSE) |
| 回归平方和 | regression sum of squares (SSR / SSReg) |
| 交互 | interaction |
| 哑变量 | dummy variable |
| 方差膨胀因子 | variance inflation factor (VIF) |
| 多重共线性 | multicollinearity |
| 杠杆值 | leverage / hat value |
| 帽子矩阵 | hat matrix |
| 标准化残差 | standardized residual (internally studentized) |
| 学生化删除残差 | studentized deleted residual (externally studentized) |
| Cook 距离 | Cook's distance |
| 影响观测 | influential observation |
| 异常值 | outlier |
| 前向选择 | forward selection |
| 后向剔除 | backward elimination |
| 逐步回归 | stepwise regression |
| 加权最小二乘 | weighted least squares (WLS) |
| 岭回归 | ridge regression |
| 套索 | LASSO |
| 主成分回归 | principal component regression (PCR) |
| 协方差分析 | ANCOVA (analysis of covariance) |
| 广义交叉验证 | generalized cross-validation (GCV) |
| 对数似然 | log-likelihood |
| 偏差 | deviance |
| Box-Cox 变换 | Box-Cox transformation |

---

## 附录: Studentized Range & Dunnett 查表口诀

- **Tukey critical**: 输入 (q, df) → 输出 $q_{q,df}(\alpha)$ → 除以 $\sqrt 2$ 得 HSD critical for $t$-scale。
- **Past paper 例**: $q_{3,12}(0.05)=3.773 \Rightarrow c_{\text{Tukey}}=3.773/\sqrt 2 = 2.668$。
- **Bonferroni vs Tukey**: 若 m=3 pairs, $t_{12}(0.05/6)=t_{12}(0.00833)\approx 2.779$;Tukey 2.668 < Bonferroni 2.779 → Tukey 更优。

---

**祝考试顺利!Good luck!** 🎯
