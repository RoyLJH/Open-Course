## 矩阵分析第三次作业

### 尼博琳 201928014628013

#### 11

（a）
$$
\begin{align}
\mathbf{B}
&=\begin{pmatrix} 2 & 0 & -1 \\ -1 & 1 & 1 \\ -1 & 2 & 1 \end{pmatrix} \\
&=\begin{pmatrix} 2 & 0 & -1 \\ -1 & 1 & 1 \\ -1 & 0 & 1 \end{pmatrix}+ \begin{pmatrix}0 & 0 & 0 \\0 & 0 & 0\\0&2&0 \end{pmatrix} \\
&=\begin{pmatrix}2 & 0 & -1 \\ -1 & 1 & 1 \\ -1 & 0 & 1 \end{pmatrix}+2\begin{pmatrix}0 \\ 0 \\ 1 \end{pmatrix} \begin{pmatrix}0 & 1 & 0 \end{pmatrix} \\ 
&= \mathbf{A} + 2\mathbf{e_3}\mathbf{e_2}^{T}
\end{align}
$$
根据 Sherman-Morrison 定理可得
$$
\begin{align}
\mathbf{B}^{-1}&=\left(\mathbf{A}+2 \mathbf{e}_{3} \mathbf{e}_{2}^{T}\right)^{-1}=\mathbf{A}^{-1}-2 \frac{\mathbf{A}^{-1} \mathbf{e}_{3} \mathbf{e}_{2}^{T} \mathbf{A}^{-1}}{1+2 \mathbf{e}_{2}^{\mathbf{T}} \mathbf{A}^{-1} \mathbf{e}_{3}} \\
&=\mathbf{A}^{-1}-2 \frac{\left[\mathbf{A}^{-1}\right]_{* 3}\left[\mathbf{A}^{-1}\right]_{2 *}}{1+2\left[\mathbf{A}^{-1}\right]_{2 3}} \\
&=\begin{pmatrix}{1} & {0} & {1} \\ {0} & {1} & {-1} \\ {1} & {0} & {2}\end{pmatrix} - 2 \frac{\begin{pmatrix}1 \\ -1 \\ 2\end{pmatrix}\begin{pmatrix} 0 & 1 & -1\end{pmatrix}}{1+2 \times (-1)} \\ 
&= \begin{pmatrix}1 & 2 & -1 \\0 & -1 & 1 \\ 1 & 4 & -2 \end{pmatrix}
\end{align}
$$
（b）

令 $\mathbf{c} = \begin{pmatrix} 0 \\ 0 \\ 1 \end{pmatrix}, d = \begin{pmatrix} 0 \\ 2 \\ 1\end{pmatrix}$

由 Sherman-Morrison 定理可得
$$
\begin{align}
\mathbf{C}^{-1} &= (\mathbf{A} + \mathbf{cd}^T)^{-1} = \mathbf{A}^{-1} - \frac{\mathbf{A}^{-1} \mathbf{cd}^T\mathbf{A}^{-1} }{1+\mathbf{d}^T\mathbf{A}^{-1}\mathbf{c}} \\
&=\begin{pmatrix}1 & 0 & 1 \\ 0 & 1 & -1\\1 & 0 & 2 \end{pmatrix} - \begin{pmatrix} 1 & 2 & 0 \\ -1 & -2 & 0 \\2 & 4 & 0\end{pmatrix} \\
&= \begin{pmatrix} 0 & -2 & 1 \\
1 & 3 & -1 \\ -1 & -4 & 2\end{pmatrix}
\end{align}
$$

#### 12

（a）
$$
\mathbf{L} = \begin{pmatrix}1 & 0 & 0 \\4 & 1 & 0 \\3 & 2 & 1 \end{pmatrix}, \mathbf{U} = \begin{pmatrix}1 & 4 & 5 \\0 & 2 & 6 \\0 & 0 & 3 \end{pmatrix}
$$
（b）

令 $\mathbf{y} = \mathbf{U} \mathbf{x_1}$，则根据 $\mathbf{Ly_1} = \mathbf{b_1}$ ,解得 $\mathbf{y_1} =\begin{pmatrix} 6 \\ -24 \\ 24\end{pmatrix}$，进而解得 $\mathbf{x_1} = \begin{pmatrix}110 \\-36\\8 \end{pmatrix}$

令 $\mathbf{y_2} = \mathbf{Ux_2}$，则根据 $\mathbf{Ly_2} = \mathbf{b_2}$,解得 $\mathbf{y_2} = \begin{pmatrix} 6 \\ -18 \\ 30\end{pmatrix}$，进而解得 $\mathbf{x_2} = \begin{pmatrix} 112 \\ -39 \\10\end{pmatrix}$

