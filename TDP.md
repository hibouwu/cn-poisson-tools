# Méthodes directes et itératives pour l'équation de la chaleur 1D / 一维稳态热方程的直接与迭代法

T. Dufaud [thomas.dufaud@uvsq.fr](thomas.dufaud@uvsq.fr), J. Gurhem [jgurhem@aneo.fr](jgurhem@aneo.fr) — M1 CHPS —  
作者与课程信息

## 1. Résolution de l'équation de la chaleur 1D stationnaire / 一维稳态热方程的求解

Objectif : appliquer les algorithmes vus en cours/TD pour résoudre un système linéaire issu de la discrétisation par différences finies de l'équation de la chaleur 1D stationnaire. Implémentations en C avec BLAS et LAPACK. Un code à trous est fourni. Validation possible via les implémentations Scilab réalisées en TD.  
目标：将课上/TD 的算法用于求解由一维稳态热方程差分离散得到的线性系统；用 C + BLAS/LAPACK 实现，提供模板代码，可用 TD 的 Scilab 实现做验证。

Demande : analyse critique des résultats (complexités temps/espace). Pour chaque algorithme, mesurer le temps d'exécution selon la taille de la matrice et tracer des courbes de performance. Un dépôt doit être fourni avant la troisième séance (correction partie 3). Un rapport est attendu après les vacances de Noël sur les parties indiquées.  
要求：对时间/空间复杂度进行分析；各算法测量随矩阵规模的运行时间并绘制性能曲线。第三次课前提交仓库（含第 3 部分修正），圣诞假期后提交指定部分的报告。

Environnement : Langage C, BLAS/LAPACK, Gnuplot.  
环境：C 语言，BLAS/LAPACK，Gnuplot。

## 2. Plan des séances / 课程安排

- Séance 1 : rappels Docker ; mise en place env. C + BLAS/LAPACK avec Docker ; travail préliminaire : cas de test ; DS3.  
  第 1 次课：Docker 回顾；用 Docker 搭建 C+BLAS/LAPACK 环境；预备工作：测试用例；DS3。  
- Séance 2 : méthode directe et stockage bande ; DS4.  
  第 2 次课：直接法与带状存储；DS4。  
- Séance 3 : rendu de la partie 1 (commit/push) avant la séance ; méthodes de résolution itérative : Richardson, Jacobi, Gauss-Seidel.  
  第 3 次课：课前提交第 1 部分（commit/push）；迭代法：Richardson、Jacobi、Gauss-Seidel。

## 3. Travail préliminaire : cas de test / 前期工作：测试用例

Équation de la chaleur 1D stationnaire dans un milieu immobile, linéaire, homogène, isotrope, avec terme source :  
一维稳态热方程（静止、线性、均匀、各向同性，含源项）：

$$
\left\{
\begin{array}{l}
-k \frac{\partial^{2} T}{\partial x^{2}} = g, \quad x \in ]0, 1[ \\
T(0) = T_0 \\
T(1) = T_1
\end{array}
\right.
$$

où $g$ est le terme source, $k>0$ la conductivité, et $T_0 < T_1$ les températures aux bords.  
其中 $g$ 为源项，$k>0$ 为导热系数，$T_0<T_1$ 为边界温度。  
Discrétisation par schéma centré d'ordre 2 sur $n+2$ noeuds $x_i$, $i=0,\dots,n+1$, pas $h$ constant.  
用二阶中心差分离散，在 $n+2$ 个节点 $x_i$（$i=0,\dots,n+1$）上，步长 $h$ 恒定。

Équation discrète en chaque noeud :  
每个节点的离散方程：
$$
-k \left(\frac{\partial^{2} T}{\partial x^{2}}\right)_i = g_i
$$

Système global :  
整体系统：
$$
A u = f, \quad A \in \mathbb{R}^{n \times n}, \; u,f \in \mathbb{R}^n
$$

Cas $g=0$ (sans source). Solution analytique :  
无源项（$g=0$）时的解析解：
$$
T(x) = T_0 + x (T_1 - T_0)
$$

### ▷ Exercice 1 (fait en TD) / 练习 1（已在 TD 完成）

1. Approximer $\partial^2 T / \partial x^2$ par un schéma centré d'ordre 2.  
   用二阶中心差分近似二阶导数。  
2. Écrire le système linéaire de dimension $n$ correspondant au problème.  
   写出对应的 $n$ 维线性系统。  
3. Présenter le problème et sa discrétisation dans l'introduction du rapport.  
   在报告引言中介绍问题与离散过程。

Objectif : développer un code C utilisant BLAS et LAPACK pour résoudre efficacement.  
目标：用 BLAS/LAPACK 的 C 代码高效求解。

## 4. Exercice 2 : mise en place de l'environnement (C + BLAS/LAPACK) / 练习 2：环境搭建

1. Créer un projet versionné avec git et un dépôt public.  
   创建 git 项目并建公开仓库。  
2. Vérifier l'installation des bibliothèques CBLAS et CLAPACK (`/usr/local/lib`, `/usr/local/include`). Sous Ubuntu : `apt-get install libblas-dev liblapack-dev`. Le code fourni peut aussi être compilé via Docker (`docker/Dockerfile`) avec : `docker build -f docker/Dockerfile --progress plain -t tp-cn:latest .`  
   检查 CBLAS/CLAPACK 安装（路径 `/usr/local/lib`, `/usr/local/include`）。Ubuntu 可用 `apt-get install libblas-dev liblapack-dev`。也可用 Dockerfile 构建：`docker build -f docker/Dockerfile --progress plain -t tp-cn:latest .`  
3. Écrire un Makefile et un code de test (ex. `testenv.c`) pour valider l'environnement.  
   编写 Makefile 和测试代码（如 `testenv.c`）验证环境。

Conseil : créer un fichier `lib_poisson1D.c` regroupant toutes les fonctions pour le problème de Poisson 1D et des fichiers séparés pour les codes de résolution.  
建议：创建 `lib_poisson1D.c` 汇总一维泊松问题函数，解法代码分文件存放。

## 5. Méthode directe et stockage bande / 直接法与带状存储

On considère des matrices tridiagonales ; on étudie le stockage par bande (General Band).  
处理三对角矩阵，研究带状存储（GB）。

Une matrice $m \times n$ avec $kl$ sous-diagonales et $ku$ sur-diagonales se stocke dans un tableau $(kl+ku+1) \times n$ (colonnes de la matrice dans les colonnes du tableau, diagonales dans les lignes). Utilisation pertinente si $kl, ku \ll \min(m,n)$. Dans LAPACK, les matrices de ce type ont un nom se terminant par `B`. L'élément $a_{ij}$ est stocké en $AB(ku+1+i-j,\, j)$.  
含 $kl$ 条下对角、$ku$ 条上对角的 $m \times n$ 矩阵可存入 $(kl+ku+1) \times n$ 的二维数组：矩阵列对应数组列，对角线对应数组行；当 $kl, ku \ll \min(m,n)$ 时适用。LAPACK 中此类矩阵名以 `B` 结尾，元素 $a_{ij}$ 存于 $AB(ku+1+i-j,\, j)$。

Exemple ($m=n=5$, $kl=2$, $ku=1$) :  
例子（$m=n=5$, $kl=2$, $ku=1$）：
$$
A = \begin{pmatrix}
a_{11} & a_{12} & & & \\
a_{21} & a_{22} & a_{23} & & \\
a_{31} & a_{32} & a_{33} & a_{34} & \\
& a_{42} & a_{43} & a_{44} & a_{45} \\
& & a_{53} & a_{54} & a_{55}
\end{pmatrix}
$$

Stockage GB :  
GB 存储：
$$
AB = \begin{pmatrix}
* & a_{12} & a_{23} & a_{34} & a_{45} \\
a_{11} & a_{22} & a_{33} & a_{44} & a_{55} \\
a_{21} & a_{32} & a_{43} & a_{54} & * \\
a_{31} & a_{42} & a_{53} & * & *
\end{pmatrix}
$$

Les éléments `*` sont généralement remplis à 0.  
`*` 处通常填 0。

## ▷ Exercice 3. Référence et utilisation de BLAS/LAPACK / 练习 3：BLAS/LAPACK 参考

1. En C, comment déclarer et allouer une matrice pour BLAS/LAPACK ?  
   在 C 中如何声明和分配矩阵以供 BLAS/LAPACK 使用？  
2. Signification de la constante `LAPACK_COL_MAJOR` ?  
   `LAPACK_COL_MAJOR` 常量含义？  
3. À quoi correspond la leading dimension (`ld`) ?  
   主维度 (`ld`) 指什么？  
4. Que fait la fonction `dgbmv` ? Quelle méthode implémente-t-elle ?  
   `dgbmv` 做什么，实现哪种运算/方法？  
5. Que fait la fonction `dgbtrf` ?  
   `dgbtrf` 功能？  
6. Que fait la fonction `dgbtrs` ?  
   `dgbtrs` 功能？  
7. Que fait la fonction `dgbsv` ?  
   `dgbsv` 功能？  
8. Comment calculer la norme du résidu relatif avec des appels BLAS ?  
   如何用 BLAS 调用计算相对残差范数？

## ▷ Exercice 4. Stockage GB et appel à `dgbmv` / 练习 4：GB 存储与 `dgbmv`

1. Écrire le stockage GB en priorité colonne pour la matrice de Poisson 1D.  
   写出一维泊松矩阵的列优先 GB 存储。  
2. Utiliser la fonction BLAS `dgbmv` avec cette matrice.  
   用该矩阵调用 `dgbmv`。  
3. Proposer une méthode de validation.  
   给出验证方法。

## ▷ Exercice 5. `dgbtrf`, `dgbtrs`, `dgbsv` (à rendre) / 练习 5（需提交）

1. Résoudre le système linéaire par une méthode directe avec LAPACK.  
   用 LAPACK 直接法求解线性系统。  
2. Évaluer les performances. Discuter la complexité des méthodes appelées.  
   评估性能并讨论所用方法的复杂度。

## ▷ Exercice 6. LU pour matrices tridiagonales (à rendre) / 练习 6（需提交）

1. Implémenter la factorisation LU pour matrices tridiagonales au format GB.  
   在 GB 格式下实现三对角矩阵的 LU 分解。  
2. Proposer une méthode de validation.  
   提出验证方法。  
3. Évaluer les performances et comparer. Discuter les complexités.  
   评估并比较性能，讨论复杂度。

## 6. Méthodes itératives / 迭代方法

On souhaite résoudre l'équation de la chaleur par méthodes itératives. Étude théorique puis maquettes Scilab, puis code C avec BLAS.  
用迭代法求解热方程：先理论与 Scilab 验证，再用 C + BLAS 编码。

### ▷ Exercice 7. Implémentation C — Richardson / 练习 7：C 实现 Richardson

1. Implémenter Richardson en C avec matrices au format GB ; sauvegarder le résidu dans un vecteur.  
   用 C/GB 格式实现 Richardson，并保存残差向量。  
2. Calculer l'erreur vs solution analytique.  
   计算相对解析解的误差。  
3. Analyser la convergence (tracer l'historique).  
   分析收敛性（绘制残差/误差历史）。

### ▷ Exercice 8. Implémentation C — Jacobi / 练习 8：C 实现 Jacobi

1. Écrire une fonction Jacobi pour matrice tridiagonale GB.  
   为 GB 三对角矩阵实现 Jacobi。  
2. Calculer l'erreur vs solution analytique.  
   计算相对解析解的误差。  
3. Analyser la convergence (historique).  
   分析收敛性（历史曲线）。

### ▷ Exercice 9. Implémentation C — Gauss-Seidel / 练习 9：C 实现 Gauss-Seidel

1. Écrire une fonction Gauss-Seidel pour matrice tridiagonale GB ; discuter le calcul de $(D - E)^{-1} r^k$.  
   为 GB 三对角矩阵实现 Gauss-Seidel；讨论 $(D - E)^{-1} r^k$ 的计算。  
2. Calculer l'erreur vs solution analytique.  
   计算相对解析解的误差。  
3. Analyser la convergence (historique).  
   分析收敛性（历史曲线）。

## 7. Autres formats de stockage / 其他存储格式

Après LU/Richardson/Jacobi/Gauss-Seidel en GB, compléter pour CSR ou CSC.  
在 GB 完成 LU、Richardson、Jacobi、Gauss-Seidel 后，扩展到 CSR/CSC。

### ▷ Exercice 10. Formats CSR / CSC / 练习 10：CSR/CSC

1. Écrire le stockage CSR pour Poisson 1D.  
   写出一维泊松的 CSR 存储。  
2. Écrire le stockage CSC pour Poisson 1D.  
   写出一维泊松的 CSC 存储。  
3. Écrire `dcsrmv` et `dcsmv` (produit matrice-vecteur en CSR/CSC).  
   编写 `dcsrmv` 与 `dcsmv`（CSR/CSC 矩阵乘向量）。  
4. Adapter les algorithmes précédents à ces formats.  
   将前述算法适配到这些格式。
