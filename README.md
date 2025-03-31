# 1D-and-2D-WENO
给出了1维和2维双曲守恒律方程的WENO差分格式
WENO格式中带有迎风Lax-Friedrishs分解
对时间变量采用了三阶TVD-RK
只需要把对应的方程和初值改成需要的即可，代码中以线性方程和Burgers方程为例

主要参考Jiang G S, Shu C W. Efficient implementation of weighted ENO schemes[J]. Journal of Computational Physics, 1996, 126(1): 202-228.
