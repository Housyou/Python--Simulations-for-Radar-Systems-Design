**解：**  
平滑后的电压表达式为
> when $t\lt 0$

$$v(t)=K+\frac{K+1}{\tau}t\qquad $$

> when $t\ge 0$

$$v(t)=K-\frac{K+1}{\tau}t\qquad $$

由于$P=\frac{U^2}{R}\propto U^2$，且$P_{real}=L_s\cdot P$，故

$$L_s=\frac{U_{real}^2}{U_0^2}=\frac{v^2}{K^2}$$

在间隔$\{0,\tau\}$上时

$$L_s=1-2\frac{K+1}{K\tau}t+(\frac{K+1}{K\tau})^2t^2$$

实际需考虑的时间为$\{-\frac{\tau}{2},\frac{\tau}{2}\}$，由于对称性只需计算正半部分，平均功率损失为

$$\overline L_s=\frac{1}{\frac{\tau}{2}}\int_0^\frac{\tau}{2}(1-2\frac{K+1}{K\tau}t+(\frac{K+1}{K\tau})^2t^2)dt$$

即

$$\overline L_s=\frac{2}{\tau}(t-\frac{K+1}{K\tau}t^2+(\frac{K+1}{K\tau})^2\frac{t^3}{3})\bigg|_0^\frac{\tau}{2}$$

即

$$\overline L_s=1-\frac{K+1}{2K}+\frac{(K+1)^2}{12K^2}$$