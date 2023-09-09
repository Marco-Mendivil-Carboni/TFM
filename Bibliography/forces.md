# Forces

$$
\bar{b}_i = \bar{r}_{i+1}-\bar{r}_i \quad , \quad d_i = \left\|\bar{b}_i\right\| \quad , \quad \cos(\theta_i) = \frac{\bar{b}_i\cdot\bar{b}_{i-1}}{d_id_{i-1}} \quad , \quad d_{ij} = \left\|\bar{r}_i-\bar{r}_j\right\|
$$

## Elastic potential

$$
V=\sum_{i=0}^{N-2}\frac{1}{2}k_e(d_i-l_0)^2
$$

$$
-\frac{\partial V}{\partial \bar{r}_i} = k_e((1-l_0/d_i)\bar{b}_i-(1-l_0/d_{i-1})\bar{b}_{i-1})
$$

## Bending potential

$$
V=\sum_{i=1}^{N-2} k_b (1-\cos(\theta_i))
$$

$$
-\frac{\partial V}{\partial \bar{r}_i} = k_b \frac{\partial}{\partial \bar{r}_i}\left(\cos(\theta_{i-1})+\cos(\theta_i)+\cos(\theta_{i+1})\right)
$$
$$
\frac{\partial}{\partial \bar{r}_i}\cos(\theta_{i-1}) = \frac{\bar{b}_{i-2}}{d_{i-1}d_{i-2}}-\cos(\theta_{i-1})\bar{b}_{i-1}/d_{i-1}^2
$$
$$
\frac{\partial}{\partial \bar{r}_i}\cos(\theta_i) = \frac{\bar{b}_i-\bar{b}_{i-1}}{d_id_{i-1}}+\cos(\theta_i)\left(\bar{b}_i/d_i^2-\bar{b}_{i-1}/d_{i-1}^2\right)
$$
$$
\frac{\partial}{\partial \bar{r}_i}\cos(\theta_{i+1}) = \frac{-\bar{b}_{i+1}}{d_{i+1}d_{i}}+\cos(\theta_{i+1})\bar{b}_{i}/d_{i}^2
$$

## Lennard-Jones potential

$$
V=\sum_{i=0}^{N-1}\sum_{j=i+1}^{N-1}4\varepsilon\left[\left(\frac{\sigma}{d_{ij}}\right)^{12}-\left(\frac{\sigma}{d_{ij}}\right)^{6}\right]
$$

$$
-\frac{\partial V}{\partial \bar{r}_i} = \sum_{j \neq i}4\varepsilon\left[ 12\frac{\sigma^{12}}{d_{ij}^{14}}-6\frac{\sigma^{6}}{d_{ij}^{8}}\right] (\bar{r}_i-\bar{r}_j)
$$
