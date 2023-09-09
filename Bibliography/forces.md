# Forces
Let

$$
\bar{b} _i = \bar{r} _{i+1} - \bar{r} _i \quad , \quad d _i = \left\|\bar{b} _i\right\| \quad , \quad \cos(\theta _i) = \frac{\bar{b} _i \cdot \bar{b} _{i-1}}{d _i d _{i-1}} \quad , \quad d _{ij} = \left\|\bar{r} _i - \bar{r} _j\right\|
$$

## Elastic potential

$$
V = \sum _{i=0} ^{N-2} \frac{1}{2} k _e (d _i - l _0) ^2
$$

$$
-\frac{\partial V}{\partial \bar{r} _i} = k _e ((1 - l _0 / d _i) \bar{b} _i - (1 - l _0 / d _{i-1}) \bar{b} _{i-1})
$$

## Bending potential

$$
V = \sum _{i=1} ^{N-2} k _b (1 - \cos(\theta _i))
$$

$$
-\frac{\partial V}{\partial \bar{r} _i} = k _b \frac{\partial}{\partial \bar{r} _i}\left(\cos(\theta _{i-1}) + \cos(\theta _i) + \cos(\theta _{i+1})\right)
$$

$$
\frac{\partial}{\partial \bar{r} _i}\cos(\theta _{i-1}) = \frac{\bar{b} _{i-2}}{d _{i-1} d _{i-2}} - \cos(\theta _{i-1}) \bar{b} _{i-1} / d _{i-1} ^2
$$

$$
\frac{\partial}{\partial \bar{r} _i}\cos(\theta _i) = \frac{\bar{b} _i - \bar{b} _{i-1}}{d _id _{i-1}} + \cos(\theta _i)\left(\bar{b} _i / d _i ^2 - \bar{b} _{i-1} / d _{i-1} ^2\right)
$$

$$
\frac{\partial}{\partial \bar{r} _i}\cos(\theta _{i+1}) = \frac{-\bar{b} _{i+1}}{d _{i+1} d _{i}} + \cos(\theta _{i+1})\bar{b} _{i} / d _{i} ^2
$$

## Lennard-Jones potential

$$
V = \sum _{i=0} ^{N-1} \sum _{j=i+1} ^{N-1} 4 \varepsilon \left[\left(\frac{\sigma}{d _{ij}}\right) ^{12} - \left(\frac{\sigma}{d _{ij}}\right) ^{6}\right]
$$

$$
-\frac{\partial V}{\partial \bar{r} _i} = \sum _{j \neq i} 4 \varepsilon \left[12 \frac{\sigma ^{12}}{d _{ij} ^{14}} - 6 \frac{\sigma ^{6}}{d _{ij} ^{8}}\right] (\bar{r} _i - \bar{r} _j)
$$
