
# Forces

Let

$$
\bar{b} _i = \bar{r} _{i+1} - \bar{r} _i \ , \ d _i = \left\|\bar{b} _i\right\| \ , \ \cos(\theta _i) = \frac{\bar{b} _i \cdot \bar{b} _{i-1}}{d _i d _{i-1}} \ , \ d _{ij} = \left\|\bar{r} _i - \bar{r} _j\right\| \ , \ {d'} _{ij} = \left\|\bar{r} _i - {\bar{r}'} _j\right\|
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

## Wang-Frenkel potential

$$
V = \sum _{i=0} ^{N-1} \sum _{j=i+1} ^{N-1} \varepsilon \left(\frac{1}{d _{ij} ^{2}} - 1\right)\left(\frac{4}{d _{ij} ^{2}} - 1\right) ^{2}
$$

$$
-\frac{\partial V}{\partial \bar{r} _i} = \sum _{j \neq i} \varepsilon \left(\frac{18 d _{ij} ^{4} - 96 d _{ij} ^{2} + 96}{d _{ij} ^{8}}\right) (\bar{r} _i - \bar{r} _j)
$$

$$
d _m = \frac{2}{\sqrt{3}}
$$

## Soft-Repulsive potential

$$
V = \sum _{i=0} ^{N-1} \sum _{j=i+1} ^{N-1} 128 \left(d _{ij} + d _m/2\right) \left(d _{ij} - d _m\right) ^{2}
$$

$$
-\frac{\partial V}{\partial \bar{r} _i} = \sum _{j \neq i} 128 \left(3 d _{ij} - 3 d _m\right) (\bar{r} _i - \bar{r} _j)
$$

## Binding potential

$$
V = \sum _{i=0} ^{N-1} \sum _{j=0} ^{n _l - 1} \varepsilon \left(\frac{4}{3} \frac{{d'} _{ij} ^{2}}{r _c ^{2}} - \frac{1}{3} \frac{{d'} _{ij} ^{8}}{r _c ^{8}} - 1\right)
$$

$$
-\frac{\partial V}{\partial \bar{r} _i} = \sum _{j} \varepsilon \frac{8}{3 r _c ^{2}} \left( \frac{{d'} _{ij} ^{6}}{r _c ^{6}} - 1\right) (\bar{r} _i - {\bar{r}'} _j)
$$
