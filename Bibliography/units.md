
# Units

## Parameters and constants

+ According to Robinson et al. there are 11 nucleosomes for each 11 nm of heterochromatin fiber, which is 33 nm wide.

+ We will consider that each one of our particles has 33 nucleosomes and thus a diameter of 33 nm.

+ Although in the euchromatic state nucleosomes are not organized in a well-defined and stable manner they interact strongly in various ways forming the so called TADs (Topologically Associating Domains). Inside this domains the nuclesomes are presumably very close together so assuming that they ocuppy a similar volume as in the heterochromatic state is justified.

+ As Bajpai et al. mention, we can more or less associate LADs to heterochromatin and non-LADs to euchromatin so in our model we will use only two types of particles for simplicity:
    + LADh (A): LAD heterochromatin
    + LNDe (B): non-LAD euchromatin

+ Despite what one may think at first according to Bajpai at al. the exact sequence of particles used in this kind of models does not change significantly the results, as long as the statistical properties of the sequences are similar. Nonetheless, we will base the sequence of our coarse grained model on experimental data. We will use the lamina-associated domains in the C. elegans genome reported by Ho et al. to determine the LADh particles of our polymer. There are several good reasons to use this organism's sequence in particular: it is a model organism, its genome is specially short and its chromosomes do not have pericentromeric constitutive heterochromatin, which heavily encourages the use of models like the one of Falk et al. with three types of beads.

+ We will also separate the polymer at the end of each one of the 6 chromosomes, whose lengths we have obtained from the [GEO NCBI database](https://www.ncbi.nlm.nih.gov/gdv/browser/genome/?id=GCF_000002985.6).

+ Because C. elegans cells (like most eukaryotes) are diploids inside the nucleus there will be two copies of each chromosome, making a total of 12 chromosomes.

+ According to Bystricky et al. heterochromatin has a persistence length of $\sim$ 200 nm which in our model translates to $k_b=6$. However, since heterochromatin only represents $\sim$ 8% of the whole sequence and belongs in most cases to some of the LAD domains (which represent almost half of the sequence) we will average the persistence length inside the LAD domains and make LADh particle have $k_b=1$.

+ While euchromatin has virtually no persistence length at this scale so LNDe particles will have $k_b=0$.

+ As explained by Camara the compaction and rigidity of heterochromatin makes it self-attractive due to depletion attraction. Meanwhile euchromatin is much softer and does not interact with itself through this mecanism. Camara also mentions that inactive regions (that is heterochromatin like domains) have higher self-affinity. Thus we will consider that LNDe particles are purely repulsive while LADh particles are self-attractive.

+ Thanks to the wonderful work of Falk et al. we can make an educated guess of the interaction strengths between our two types of particles. In their more elaborate model with three types of particles (one for euchromatin and two for heterochromatin) they found that euchromatin interacted very little with itself and the other types of particles so we can safely make it purely repulsive for simplicity. The interaction strength between heterochromatin particles was different for the two types of heterochromatin particles but their average was of the order of $0.5k_BT$ which is the value we will use.

+ Chromatin has other distinguishable regions like centromeres, telomeres and NORs (nucleolus organizer regions) which we will not consider since they consitute smaller fractions of the genome than general euchromatin and heterochromatin. For example, according to [BioNumbers](https://bionumbers.hms.harvard.edu/bionumber.aspx?id=113849&ver=3&trm=nucleolus&org=) the nucleolus has an average radius of 0.9 $\text{um}$ and the nucleus of 2.7 $\text{um}$ so NORs occupy less than 4% of the nucleus volume despite their crucial biological function.

+ The exact value of the elastic constant is not really relevant in this context as long as it is high enough to keep the distance between bonded particles from fluctuating too much. According to Falo et al. the statiscal properties of an elastic freely jointed chain only depend on the elastic constant through the fraction $2/(1+\beta k_el_0^2)$ so taking $k_e=128e_u/l_u^2$ should suffice.

+ In a similar but simpler way to Bajpai et al. we will add potential wells throughout the lamina (at random positions) which we will call lamina binding sites. This sites represent the places of the lamina to which a LAD (a LADh particle in our model) can bind via the so called lamina binding proteins. We will make this wells as wide as our particles and $8k_BT$ deep since according to Maji et al. using a very similar model the lamina-chromatin interaction strengths that best reproduced experimental observations were above $2k_BT$ and Bajpai et al. used a $10k_BT$ deep potential.

+ With all this in mind we have a very idealized but relatively complex model of the chromatin fiber confined inside a spherical nucleus. By modulating the unknown parameters of this model we wish to investigate the different conformations that chromatin can adopt and hopefully reproduce the features seen on experiments.

+ The modifiable parameters of the model will be:
    + The nucleus radius and bleb geometry ($R_n$, $R_o$, $R_b$).
    + The number of lamina binding sites ($N'$).

+ For our first simulation we will give this parameters some reasonable values under normal conditions. Ikegami et al. state in their article that C. elegans nuclei have a radius of approximately 1 $\text{um}$.

| name            | value (PU)            | value (SI)                       |
|-----------------|-----------------------|----------------------------------|
| $\xi$           | $1m_u/t_u$            | $3\pi\eta d=0.28\text{ pg / us}$ |
| $k_B$           | $(1/298)e_u/\text{K}$ | $0.014\text{ pN nm / K}$         |
| $l_0$           | $1l_u$                | $33\text{ nm}$                   |
| $k_e$           | $128e_u/l_u^2$        | $0.48\text{ pN / nm}$            |
| $k_b^A$         | $1e_u$                | $4.11\text{ pN nm}$              |
| $k_b^B$         | $0e_u$                | $0\text{ pN nm}$                 |
| $\sigma$        | $1l_u$                | $33\text{ nm}$                   |
| $\epsilon^{AA}$ | $0.5e_u$              | $2.06\text{ pN nm}$              |
| $\epsilon^{AB}$ | $0e_u$                | $0\text{ pN nm}$                 |
| $\epsilon^{BB}$ | $0e_u$                | $0\text{ pN nm}$                 |
| $\epsilon'$     | $8e_u$                | $32.9\text{ pN nm}$              |
| $dt$            | $(1/2048)t_u$         | $36\text{ ns}$                   |

+ We take the viscosity of water at 25ºC to be $\eta=0.0009\text{ Pa s}$ and the hydrodynamic radius to be $d=\sigma$.

+ The Reynolds number is $\text{Re}=\rho v_ul_u/\eta=1.7\cdot 10^{-5}$.

## Program units

| quantity    | unit                            |
|-------------|---------------------------------|
| length      | $\sigma=33\text{ nm}$           |
| energy      | $k_BT=4.11\text{ pN nm}$        |
| time        | $\xi\sigma^2/k_BT=74\text{ us}$ |
| force       | $e_u/l_u=0.12\text{ pN}$        |
| velocity    | $l_u/t_u=0.46\text{ nm / us}$   |
| temperature | $1\text{ K}$                    |
