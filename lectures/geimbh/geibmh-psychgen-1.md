# Genetics of common complex psychiatric disorders
Mark Adams

Mark Adams  
Division of Psychiatry  
`mark.adams@ed.ac.uk`  
*Genetics and Environmental Influences on Behaviour and Mental Health*

# Topics:

- **Quantitative genetics and heritability**
- **Candidate gene studies**
- **Genome-wide studies**
- **Prediction**
- **Causality**

# What is a *“common”*, *“complex”* psychiatric disorder?

**Common:** Affects 1% or more of the population  
**Complex:** Inheritance cannot be explained by a single gene

<div class="notes">

Psychiatric disorders are defined by disruption to higher-order brain
functions of moods, perceptions, thoughts, beliefs, and behaviours but
usually in the absence of major neurological impairments (consciousness,
senses, memory). Psychiatric disorders include depressive and anxiety
disorders (major depressive disorder, panic disorder), manic and
psychotic disorders (bipolar disorder, schizophrenia),
obsessive-compulsive disorders, eating disorders (anorexia nervosa,
bulimia nervosa), substance-use disorders and personality disorders.
Childhood conditions like attention-deficit/hyperactivity and autism can
also be included, but only when they lead to clinically-salient
impairment or distress. There are also many shades of sadness,
hallucinations, eccentricities, mood swings, body-image preoccupations,
recreational substances use, personalities, etc that are not psychiatric
disorders but may still be informative to study from an aetiological and
genetic standpoint.

- Sullivan PF and Geschwind DH (2019) [Defining the Genetic, Genomic,
  Cellular, and Diagnostic Architectures of Psychiatric
  Disorders](https://doi.org/10.1016/j.cell.2019.01.015). *Cell*
  doi:10.1016/j.cell.2019.01.015

</div>

## 🧬👪🚬💢🏡💞🩻🏫

- Depression: 3% in a week
- Schizophrenia: 1% in lifetime
- Bipolar disorder: 2% in lifetime
- Anxiety disorder: 6% in a week

<div class="notes">

Psychiatric disorders have many causes, correlates, and consequences
(genetics, environment, family life, substance use, relationships)

Incidence of psychiatric disorders range from the common (depression,
anxiety) to the rare (bipolar disorder, schizophrenia).

</div>

# Why genetics?

Why use genetics to study mental health and psychiatric disorders?

- Biological understanding of genes, pathways
- Shared aetiology with other disorders
- Risk prediction
- Drug retargeting
- Causal analysis of environmental risk factors

## Genetics of categorical traits

![Diagram showing the seven “characters” observed by
Mendel](assets/Mendel_seven_characters.svg)

<div class="notes">

Gregor Mendel (1822–1884), working in what is now Czechia, discovered
the transmission of traits from parents to offspring could be explained
by the inheritance of two “elements”, which we now call alleles. Mendel
was concerned with discrete or categorical phenotypes.

[Mendel pea plant
figure](https://commons.wikimedia.org/wiki/File:Mendel_seven_characters.svg)
by Mariana Ruiz (LadyofHats) \[public domain\]

</div>

## Genetics of continuous traits

![](assets/Galton-Regression-PlateIX.png)

<div class="notes">

Separately, Francis Galton (1822–1911), was studying the inheritance of
continuous or metric phenotypes. He noticed the parents who were tall
tended to have children that were slightly shorter than themselves (and
vice versa). This was termed “regression to the mean” from which the
name of the statistical method “regression” is derived.

> “In their search for universal hereditary laws, Galton and Pearson
> were driven by the linear model and the normal distribution because
> the associated parameters had scientific meaning for them that went
> beyond mere description.” - Wachsmuth, A., Wilkinson, L., & Dallal, G.
> E. (2003). [Galton’s Bend](https://dx.doi.org/10.1198/0003130031874).
> The American Statistician, 57(3), 190–192. doi:10.1198/0003130031874

For more on Galton’s legacy, see <https://adelphigenetics.org/history/>

</div>

## Reconciling categorical + continuous genetics = quantitative genetics

![](assets/Fisher-Supposition.png)

<div class="notes">

Ronald Fisher reconciled the inheritance of continuous and categorical
phenotypes by showing that a continuous phenotype could be made from the
inheritance of a large number (dozens, hundreds, or thousands)
categorical genes. The term “variance” comes from Fisher’s discoveries.

- Fisher, R. A. (1918). [XV.—The correlation between relatives on the
  supposition of Mendelian
  inheritance](https://dx.doi.org/10.1017/S0080456800012163).
  *Transactions of the Royal Society of Edinburgh*.
  doi:10.1017/S0080456800012163
- Charlesworth and Edwards (2018). [A century of
  variance](https://doi.org/10.1111/j.1740-9713.2018.01170.x).
  *Significance* 15(4). doi:10.1111/j.1740-9713.2018.01170.x
- Bodmer et al (2021) [The outstanding scientist, R.A. Fisher: his views
  on eugenics and race](https://dx.doi.org/10.1038/s41437-020-00394-6).
  *Heredity* doi:10.1038/s41437-020-00394-6

</div>

## Polygenic traits are quantitative traits

Adding up effects from a large number of genetic effects to make a
continuous phenotype is related to the Central Limit Theorem.

``` r
require(ggplot2)
```

    Loading required package: ggplot2

``` r
require(dplyr)
```

    Loading required package: dplyr


    Attaching package: 'dplyr'

    The following objects are masked from 'package:stats':

        filter, lag

    The following objects are masked from 'package:base':

        intersect, setdiff, setequal, union

``` r
# calculate expected genotype frequency for number of increasing alleles
n_allele_freq <- function(n) {
  alleles <- seq(from = 0, to = 2*n)
  freq <- dbinom(alleles, size = 2*n, prob = 0.5)
  return(data.frame(alleles, freq, loci = n))
}

number_of_loci <- c(1, 2, 3, 10)
loci_freq <- bind_rows(lapply(number_of_loci, n_allele_freq))

ggplot(loci_freq, aes(x = alleles, y = freq)) + geom_col() +
facet_grid(
  . ~ loci,
  scales = "free_x",
  space = "free_x"
) +
scale_x_continuous("Number of phenotype increasing alleles") +
scale_y_continuous("Genotype frequency")
```

![](geibmh-psychgen-1_files/figure-commonmark/polygenic_quantitative-1.png)

<div class="notes">

> “R. A. Fisher’s 1918 paper, ‘The correlation between relatives on the
> supposition of Mendelian inheritance’, resolved the often bitter
> conflict between biometricians and Mendelians, which raged for a
> decade following the rediscovery of Mendel’s work. Fisher showed that
> a complex quantitative trait could be explained by Mendelian
> inheritance if several genes affect the trait.”Because he crossed
> true-breeding plants, Mendel’s experiments showed that a single locus
> with two alleles of equal frequency results in three genotypes (see
> the figure, part a). If the allelic effects are additive, the three
> genotypes produce three phenotypes; in the case of Mendel’s
> qualitative traits, the allelic effects showed complete dominance, so
> only two phenotypes were observed. However, assuming equal and
> additive effects, 2 genes yield 9 genotypes and 5 phenotypes (part b)
> and 3 genes yield 27 genotypes and 7 phenotypes (part c). With unequal
> and non-additive allelic effects and some environmental influence,
> three genes would result in a normal bell-shaped curve of continuous
> variation (part d). This logic assumes common alleles; rare alleles
> will skew the distribution. Genome-wide association research suggests
> that many more than three genes affect most traits, which underscores
> the expectation that polygenic traits are quantitative traits.”

- Plomin, R., Haworth, C. & Davis, O. [Common disorders are quantitative
  traits](https://dx.doi.org/10.1038/nrg2670). *Nat Rev Genet* **10**,
  872-878 (2009). doi:10.1038/nrg2670

</div>
