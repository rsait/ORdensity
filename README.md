# ORdensity

In this repository is located the R package ORdensity that implements the statistical method presented in the paper [Identification of differentially expressed genes by means of outlier detection](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2318-8) by Irigoien and Arenas, BMC Bioinformatics 2018 (\[1\]).

## Introduction

An important issue in microarray data is to select, from thousands of genes, a small number of informative differentially expressed (DE) genes that may be key elements for a disease. If each gene is analyzed individually, there is a big number of hypotheses to test and a multiple comparison correction method must be used. Consequently, the resulting cut-off value may be too small. Moreover, an important issue is the selection’s replicability of the DE genes. The package ORdensity is designed to obtain a reproducible selection of DE genes by the method presented in \[1\], which is not a gene-by-gene approach. The core function 'findDEgenes' provides three measures related to the concepts of outlier and density of false positives in a neighbourhood, which allow identify the DE genes with high classification accuracy. The first measure is an index called OR and previously introduced in \[2, 3\]; the other two measures called FP and dFP were introduced in \[1\]. Additional functions provided in this package like 'preclusteredData' and 'plot' facilitate exploring and understanding the results.  As, working with large datasets, long execution times and great computational efforts are required, parallelization strategies were used to perform the analysis in a short time.

## References

\[1\] Irigoien I, Arenas C. [Identification of differentially expressed genes by means of outlier detection](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2318-8). 19(1): 317:1–317:20 (2018)

\[2\] Arenas C, Toma C, Cormand B, Irigoien I. [Identifying extreme observations, outliers and noise in clinical and genetic data.](http://www.eurekaselect.com/142998/article) Current Bioinformatics 2017;12(2):101–17.

\[3\] Arenas C, Irigoien I, Mestres F, Toma C, Cormand B. [Extreme observations in biomedical data.](https://link.springer.com/chapter/10.1007/978-3-319-55639-0_1) In: Ainsbury EA, Calle ML, Cardis E, et al., editors. Extended Abstracts Fall 2015. Trends in Mathematics vol 7. Birkhäuser, Cham: Springer; 2017. p. 3–8.

## Installation

To install the package from this repository, just run the following code

```
library('devtools')
install_github('rsait/ORdensity')
```

This package requires the `cluster` library to be installed; otherwise it will automatically install and load it. Likewise, the `foreach` library is used for parallelization.

To start working with the package, just load it in the R enviroment with the following command

```
library('ORdensity')
```

## Example

There is an example dataframe called `simexpr` shipped with the package. This data is the result of a simulation of 100 differentially expressed genes in a pool of 1000 genes. It contains 1000 observations of 62 variables. Each row correspond to a gene and contains 62 values: DEgen, gap and the values for the gene expression in 30 positive cases and in 30 negative cases. The DEgen field value is 1 for differentially expressed genes and 0 for those which are not.

First, let us extract the samples from each experimental condition from the `simexpr` database.

```
x <- simexpr[, 3:32]
y <- simexpr[, 33:62]
EXC.1 <- as.matrix(x)
EXC.2 <- as.matrix(y)
```
To create an S4 object to perform the analysis, follow this command

```
myORdensity <- new("ORdensity", Exp_cond_1 = EXC.1, Exp_cond_2 = EXC.2)
```
By default, no parallelizing is enabled. To enable it, just run instead

```
myORdensity <- new("ORdensity", Exp_cond_1 = EXC.1, Exp_cond_2 = EXC.2, parallel = TRUE)
```
The number of processes launched is equal to the number of processors detected in the machine. That number can be changed with the `nprocs` parameter

```
myORdensity <- new("ORdensity", Exp_cond_1 = EXC.1, Exp_cond_2 = EXC.2, parallel = TRUE, nprocs = 3)
```
It is also possible to enable or disable replicability, and to pass the seed to the pseudorandom number generator. The default values are 

```
myORdensity <- new("ORdensity", Exp_cond_1 = EXC.1, Exp_cond_2 = EXC.2, replicable = TRUE, seed = 0)
```
with the function using the given seed to set the random generator. If replicable = FALSE, no seed is used.

A summary of the object can be generated with the `summary` function.

```
> summary(myORdensity)
The ORdensity method has found that the optimal clustering of the data consists of 2 clusters, computed from a maximum of 10 when the ORdensity object was created
$Cluster1
$Cluster1$numberOfGenes
[1] 83

$Cluster1$CharacteristicsCluster
           OR        FP      dFP
mean 63.58153 0.6361446 2.090888
sd   40.60754 1.0377718 3.577757

$Cluster1$genes
 [1] "Gene1"   "Gene10"  "Gene100" "Gene11"  "Gene12"  "Gene13"  "Gene14" 
 [8] "Gene15"  "Gene16"  "Gene17"  "Gene2"   "Gene20"  "Gene21"  "Gene22" 
[15] "Gene23"  "Gene24"  "Gene25"  "Gene26"  "Gene27"  "Gene28"  "Gene29" 
[22] "Gene3"   "Gene30"  "Gene32"  "Gene34"  "Gene35"  "Gene36"  "Gene37" 
[29] "Gene39"  "Gene4"   "Gene40"  "Gene42"  "Gene47"  "Gene48"  "Gene49" 
[36] "Gene5"   "Gene50"  "Gene52"  "Gene53"  "Gene54"  "Gene55"  "Gene56" 
[43] "Gene57"  "Gene58"  "Gene59"  "Gene6"   "Gene60"  "Gene61"  "Gene62" 
[50] "Gene63"  "Gene64"  "Gene65"  "Gene66"  "Gene67"  "Gene68"  "Gene69" 
[57] "Gene7"   "Gene70"  "Gene71"  "Gene72"  "Gene73"  "Gene76"  "Gene79" 
[64] "Gene8"   "Gene80"  "Gene81"  "Gene82"  "Gene83"  "Gene84"  "Gene86" 
[71] "Gene87"  "Gene88"  "Gene89"  "Gene9"   "Gene90"  "Gene91"  "Gene92" 
[78] "Gene93"  "Gene95"  "Gene96"  "Gene97"  "Gene98"  "Gene99" 


$Cluster2
$Cluster2$numberOfGenes
[1] 24

$Cluster2$CharacteristicsCluster
            OR       FP      dFP
mean 11.114245 7.504167 42.69873
sd    3.318155 1.964794 20.47935

$Cluster2$genes
 [1] "Gene104" "Gene18"  "Gene19"  "Gene277" "Gene31"  "Gene33"  "Gene38" 
 [8] "Gene399" "Gene41"  "Gene43"  "Gene44"  "Gene45"  "Gene46"  "Gene51" 
[15] "Gene598" "Gene618" "Gene651" "Gene670" "Gene74"  "Gene75"  "Gene78" 
[22] "Gene85"  "Gene94"  "Gene946"
```
The summary tells us the estimated optimal clustering of the data, and the number of genes in each cluster, along with their names. The clusters are ordered in decreasig order according to the value of the mean of the `OR` statistic. We see that the mean is higher in the first cluster (61.8986) than in the second one (10.510895), which means that the first cluster is more likely composed of true differentially expressed genes, and the second one to be composed of false positives. With more clusters, the last ones are likely false negatives.

If the researcher just wants to extract the differentially expressed genes detected by the ORdensity method, a call to `findDEgenes` will return a list with the clusters found, along with the values of the OR statistic corresponding to each gene, and an indicator showing if the gene fulfil the strong and/or relaxed selection requirements. Following [1], two types of differentially expressed gene selection can be made:

* **ORdensity strong selection:** take as differentially expressed genes those with a largeOR value and with FP and dFP equal to 0.
* **ORdensity relaxed selection:** take as differentially expressed genes those with a largeOR value and with small FP and dFP values. As a reference to look for small values the expected number of false positive neighbours is computed.

The motivation of the clustering is to distinguish those false positives that score high in OR and low in meanFP and density, but are similar to other known false positives obtained by boostrapping. The procedure is detailed in [1] and it uses the PAM cluster procedure.

After running this code

```
result <- findDEgenes(myORdensity)
```

the method indicated that the optimal clustering consists of just two clusters, 

```
The ORdensity method has found that the optimal clustering of the data consists of 2 clusters
```

we could then look the results

```
> result
$neighbours
[1] 10

$expectedFalsePositiveNeighbours
[1] 8.237232

$clusters
$clusters[[1]]
        id        OR  FP        dFP Strong Relaxed
62  Gene62 175.04323 0.0  0.0000000      S       R
50  Gene50 172.28779 0.0  0.0000000      S       R
61  Gene61 155.53626 0.0  0.0000000      S       R
70  Gene70 152.54705 0.0  0.0000000      S       R
2    Gene2 149.65335 0.0  0.0000000      S       R
68  Gene68 148.87294 0.0  0.0000000      S       R
52  Gene52 134.87709 0.0  0.0000000      S       R
34  Gene34  77.32266 0.0  0.0000000      S       R
71  Gene71  74.37631 0.0  0.0000000      S       R
53  Gene53  69.74061 0.0  0.0000000      S       R
28  Gene28  62.75227 0.0  0.0000000      S       R
80  Gene81  54.30241 0.0  0.0000000      S       R
40  Gene40 104.99403 0.1  0.1991604      -       R
81  Gene82  93.45266 0.1  0.2892853      -       R
15  Gene15  93.16311 0.1  0.2956274      -       R
64  Gene64  92.59491 0.1  0.2513645      -       R
73  Gene73  91.67733 0.1  0.2812384      -       R
67  Gene67  88.44531 0.1  0.3110092      -       R
29  Gene29  88.43950 0.1  0.3097112      -       R
6    Gene6  84.47275 0.1  0.3047120      -       R
23  Gene23  76.09034 0.1  0.2390754      -       R
90  Gene91  52.91574 0.1  0.3281189      -       R
21  Gene21  50.88096 0.1  0.3199397      -       R
96  Gene97  48.44279 0.1  0.2908925      -       R
42  Gene42  44.77756 0.1  0.5999895      -       R
1    Gene1  40.13427 0.1  0.3922824      -       R
82  Gene83  39.63939 0.1  0.4893582      -       R
14  Gene14  38.54946 0.1  0.4340459      -       R
78  Gene79  35.04591 0.1  0.3481011      -       R
32  Gene32 144.80790 0.2  0.1803355      -       R
7    Gene7 143.42201 0.2  0.1889123      -       R
10  Gene10 130.61844 0.2  0.2217622      -       R
36  Gene36 120.82917 0.2  0.2507271      -       R
92  Gene93 103.66296 0.2  0.3341376      -       R
65  Gene65  98.36185 0.2  0.3128623      -       R
24  Gene24  92.20146 0.2  0.3046823      -       R
17  Gene17  84.10083 0.2  0.3811577      -       R
11  Gene11  68.71638 0.2  0.4969852      -       R
88  Gene89  64.38506 0.2  0.4804907      -       R
54  Gene54  59.93348 0.2  0.5642517      -       R
26  Gene26  58.38549 0.2  0.6215621      -       R
25  Gene25  51.70273 0.2  0.7066071      -       R
3    Gene3  46.00194 0.2  0.8028123      -       R
91  Gene92  45.90386 0.2  0.4838004      -       R
58  Gene58  45.50605 0.2  0.7348415      -       R
16  Gene16  43.84997 0.2  0.9786969      -       R
89  Gene90  42.72040 0.2  0.8680919      -       R
9    Gene9  42.01739 0.2  0.9558635      -       R
48  Gene48  41.82013 0.2  1.0728645      -       R
22  Gene22  40.45328 0.2  0.8163934      -       R
76  Gene76  38.54992 0.2  0.8814084      -       R
87  Gene88  35.14701 0.2  0.8821751      -       R
79  Gene80  34.41224 0.2  0.6539720      -       R
56  Gene56  34.25862 0.2  0.6983485      -       R
37  Gene37  64.09251 0.3  0.7072218      -       R
98  Gene99  48.72145 0.3  0.8812622      -       R
86  Gene87  44.27434 0.3  1.1906371      -       R
4    Gene4  38.79233 0.3  1.0669108      -       R
83  Gene84  31.73809 0.3  0.9771723      -       R
5    Gene5  73.18440 0.4  0.8496894      -       R
8    Gene8  35.75174 0.4  1.1839906      -       R
85  Gene86  31.39572 0.4  1.1725224      -       R
72  Gene72  54.06801 0.5  1.1406095      -       R
39  Gene39  44.21818 0.5  1.6425369      -       R
99 Gene100  30.78572 0.5  1.9949584      -       R
63  Gene63  65.35801 0.6  1.4317514      -       R
35  Gene35  34.61745 0.6  2.0142300      -       R
13  Gene13  33.95962 0.7  2.2632251      -       R
20  Gene20  29.47172 0.8  2.3888816      -       R
47  Gene47  38.60670 1.0  2.7468553      -       R
66  Gene66  24.24463 1.8  5.6432982      -       R
60  Gene60  23.04710 2.2  7.7529417      -       R
94  Gene95  19.42965 2.3  8.4947719      -       R
59  Gene59  18.74276 2.4  9.8450535      -       R
27  Gene27  24.79342 3.0  8.4300048      -       R
55  Gene55  26.19269 3.1  9.9946725      -       R
12  Gene12  20.62520 3.1 10.0703651      -       R
97  Gene98  18.32822 3.1 12.8366699      -       R
49  Gene49  23.11792 3.3 10.2350192      -       R
30  Gene30  18.37483 3.3  7.4443107      -       R
69  Gene69  22.76689 3.5 11.4504084      -       R
57  Gene57  18.91166 3.5 12.7353833      -       R
95  Gene96  16.85914 3.6 14.4007254      -       R

$clusters[[2]]
         id        OR   FP      dFP Strong Relaxed
46   Gene46 14.632080  4.3 15.79809      -       R
31   Gene31 17.597946  4.9 19.29094      -       R
45   Gene45 14.734285  4.9 20.34821      -       R
74   Gene74 17.892082  5.2 16.30316      -       R
84   Gene85 13.488934  5.2 21.94855      -       R
43   Gene43 14.577235  5.8 26.09057      -       R
33   Gene33 11.898146  5.9 22.86455      -       R
77   Gene78 10.339746  5.9 37.95774      -       R
75   Gene75 14.767943  6.5 30.38618      -       R
18   Gene18  9.955955  6.7 39.78764      -       R
19   Gene19 11.400054  6.8 32.18989      -       R
51   Gene51 11.009625  7.1 42.79650      -       R
93   Gene94  8.452701  8.0 58.82594      -       R
102 Gene399  7.992158  8.1 63.43035      -       R
100 Gene104  9.543673  8.3 55.43697      -       -
44   Gene44  9.085678  8.6 60.68339      -       -
38   Gene38  8.201157  8.7 75.19811      -       -
41   Gene41 14.155429  9.3 32.65473      -       -
103 Gene598  7.385797  9.9 45.75472      -       -
101 Gene277  8.133861 10.0 35.83662      -       -
104 Gene618  8.115441 10.0 83.18225      -       -
107 Gene946  8.108669 10.0 55.55365      -       -
105 Gene651  7.659571 10.0 50.64723      -       -
106 Gene670  7.613718 10.0 81.80362      -       -
```

As a rule of thumb, differentially expressed genes are expected to present high values of OR and low values of meanFP and density. We could also analyze each gene individually inside each cluster. The motivation of the clustering is to distinguish those false positives that score high in OR and low in meanFP and density, but are similar to other known false positives obtained by boostrapping. The procedure is detailed in the paper referenced above. 

If the researcher is interested in a more thorough analysis, other functions are at their service.

The data before being clustered can be obtained with the following function

```
> preclusteredData(myORdensity)
Columns "Strong" and "Relaxed" show the genes identified as DE genes
They denote the strong selection (FP=0) with S and the relaxed selection (FP < expectedFalsePositives) with F
         id         OR   FP        dFP Strong Relaxed
62   Gene62 175.043226  0.0  0.0000000      S       R
50   Gene50 172.287790  0.0  0.0000000      S       R
61   Gene61 155.536262  0.0  0.0000000      S       R
70   Gene70 152.547051  0.0  0.0000000      S       R
2     Gene2 149.653354  0.0  0.0000000      S       R
68   Gene68 148.872937  0.0  0.0000000      S       R
52   Gene52 134.877088  0.0  0.0000000      S       R
34   Gene34  77.322660  0.0  0.0000000      S       R
71   Gene71  74.376305  0.0  0.0000000      S       R
53   Gene53  69.740608  0.0  0.0000000      S       R
28   Gene28  62.752272  0.0  0.0000000      S       R
80   Gene81  54.302410  0.0  0.0000000      S       R
40   Gene40 104.994030  0.1  0.1991604      -       R
81   Gene82  93.452660  0.1  0.2892853      -       R
15   Gene15  93.163106  0.1  0.2956274      -       R
64   Gene64  92.594912  0.1  0.2513645      -       R
73   Gene73  91.677327  0.1  0.2812384      -       R
67   Gene67  88.445307  0.1  0.3110092      -       R
29   Gene29  88.439498  0.1  0.3097112      -       R
6     Gene6  84.472749  0.1  0.3047120      -       R
23   Gene23  76.090339  0.1  0.2390754      -       R
90   Gene91  52.915740  0.1  0.3281189      -       R
21   Gene21  50.880955  0.1  0.3199397      -       R
96   Gene97  48.442795  0.1  0.2908925      -       R
42   Gene42  44.777561  0.1  0.5999895      -       R
1     Gene1  40.134270  0.1  0.3922824      -       R
82   Gene83  39.639390  0.1  0.4893582      -       R
14   Gene14  38.549462  0.1  0.4340459      -       R
78   Gene79  35.045911  0.1  0.3481011      -       R
32   Gene32 144.807897  0.2  0.1803355      -       R
7     Gene7 143.422006  0.2  0.1889123      -       R
10   Gene10 130.618437  0.2  0.2217622      -       R
36   Gene36 120.829166  0.2  0.2507271      -       R
92   Gene93 103.662961  0.2  0.3341376      -       R
65   Gene65  98.361847  0.2  0.3128623      -       R
24   Gene24  92.201464  0.2  0.3046823      -       R
17   Gene17  84.100834  0.2  0.3811577      -       R
11   Gene11  68.716381  0.2  0.4969852      -       R
88   Gene89  64.385062  0.2  0.4804907      -       R
54   Gene54  59.933478  0.2  0.5642517      -       R
26   Gene26  58.385494  0.2  0.6215621      -       R
25   Gene25  51.702726  0.2  0.7066071      -       R
3     Gene3  46.001940  0.2  0.8028123      -       R
91   Gene92  45.903855  0.2  0.4838004      -       R
58   Gene58  45.506048  0.2  0.7348415      -       R
16   Gene16  43.849974  0.2  0.9786969      -       R
89   Gene90  42.720399  0.2  0.8680919      -       R
9     Gene9  42.017391  0.2  0.9558635      -       R
48   Gene48  41.820128  0.2  1.0728645      -       R
22   Gene22  40.453284  0.2  0.8163934      -       R
76   Gene76  38.549923  0.2  0.8814084      -       R
87   Gene88  35.147011  0.2  0.8821751      -       R
79   Gene80  34.412237  0.2  0.6539720      -       R
56   Gene56  34.258621  0.2  0.6983485      -       R
37   Gene37  64.092508  0.3  0.7072218      -       R
98   Gene99  48.721454  0.3  0.8812622      -       R
86   Gene87  44.274336  0.3  1.1906371      -       R
4     Gene4  38.792334  0.3  1.0669108      -       R
83   Gene84  31.738086  0.3  0.9771723      -       R
5     Gene5  73.184401  0.4  0.8496894      -       R
8     Gene8  35.751735  0.4  1.1839906      -       R
85   Gene86  31.395715  0.4  1.1725224      -       R
72   Gene72  54.068008  0.5  1.1406095      -       R
39   Gene39  44.218182  0.5  1.6425369      -       R
99  Gene100  30.785718  0.5  1.9949584      -       R
63   Gene63  65.358008  0.6  1.4317514      -       R
35   Gene35  34.617447  0.6  2.0142300      -       R
13   Gene13  33.959617  0.7  2.2632251      -       R
20   Gene20  29.471715  0.8  2.3888816      -       R
47   Gene47  38.606697  1.0  2.7468553      -       R
66   Gene66  24.244631  1.8  5.6432982      -       R
60   Gene60  23.047097  2.2  7.7529417      -       R
94   Gene95  19.429654  2.3  8.4947719      -       R
59   Gene59  18.742761  2.4  9.8450535      -       R
27   Gene27  24.793415  3.0  8.4300048      -       R
55   Gene55  26.192693  3.1  9.9946725      -       R
12   Gene12  20.625202  3.1 10.0703651      -       R
97   Gene98  18.328224  3.1 12.8366699      -       R
49   Gene49  23.117921  3.3 10.2350192      -       R
30   Gene30  18.374833  3.3  7.4443107      -       R
69   Gene69  22.766890  3.5 11.4504084      -       R
57   Gene57  18.911663  3.5 12.7353833      -       R
95   Gene96  16.859135  3.6 14.4007254      -       R
46   Gene46  14.632080  4.3 15.7980896      -       R
31   Gene31  17.597946  4.9 19.2909380      -       R
45   Gene45  14.734285  4.9 20.3482096      -       R
74   Gene74  17.892082  5.2 16.3031581      -       R
84   Gene85  13.488934  5.2 21.9485478      -       R
43   Gene43  14.577235  5.8 26.0905728      -       R
33   Gene33  11.898146  5.9 22.8645499      -       R
77   Gene78  10.339746  5.9 37.9577425      -       R
75   Gene75  14.767943  6.5 30.3861762      -       R
18   Gene18   9.955955  6.7 39.7876362      -       R
19   Gene19  11.400054  6.8 32.1898937      -       R
51   Gene51  11.009625  7.1 42.7965016      -       R
93   Gene94   8.452701  8.0 58.8259400      -       R
102 Gene399   7.992158  8.1 63.4303484      -       R
100 Gene104   9.543673  8.3 55.4369711      -       -
44   Gene44   9.085678  8.6 60.6833925      -       -
38   Gene38   8.201157  8.7 75.1981087      -       -
41   Gene41  14.155429  9.3 32.6547299      -       -
103 Gene598   7.385797  9.9 45.7547190      -       -
101 Gene277   8.133861 10.0 35.8366225      -       -
104 Gene618   8.115441 10.0 83.1822532      -       -
107 Gene946   8.108669 10.0 55.5536499      -       -
105 Gene651   7.659571 10.0 50.6472306      -       -
106 Gene670   7.613718 10.0 81.8036205      -       -
```

A plot with a representation of the potential genes based on OR (vertical axis), FP (horizontal axis) and dFP (size of the circle is inversely proportional to its value) can also be obtained. Genes that fulfil the relaxed criterion are drawn with triangles. The resulting plot is similar to Fig.3b in \[1\].

```
plot(myORdensity)
```

![plot1](/images/plot.png)

By default, the number of clusters computed by the ORdensity method is used. Other values for the number of clusters can be specified.

```
plot(myORdensity, k = 5)
```

![plot2](/images/plot5.png)



The plot of k values against the silhouette measure is also provided.

```
silhouetteAnalysis(myORdensity)
```

![plot2](/images/silhouetteAnalysis.png)



