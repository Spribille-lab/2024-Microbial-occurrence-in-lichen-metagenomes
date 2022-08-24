# Reanalyzing metagenomes

## Yeast presence

### Dothist
1. Plotted all metagenomes and their yeast detection as dothist in `results/figures/metagenome_reanalysis_dothist.png`

Reads have higher detection rate, coverage seems to affect. Next, tried to slice the dataset exploring what are the trends.

2. Yeast detection juxtaposed with architecture (macro/crust) in Lecanorales as dothist in `results/figures/metagenome_reanalysis_dothist_lecanorales_arcitechture.png`

Yeasts are way more prevalent in macro than in crusts given similar sequencing depth.

3. Yeast detection juxtaposed with architecture (macro/crust) in ALL as dothist in `results/figures/metagenome_reanalysis_dothist_arcitechture.png`

Same pattern as before

4. Saved this as a table in `results/tables/summary_yeast_architecture.txt` and `results/tables/summary_yeast_architecture_lecanorales.txt`


```
> df_summary2
# A tibble: 4 x 5
# Groups:   architecture, yeast [4]
architecture yeast    absent `present in assembly and reads` `present in reads`
<fct>        <fct>     <dbl>                           <dbl>              <dbl>
1 crust        cypho      81.1                            3.88               15.0
2 crust        tremella   49.5                           16.5                34.0
3 macro        cypho      72.9                            8.60               18.6
4 macro        tremella   40.3                           19.9                39.8

> df_summary2_lecanorales
# A tibble: 4 x 5
# Groups:   architecture, yeast [4]
architecture yeast    absent `present in reads` `present in assembly and reads`
<fct>        <fct>     <dbl>              <dbl>                           <dbl>
1 crust        cypho      84.3               15.7                            NA
2 crust        tremella   51.0               29.4                            19.6
3 macro        cypho      64.5               23.9                            11.6
4 macro        tremella   26.1               46.4                            27.5
```

Can see that:
* In Lecanorales crusts prevalence is the same as in all crusts, in Lecanorales macro it's higher (not surprising, a bulk of non-Lecanorales macro are Peltigerales that don't have much)
* In all macro had higher prevalence of both yeasts
* In Tremella the difference is actually bigger (though this might be caused by the extremely low coverage metagenomes that exist in crusts but not in macro)
```
> (84.3-64.5)/84.3
[1] 0.2348754
> (51-26.1)/51
[1] 0.4882353
```
5. Plotted Peltigerales in `results/figures/metagenome_reanalysis_dothist_peltigerales.png`. Not much ging on, Yeasts mostly appear at very high depth (and could be contamination). A contrast with Lecanorales

6. Saved parmeliaceae graph in `results/figures/metagenome_reanalysis_dothist_parmeliaceae.png`. Cypho is present in all metagenomes >5Gb, but missing in a bunch of low coverages. Tremella is way more prevalent, but missing from one high coverage
```
> df[df$family=="Parmeliaceae",c(3,8)] %>% group_by(yeast,presence)%>%summarise(n=n(),freq_percent=n*100/length(unique(df[df$family=="Parmeliaceae",1])))
`summarise()` has grouped output by 'yeast'. You can override using the `.groups` argument.
# A tibble: 6 x 4
# Groups:   yeast [2]
yeast    presence                          n  freq_percent
<fct>    <chr>                         <int> <dbl>
1 cypho    absent                           61  63.5
2 cypho    present in assembly and reads    16  16.7
3 cypho    present in reads                 19  19.8
4 tremella absent                           24  25
5 tremella present in assembly and reads    32  33.3
6 tremella present in reads                 40  41.7
```

7. Saved lecanorales graph in `metagenome_reanalysis_dothist_lecanorales.png`. Tremella prevalence is high, but hard to know whether these are lichen-associated tremellas

8. What are non-lecanorales metagenomes with cyphobasidium?
```
nonlecanorales_cypho_table<-df %>% filter(order!="Lecanorales",yeast=="cypho",presence!="absent") %>% select(c(Lichen.metagenomes,bp,presence,architecture,family, order, class))
write.table(nonlecanorales_cypho_table,"results/tables/nonlecanorales_cypho_table.txt",row.names = F,quote = F)
```

## Metagenome by source
R script in `code/metagenome_by_source_fig.R`, resulting figure in `results/figures/metagenome_dothist_by_source.png`. Nicely shows the distribution of coverage in all metagenomes and the paper the metagenomes came from. Can see that Lendemer's data is a big portion of the whole data set, and also is low coverage

## Size of assembly
Needed to know how many metagenomes fail to produce a functional assembly. Added a rule to snakemake, the R script handling the result is in 
