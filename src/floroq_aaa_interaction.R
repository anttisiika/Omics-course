# 
library(tidyverse)
library(broom)
library(stringr)

# Load drug interactions file. 
read_csv('drug_interactions.csv') %>%  
    filter(match == 'TRUE')  %>% 
    select(-X1, -match) -> 
    drug_interactions

read_csv('../data/beta.csv') %>%  
    filter(match == 'TRUE')  %>% 
    select(-X1, -match) -> 
    beta_lactams


### ggplot2 code-snippets 

remove_x_axes <- theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

remove_y_axes <- theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())



## load gene expr data

gene_expr <- read_delim('../data/aaa_array.csv', delim = ';')
colnames(gene_expr)[1] <- 'ProbeSetId'
colnames(gene_expr)[2:43] <- str_split_fixed(colnames(gene_expr)[2:43], '-', 2)[,1] %>% 
        substr(5,16) # process sample-names into ids + condition. 

## Load features and delete all ProbeSetIds that are not associated with
## any geneSymbol (i.e '---').

features <- read_delim('../data/features.csv', delim = ';', na = c('', '---'))  %>%  
    select(`Probe Set ID`, `Gene Symbol`)  %>% 
    na.omit() # omit all mRNAs that are not catalogued
colnames(features) <- c('ProbeSetId', 'geneName')
features$geneName <- str_split_fixed(features$geneName, ' ', 2)[,1] # select only the first gene name if gene is 'XXXX // YYYY // ZZZZ'
         
## Now we join the gene-expr to the features to match geneNames, here we want
# a df with equal n_rows as set of features and all columns from the gene
# expr-set (+1 for the geneNames). 

df <- left_join(features, gene_expr, by = 'ProbeSetId')  %>% 
    select(-ProbeSetId) 

## 'df' contains more than one entry per geneName, so here we simply compute the mean-value for each geneName. 
df %>% 
    group_by(geneName) %>%
    summarize_all(funs(mean)) -> df

## put into long format for easier data manipulation.
df  %>% 
    reshape2::melt(id.vars = 'geneName') ->
    df_long

#extract disease status and adv/med condition and id. 
df_long$status <- str_sub(df_long$variable, 1, 1) 
df_long$id <- str_sub(df_long$variable, 1,5)
df_long$tissue <- str_sub(df_long$variable, 10,13)

### Ranked gene list by t.test 
# note this seems pretty slow. 
t.ranked_diff_genes_med <- df_long %>%
    filter(tissue == 'med')  %>% 
    group_by(geneName) %>%
    do(tidy(t.test(.$value ~ .$status)))  %>% 
    mutate(fdr = p.adjust(p.value, method = "fdr"))  %>% 
    arrange(desc(estimate))

t.ranked_diff_genes_adv <- df_long %>%
    filter(tissue == 'adv')  %>% 
    group_by(geneName) %>%
    do(tidy(t.test(.$value ~ .$status)))  %>% 
    mutate(fdr = p.adjust(p.value, method = "fdr"))  %>% 
    arrange(desc(estimate))