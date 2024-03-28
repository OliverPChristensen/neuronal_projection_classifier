library('tidyverse')
library('Seurat')
library('qs')

neurons <- qread('./raw_data/7plex_neurons_filtered_labelled.qs')

head(neurons)

sum(neurons$infected == FALSE)
sum(neurons$infected == TRUE)

neurons$v147

#CR: cell ranger
#TS: Tap seq
#Lane: pos or neg

unique(neurons$sort)

table(neurons$sort)
neurons@meta.data %>% ggplot() + geom_density(aes(x = as.numeric(sort)))

neurons$lane[10000 + 1:100]
neurons$sort[10000 + 1:100]


ord <- order(unique(neurons@meta.data[,c("area","sort")])[,2])

unique(neurons@meta.data[,c("area","sort")])[ord,]


t <- as_tibble(table(neurons@meta.data[,c("area","sort")]))

t %>% filter(n > 0) %>% arrange(sort,n) -> t

print(t,n = 100)

#CEA High and low prob
#BST High and low prob
#Måske PVH, men ikke så mange gode negativer

neurons@meta.data %>% filter(area == 'CEA', sort == 1) %>% ggplot() + geom_histogram(aes(x = v147))
neurons@meta.data %>% filter(area == 'PVH', sort == 1) %>% ggplot() + geom_histogram(aes(x = v147))
