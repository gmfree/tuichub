##stacked bar plots

# Using physeq.rare to make community comp figure 

# below is aggregating ASVs by phylum 
phyla_counts_tab <- otu_table(tax_glom(ps.rare, taxrank="Phylum")) 

# aggregating ASVs by class 
class_counts_tab <- otu_table(tax_glom(ps.rare, taxrank="Class"))
class_tax_vec <- as.vector(tax_table(tax_glom(ps.rare, taxrank="Class"))[,"Class"]) 
rownames(class_counts_tab) <- as.vector(class_tax_vec)

# making rownames phyla 
phyla_tax_vec <- as.vector(tax_table(tax_glom(ps.rare, taxrank="Phylum"))[,"Phylum"]) 
colnames(phyla_counts_tab) <- as.vector(phyla_tax_vec)

# now we have a dataframe with the number of reads corresponding for each phylum and class in each sample 

# making a relative abundance df at the phylum level 
phylum_proportions <- apply(phyla_counts_tab, 1, function(x) x/sum(x))
phylum_proportions <- as.data.frame(phylum_proportions) # making df
phylum_proportions$Phyla <- row.names(phylum_proportions) # column of phyla

# making data long - sample column, phyla column, relative abundance column 
phylum_prop_long <- phylum_proportions %>% 
  pivot_longer(cols = colnames(phylum_proportions[,1:83]), 
               names_to = "id", 
               values_to = "rel_abund")

# ordering by sample ID 
phylum_prop_long <- phylum_prop_long[order(phylum_prop_long$id),] 

# now we have long df with each sample and the rel abundance of phyla  
```

merging with metadata - specifically site level and elev.cat level 
```{r}
# manipulating metadata file: 
metadata_final = data.frame (sample_data(ps.rare))
metadata_final$id <- metadata_final$no
str(metadata_final)

comp.metadata <- as.data.frame(metadata_final[, c("id", "location", "loc_type", "mesocosm")])

# merging metadata with phyla df 
phyla_comp <-  inner_join(phylum_prop_long, comp.metadata, by = "id")
```

# plotting Phylum composition 
```{r}
# basic figure of mean relative abundance for above and below treeline
phyla_location <- phyla_comp %>% 
  group_by(location, Phyla) %>% 
  summarize(mean_rel_abund = 100*mean(rel_abund)) %>% 
  ggplot(aes(x=location, y=mean_rel_abund, fill=Phyla)) + 
  geom_col() + 
  theme_classic() + 
  scale_fill_discrete(name = NULL) + 
  labs(y = "Mean Relative Abundance (%)", x=NULL) + 
  theme(legend.text = element_text(face = "italic"), 
        legend.key.size = unit(10, "pt"))
phyla_location
ggsave("figures/mean_phyla_location_bar_raw.tiff", width = 12, height = 5)










phyla_counts_tab_fin <- otu_table(tax_glom(ps.fin, taxrank="Phylum")) 
phyla_tax_vec_fin <- as.vector(tax_table(tax_glom(ps.fin, taxrank="Phylum"))[,"Phylum"]) 
colnames(phyla_counts_tab_fin) <- as.vector(phyla_tax_vec_fin)

class_counts_tab_fin <- otu_table(tax_glom(ps.fin, taxrank="Class"))
class_tax_vec_fin <- as.vector(tax_table(tax_glom(ps.fin, taxrank="Class"))[,"Class"]) 
colnames(class_counts_tab_fin) <- as.vector(class_tax_vec_fin)

phylum_proportions <- apply(phyla_counts_tab, 2, function(x) x/sum(x))
phylum_prop_transpose <- t(phylum_proportions)
phylum_prop_transpose <- as.data.frame(phylum_prop_transpose)

phylum_prop_transpose$Phyla <- row.names(phylum_prop_transpose)


# making data long - sample column, phyla column, relative abundance column 
phylum_prop_long <- phylum_prop_transpose %>% 
  pivot_longer(cols = colnames(phylum_prop_transpose[,1:83]), 
               names_to = "id", 
               values_to = "rel_abund")

# ordering by sample ID 
phylum_prop_long <- phylum_prop_long[order(phylum_prop_long$id),] phylum_prop_long
# now we have long df with each sample and the rel abundance of phyla 






####phyla_location_2percent (top 10 phyla)
# pooling everything over 2% using a logical 
phyla_pool_2 <- phyla_mean_rel_abund %>% 
  group_by(Phyla) %>% 
  summarize(pool = max(mean_rel_abund) < 2, 
            mean = mean(mean_rel_abund)) # this creates a logical of if a phyla meets the condition (is less than 2%)

# join with df and plot! 
inner_join(phyla_mean_rel_abund, phyla_pool_2, by = "Phyla") %>% 
  mutate(Phyla = if_else(pool, "Other", Phyla)) %>%
  group_by(location, Phyla) %>% 
  summarize(mean_rel_abund = sum(mean_rel_abund),
            mean = sum(mean)) %>% 
  mutate(Phyla = factor(Phyla), 
         Phyla = fct_reorder(Phyla, mean, .desc = TRUE), 
         Phyla = fct_shift(Phyla, n=1)) %>%
  ggplot(aes(x=location, y=mean_rel_abund, fill=Phyla)) + 
  geom_col() + 
  theme_classic() + 
  scale_fill_manual(name = NULL, 
                    breaks = c("Other", "Firmicutes", "Proteobacteria", "Fusobacteriota", "Cyanobacteria", "Actinobacteriota", "Spirochaetota", "Bacteroidota", "Planctomycetota", "Chloroflexi", "Verrucomicrobiota"), 
                    values = c(brewer.pal(12, "Set3"))) +
  labs(y = "Mean Relative Abundance (%)", x=NULL) + 
  theme(legend.text = element_text(face = "italic"), 
        legend.key.size = unit(12, "pt"))
ggsave("figures/phyla_location_2percent.tiff", height = 5, width = 12)

