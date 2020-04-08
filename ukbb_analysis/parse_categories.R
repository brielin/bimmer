categories <- readr::read_tsv("/gpfs/commons/projects/UKBB/sumstats/ukbb_categories.txt") %>%
  mutate(category_name = title) %>% select(category_id, category_name)
fields <- readr::read_tsv("/gpfs/commons/projects/UKBB/sumstats/ukbb_fields.txt") %>%
  mutate(trait_id = field_id, trait_name = title) %>% select(trait_id, trait_name, main_category)
fields <- dplyr::left_join(fields, categories, by = c("main_category" = "category_id"))

