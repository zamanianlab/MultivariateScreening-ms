library(tidyverse)
library(conflicted)
library(ZamanianLabThemes)
library(tidymodels)

# conflict resolution
conflict_prefer("filter", "dplyr")


# tidymodels --------------------------------------------------------------

model_data <- read_rds(here("Fig5/data/fecundity_model_data.rds")) %>%
  filter(!is.na(total), !is.na(pixel_count))

# 10-fold cross-validation
set.seed(123)
split <- initial_split(model_data,
  strata = pixel_count
)

train <- training(split)
test <- testing(split)

set.seed(345)
folds <- vfold_cv(train, v = 10)

# make recipe
lm_rec <-
  recipe(total ~ ., data = train) %>%
  update_role(plate:assay_date, new_role = "ID") %>%
  update_role(mf, pretzel, embryo, new_role = "NA") %>%
  update_role(pixel_count, new_role = "predictor")

# set the model, engine, and fit
lm_mod <-
  linear_reg() %>%
  set_engine("lm")

set.seed(234)
lm_fit <-
  lm_mod %>%
  fit(total ~ pixel_count, data = train)

# generate the worklflow
lm_wf <-
  workflow() %>%
  add_model(lm_mod) %>%
  add_recipe(lm_rec)

# perform the cross-validation
set.seed(456)
lm_fit_rs <-
  lm_wf %>%
  fit_resamples(folds)

collect_metrics(lm_fit_rs)

# evaluate on the test data
lm_testing_pred <-
  augment(lm_fit, test) %>%
  mutate(split = "Test") %>%
  bind_rows(train %>% mutate(split = "Train"))

lm_testing_pred %>%
  rsq(truth = total, .pred)

(model_performance <- lm_testing_pred %>%
  ggplot() +
  geom_line(data = . %>% drop_na(.pred), aes(x = pixel_count, y = .pred)) +
  geom_point(aes(x = pixel_count, y = total, color = split), alpha = 0.75) +
  ggtext::geom_richtext(
    data = . %>% rsq(truth = total, .pred),
    aes(label = paste0("R<sub>2</sub> = ", round(.estimate, digits = 2))),
    x = 250000, y = 1250, color = "white", fill = "black", alpha = 0.75
  ) +
  scale_color_manual(values = wesanderson::wes_palette("Darjeeling1")) +
  scale_fill_manual(values = wesanderson::wes_palette("Darjeeling1")) +
  lims(x = c(0, 1500000)) +
  labs(x = "Pixel Count", y = "Total Progeny", color = "") +
  theme_nw2() +
  theme(legend.position = "right",
        axis.text.x = element_markdown(angle = 0, hjust = 0.5)) +
  NULL
)

cowplot::save_plot(here("Fig5/supplementary/SupplementaryFig5.pdf"), model_performance, base_height = 4, base_asp = 1)
