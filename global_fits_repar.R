
library(tidyverse)
library(readxl)
library(bioinactivation2)
library(FME)
library(ggsci)
library(cowplot)

## Load data

d <- read_excel("./data/Listeria data_latest.xlsx", sheet = "Iso exp") %>%
  mutate(Strain = factor(Strain)) %>%
  rename(time = Time, logN = Log, temp = Temp) %>%
  left_join(
    .,
    tribble(
      ~Strain, ~strain_name,
      "1", "EGD-e (WT)",
      "2", "ΔsigB",
      "3", "RsbR1 only",
      "4", "rsbS S56A",
      "5", "rsbT N49A",
      "6", "rsbR1 T241A"
    )
  )

ggplot(d, aes(x = time, y = logN, colour = strain_name)) +
  geom_point(aes(shape = out)) +
  geom_line() +
  facet_wrap(Strain ~ temp, scales = "free")

d <- d %>%
  filter(is.na(out)) 

ggplot(d, aes(x = time, y = logN, colour = strain_name)) +
  geom_point() +
  # geom_line() +
  facet_wrap(temp ~ strain_name, scales = "free")

## Functions for One step fitting


#' p = c(logDref, z, logC0, logN0)
#'
make_pred <- function(p, temp, times, Tref) {
  
  p <- as.list(p)
  
  Dref <- 10^p$logDref
  kref <- log(10)/Dref
  p$logkref <- log(kref)
  p$b <- log(10)/p$z
  p$a <- p$logkref - p$b*Tref
  
  ## Secondary models
  
  SL <- pred_SL(p, temp)
  logk <- pred_k(p, temp)
  k <- exp(logk)
  
  ## Primary model
  
  logN0 <- p$logN0
  # logN0 <- 8
  
  N <- 10^logN0 * exp(-k*times) * exp(k*SL) / ( 1 + ( exp(k*SL) - 1 )*exp(-k*times) )
  
  # list(SL = SL, k = k, N = N)
  
  tibble(time = times, temp = temp, SL = SL, k = k, logN = log10(N))
  
}

pred_SL <- function(p, temp, Tref) {
  
  p <- as.list(p)
  
  p$a
  
  log(10^p$logC0 + 1)/exp(p$a + p$b*temp)
}

pred_k <- function(p, temp, Tref) {
  p <- as.list(p)
  p$a + p$b*temp
}

cost_onestep <- function(p, this_data, Tref) {
  
  p <- as.list(p)
  
  pred <- make_pred(p, this_data$temp, this_data$time, Tref)
  
  pred$logN - this_data$logN
  
}

## Fitting the reparameterized model

Tref <- 58

global_models <- d %>%
  split(.$strain_name) %>%
  map(.,
      ~ modFit(cost_onestep,
               # p = c(coef(sec_models[[1]]), logN0 = 6),
               p = c(logC0 = 0, logDref = 1, z = 8, logN0 = 6),
               this_data = .,
               lower = c(logC0 = -3, logDref = -5, z = 0, logN0 = 0),
               Tref = Tref
      )
  )

global_models %>% map(., ~ summary(.)$par)

## Table 2

global_models %>% 
  map(., ~ summary(.)$par) %>%
  map(., ~ as_tibble(., rownames = "par")) %>%
  map(., ~ mutate(., 
                  Estimate = round(Estimate, 2),
                  `Std. Error` = round(`Std. Error`, 2)
  )
  ) %>%
  map(.,
      ~ unite(., value, Estimate, `Std. Error`)
  ) %>%
  imap_dfr(., ~ mutate(.x, strain = .y)) %>%
  select(par, strain, value) %>%
  pivot_wider(names_from = par, values_from = value) %>%
  write_excel_csv2(., file = "Table_2.csv")

global_models %>% 
  map(., ~ tibble(r = residuals(.))) %>%
  map(., ~ mutate(., r2 = r^2)) %>%
  map(., ~ summarize(., ME = mean(r), MSE = mean(r2))) %>%
  map(., ~ mutate(., RMSE = sqrt(MSE)))

## Supp. Figure 1

p <- global_models %>%
  map(coef) %>%
  map(.,
      ~ make_pred(., d$temp, d$time, Tref = Tref)
  ) %>%
  imap_dfr(., ~ mutate(.x, strain_name = .y)) %>%
  # mutate(Strain = paste0("Strain: ", Strain)) %>%
  ggplot(aes(x = time, y = logN, colour = factor(temp))) +
  geom_line(linewidth = 1, linetype = 2) +
  geom_point(data = d,
             shape = 1, size = 3) +
  facet_wrap("strain_name") +#, scales = "free") +
  coord_cartesian(ylim = c(0, 5.5)) +
  labs(x = "Treatment time (min)",
       y = "Microbial concentration (log CFU/mL)") +
  # theme_bw(base_size = 14) +
  ggthemes::theme_few(base_size = 14) +
  theme(legend.position = "none") +
  scale_color_npg()

ggsave(p, filename = "supp_Figure_1.png", width = 12, height = 6)

# ## Visualizations
# 
# global_models %>%
#   map(coef) %>%
#   map(.,
#       ~ make_pred(., d$temp, d$time, Tref = Tref)
#   ) %>%
#   imap_dfr(., ~ mutate(.x, Strain = .y)) %>%
#   ggplot(aes(x = time, y = logN, colour = factor(temp))) +
#   geom_line() +
#   geom_point(data = d) +
#   facet_wrap("Strain", scales = "free")
# 
# global_models %>%
#   map(coef) %>%
#   map(.,
#       ~ make_pred(., d$temp, d$time, Tref = Tref)
#   ) %>%
#   imap_dfr(., ~ mutate(.x, Strain = .y)) %>%
#   ggplot(aes(x = time, y = logN, colour = Strain)) +
#   geom_line() +
#   geom_point(data = d) +
#   facet_wrap("temp", scales = "free")
# 
# global_models %>%
#   map(coef) %>%
#   map(.,
#       ~ make_pred(., d$temp, d$time, Tref = Tref)
#   ) %>%
#   imap_dfr(., ~ mutate(.x, Strain = .y)) %>%
#   ggplot(aes(x = time, y = logN, colour = Strain)) +
#   geom_line() +
#   geom_point(data = d) +
#   facet_wrap(Strain ~ temp, scales = "free")

## Load stationary data

d_stat <- read_excel("./data/Listeria data_latest.xlsx", sheet = "Iso sta") %>%
  mutate(Strain = factor(Strain)) %>%
  rename(time = `Time (min)`, logN = `log CFU/ml`, temp = `Temperature (ºC)`) %>%
  left_join(
    .,
    tribble(
      ~Strain, ~strain_name,
      "1", "EGD-e (WT)",
      "2", "ΔsigB",
      "3", "RsbR1 only",
      "4", "rsbS S56A",
      "5", "rsbT N49A",
      "6", "rsbR1 T241A"
    )
  )

ggplot(d_stat, aes(x = time, y = logN, colour = strain_name)) +
  geom_line() +
  geom_point(aes(shape = out)) +
  facet_wrap(Strain ~ temp, scales = "free")

d_stat <- d_stat %>%
  filter(is.na(out))

ggplot(d_stat, aes(x = time, y = logN, colour = strain_name)) +
  geom_point() +
  facet_wrap(temp ~ Strain, scales = "free")

## One step fitting

global_models_stat <- d_stat %>%
  filter(Strain %in% c("1", "2", "3")) %>%
  split(.$strain_name, drop = TRUE) %>%
  map(.,
      ~ modFit(cost_onestep,
               # p = c(coef(sec_models[[1]]), logN0 = 6),
               p = c(logC0 = .8, logDref = 0, z = 8, logN0 = 6),
               this_data = .,
               lower = c(logC0 = -1, a = -50, z = 0, logN0 = 0),
               Tref = Tref
      )
  )

global_models_stat %>% map(., ~ summary(.)$par)

## Table 3

global_models_stat %>% 
  map(., ~ summary(.)$par) %>%
  map(., ~ as_tibble(., rownames = "par")) %>%
  map(., ~ mutate(., 
                  Estimate = round(Estimate, 2),
                  `Std. Error` = round(`Std. Error`, 2)
  )
  ) %>%
  map(.,
      ~ unite(., value, Estimate, `Std. Error`)
  ) %>%
  imap_dfr(., ~ mutate(.x, strain = .y)) %>%
  select(par, strain, value) %>%
  pivot_wider(names_from = par, values_from = value) %>%
  write_excel_csv2(., file = "Table_3.csv")

global_models_stat %>% 
  map(., ~ tibble(r = residuals(.))) %>%
  map(., ~ mutate(., r2 = r^2)) %>%
  map(., ~ summarize(., ME = mean(r), MSE = mean(r2))) %>%
  map(., ~ mutate(., RMSE = sqrt(MSE)))

## Supp. Figure 2

p <- global_models_stat %>%
  map(coef) %>%
  map(.,
      ~ make_pred(., d_stat$temp, d_stat$time, Tref = Tref)
  ) %>%
  imap_dfr(., ~ mutate(.x, strain_name = .y)) %>%
  ggplot(aes(x = time, y = logN, colour = factor(temp))) +
  geom_line(linewidth = 1, linetype = 2) +
  geom_point(data = d_stat %>% 
               filter(Strain %in% c("1", "2", "3")),
             shape = 1, size = 4) +
  facet_wrap("strain_name") + #, scales = "free") +
  coord_cartesian(ylim = c(0, 6.5)) +
  labs(x = "Treatment time (min)",
       y = "Microbial concentration (log CFU/mL)") +
  ggthemes::theme_few(base_size = 14) +
  theme(legend.position = "none") +
  scale_color_npg()
p

ggsave(p, filename = "supp_Figure_2.png", width = 15, height = 5)

# ## Visualizations
# 
# global_models_stat %>%
#   map(coef) %>%
#   map(.,
#       ~ make_pred(., d_stat$temp, d_stat$time, Tref = Tref)
#   ) %>%
#   imap_dfr(., ~ mutate(.x, Strain = .y)) %>%
#   ggplot(aes(x = time, y = logN, colour = Strain)) +
#   geom_line() +
#   geom_point(data = d_stat) +
#   facet_wrap(Strain ~ temp, scales = "free")
# 
# 
# global_models_stat %>%
#   map(coef) %>%
#   map(.,
#       ~ make_pred(., d_stat$temp, d_stat$time, Tref = Tref)
#   ) %>%
#   imap_dfr(., ~ mutate(.x, Strain = .y)) %>%
#   ggplot(aes(x = time, y = logN, colour = factor(temp))) +
#   geom_line() +
#   geom_point(data = d_stat) +
#   facet_wrap("Strain", scales = "free")
# 
# global_models_stat %>%
#   map(coef) %>%
#   map(.,
#       ~ make_pred(., d_stat$temp, d_stat$time, Tref = Tref)
#   ) %>%
#   imap_dfr(., ~ mutate(.x, Strain = .y)) %>%
#   ggplot(aes(x = time, y = logN, colour = Strain)) +
#   geom_line() +
#   geom_point(data = d_stat) +
#   facet_wrap("temp", scales = "free")

## Figure 2

p1 <- global_models_stat %>%
  map(coef) %>%
  map(., ~ tibble(par = names(.), value = .)) %>%
  map(., ~ pivot_wider(., names_from = par, values_from = value)) %>%
  map(., 
      ~ mutate(., 
               b = log(10)/z,
               Dref = 10^logDref,
               kref = log(10)/Dref,
               logkref = log(kref),
               a = logkref - b*Tref
               )
      ) %>%
  map(as.list) %>%
  map(.,
      ~ tibble(temp = seq(55, 60, length = 100),
               SL = pred_SL(., temp)
      )
  ) %>%
  imap_dfr(., ~ mutate(.x, strain = .y)) %>%
  ggplot() +
  geom_line(aes(x = temp, y = SL, colour = strain),
            linetype = 1, linewidth = 1) +
  labs(x = "Treatment temperature (ºC)",
       y = "Shoulder length (min)") +
  ggthemes::theme_few(base_size = 14) +
  # theme_classic(base_size = 14) +
  theme(legend.position = "none") +
  scale_color_manual(values = pal_npg("nrc")(6)[c(1,2,6)])

p2 <- global_models %>%
  map(coef) %>%
  map(., ~ tibble(par = names(.), value = .)) %>%
  map(., ~ pivot_wider(., names_from = par, values_from = value)) %>%
  map(., 
      ~ mutate(., 
               b = log(10)/z,
               Dref = 10^logDref,
               kref = log(10)/Dref,
               logkref = log(kref),
               a = logkref - b*Tref
      )
  ) %>%
  map(as.list) %>%
  map(.,
      ~ tibble(temp = seq(55, 60, length = 100),
               SL = pred_SL(., temp)
      )
  ) %>%
  imap_dfr(., ~ mutate(.x, strain = .y)) %>%
  ggplot() +
  geom_line(aes(x = temp, y = SL, colour = strain),
            linetype = 1, linewidth = 1) +
  labs(x = "Treatment temperature (ºC)",
       y = "Shoulder length (min)") +
  ggthemes::theme_few(base_size = 14) +
  # theme(legend.position = "none",
  #       legend.title = element_blank()) +
  scale_color_npg() +
  theme(legend.title = element_blank())

p <- plot_grid(p1, p2, labels = "AUTO",
               rel_widths = c(.4, .6))
p
ggsave(p, filename = "Figure_2.png", width = 8, height = 4)

## Figure 1

p1 <- global_models_stat %>%
  map(coef) %>%
  map(., ~ tibble(par = names(.), value = .)) %>%
  map(., ~ pivot_wider(., names_from = par, values_from = value)) %>%
  map(., 
      ~ mutate(., 
               b = log(10)/z,
               Dref = 10^logDref,
               kref = log(10)/Dref,
               logkref = log(kref),
               a = logkref - b*Tref
      )
  ) %>%
  map(as.list) %>%
  map(.,
      ~ tibble(temp = seq(55, 60, length = 100),
               logk = pred_k(., temp),
               k = exp(logk),
               D = log(10)/k,
               logD = log10(D)
      )
  ) %>%
  imap_dfr(., ~ mutate(.x, strain = .y)) %>%
  ggplot() +
  geom_line(aes(x = temp, y = D, colour = strain),
            linetype = 1, linewidth = 1) +
  scale_y_log10() +
  labs(x = "Treatment temperature (ºC)",
       y = "D-value (min)") +
  ggthemes::theme_few(base_size = 14) +
  theme(legend.position = "none",
        legend.title = element_blank()) +
  scale_color_manual(values = pal_npg("nrc")(6)[c(1,2,6)])

p2 <- global_models %>%
  map(coef) %>%
  map(., ~ tibble(par = names(.), value = .)) %>%
  map(., ~ pivot_wider(., names_from = par, values_from = value)) %>%
  map(., 
      ~ mutate(., 
               b = log(10)/z,
               Dref = 10^logDref,
               kref = log(10)/Dref,
               logkref = log(kref),
               a = logkref - b*Tref
      )
  ) %>%
  map(as.list) %>%
  map(.,
      ~ tibble(temp = seq(55, 60, length = 100),
               logk = pred_k(., temp),
               k = exp(logk),
               D = log(10)/k,
               logD = log10(D)
      )
  ) %>%
  imap_dfr(., ~ mutate(.x, strain = .y)) %>%
  ggplot() +
  geom_line(aes(x = temp, y = D, colour = strain),
            linetype = 1, linewidth = 1) +
  scale_y_log10() +
  labs(x = "Treatment temperature (ºC)",
       y = "D-value (min)") +
  ggthemes::theme_few(base_size = 14) +
  # theme(legend.position = "none",
  #       legend.title = element_blank()) +
  scale_color_npg() +
  theme(legend.title = element_blank())

p <- plot_grid(p1, p2, labels = "AUTO",
               rel_widths = c(.4, .6))
p
ggsave(p, filename = "Figure_1.png", width = 8, height = 4)

################################################################################

## Load the dynamic data

d_dyna <- read_excel("./data/Listeria data_latest.xlsx", sheet = "Dynamic exp") %>%
  mutate(Strain = factor(Strain),
         Ramp = factor(Ramp)) %>%
  rename(time = Time, logN = Log) %>%
  left_join(
    .,
    tribble(
      ~Strain, ~strain_name,
      "1", "EGD-e (WT)",
      "2", "ΔsigB",
      "3", "RsbR1 only",
      "4", "rsbS S56A",
      "5", "rsbT N49A",
      "6", "rsbR1 T241A"
    )
  )

d_dyna %>%
  ggplot(aes(x = time, y = logN, colour = strain_name)) +
  geom_point() +
  facet_wrap(Ramp ~ Strain, scales = "free")

## Auxiliary stuff

my_COs <- global_models %>% map(., ~ 10^coef(.)[["logC0"]])

sec_models_dyna <- list(
  list(par = "D",
       model = "Bigelow",
       depends_on = c("temp")
  ),
  list(par = "Nres",
       model = "Bigelow",
       depends_on = c())
)

my_guess <- global_models %>%
  map(coef) %>%
  map(.,
      ~ list(guess = c(logN0 = 5,
                       D_logref = 0.5,
                       logC0 = 2
      ),
      known = c(D_temp_xref = Tref,
                Nres_logref = -1,
                D_temp_z = .[["z"]]
      )
      )
  )

## Fitting of ramp 2ºC/min

this_d <- d_dyna %>%
  filter(Ramp == 2)

logN0 <- this_d %>% filter(time == 0) %>% summarize(l = mean(logN)) %>% pull(l)

times <- seq(0, max(this_d$time), length = 1000)

env_conditions <- tibble(time = c(0, (58-30)/2, max(d$time)),
                         temp = c(30, 58, 58))

pred_geeraerd <- global_models %>%
  map(coef) %>%
  map(.,
      ~ list(
        list(par = "logD",
             model = "Bigelow",
             temp = c("xref" = Tref, z = .[["z"]]),
             ref = .[["logDref"]]
        )
      )
  ) %>%
  imap(.,
       ~ predict_inactivation(times,
                              list(model = "Geeraerd_noTail", logN0 = logN0,
                                   C0 = my_COs[[.y]]
                                   # C0 = my_C0
                              ),
                              environment = "dynamic",
                              .x,
                              env_conditions)
  )

pred_geeraerd %>%
  map(., ~.$simulation) %>%
  imap_dfr(., ~ mutate(.x, strain_name = .y)) %>%
  ggplot() +
  geom_line(aes(x = time, y = logN, colour = strain_name)) +
  geom_point(aes(x = time, y = logN, colour = strain_name),
             data = this_d) +
  facet_wrap("strain_name") +
  coord_cartesian(ylim = c(0, 5.5))

dyna_fit <- this_d %>%
  split(.$strain_name) %>%
  map(., ~ select(., time, logN)) %>%
  map2(., my_guess,
       ~ fit_inactivation("dynamic", 
                          .x, 
                          "Geeraerd",
                          .y$guess, 
                          .y$known, 
                          secondary_models = sec_models_dyna,
                          env_conditions = env_conditions,
                          algorithm = "regression"
       )
  )

dyna_fit %>% map(plot) %>% plot_grid(plotlist = .)
dyna_fit %>% map(summary)
dyna_fit %>% map(coef)

#- Parameter comparison

global_models %>% 
  map(coef) %>%
  map(., ~ tibble(par = names(.), value = .)) %>%
  map(., ~ pivot_wider(., names_from = par, values_from = value)) %>%
  # imap_dfr(., ~ mutate(.x, strain = .y)) %>%
  # map(.,
  #     ~   mutate(.,
  #                b = log(10)/z,
  #                logD0 = log10(log(10)/exp(a)),
  #                logD58 = logD0 - (58)/z,
  #                D58 = 10^logD58
  #     )
  # ) %>%
  map2(., map(dyna_fit, coef),
       ~ mutate(.x, logD_dyna = .y[["D_logref"]], logC0_dyna =.y[["logC0"]])
  ) %>%
  imap_dfr(., ~ mutate(.x, strain = .y)) %>%
  mutate(
    C0 = 10^logC0,
    C0_dyna = 10^logC0_dyna,
    D58_dyna = 10^logD_dyna,
    D58_iso = 10^logDref
  ) %>%
  mutate(rel_D58 = D58_dyna/D58_iso,
         relC0 = logC0_dyna - logC0) 

## Plot Ramp 5

this_d <- d_dyna %>%
  filter(Ramp == 5)

logN0 <- this_d %>% filter(time == 0) %>% summarize(l = mean(logN)) %>% pull(l)

times <- seq(0, max(this_d$time), length = 1000)

env_conditions <- tibble(time = c(0, (58-30)/5, max(d$time)),
                         temp = c(30, 58, 58))

pred_geeraerd_r5 <- global_models %>%
  map(coef) %>%
  map(.,
      ~ list(
        list(par = "logD",
             model = "Bigelow",
             temp = c("xref" = Tref, z = .[["z"]]),
             ref = .[["logDref"]]
        )
      )
  ) %>%
  imap(.,
       ~ predict_inactivation(times,
                              list(model = "Geeraerd_noTail", logN0 = logN0,
                                   C0 = my_COs[[.y]]
                                   # C0 = my_C0
                              ),
                              environment = "dynamic",
                              .x,
                              env_conditions)
  )

pred_geeraerd_r5 %>%
  map(., ~.$simulation) %>%
  imap_dfr(., ~ mutate(.x, strain_name = .y)) %>%
  ggplot() +
  geom_line(aes(x = time, y = logN, colour = strain_name)) +
  geom_point(aes(x = time, y = logN, colour = strain_name),
             data = this_d) +
  facet_wrap("strain_name") +
  coord_cartesian(ylim = c(0, 5.5))

dyna_fit_r5 <- this_d %>%
  split(.$strain_name) %>%
  map(., ~ select(., time, logN)) %>%
  map2(., my_guess,
       ~ fit_inactivation("dynamic", 
                          .x, 
                          "Geeraerd",
                          .y$guess, 
                          .y$known, 
                          secondary_models = sec_models_dyna,
                          env_conditions = env_conditions,
                          algorithm = "regression"
       )
  )

dyna_fit_r5 %>% map(plot) %>% plot_grid(plotlist = .)
dyna_fit_r5 %>% map(summary)
dyna_fit_r5 %>% map(coef)

#- Parameter comparison

global_models %>% 
  map(coef) %>%
  map(., ~ tibble(par = names(.), value = .)) %>%
  map(., ~ pivot_wider(., names_from = par, values_from = value)) %>%
  # imap_dfr(., ~ mutate(.x, strain = .y)) %>%
  # map(.,
  #     ~   mutate(.,
  #                b = log(10)/z,
  #                logD0 = log10(log(10)/exp(a)),
  #                logD58 = logD0 - (58)/z,
  #                D58 = 10^logD58
  #     )
  # ) %>%
  map2(., map(dyna_fit_r5, coef),
       ~ mutate(.x, logD_dyna = .y[["D_logref"]], logC0_dyna =.y[["logC0"]])
  ) %>%
  imap_dfr(., ~ mutate(.x, strain = .y)) %>%
  mutate(
    C0 = 10^logC0,
    C0_dyna = 10^logC0_dyna,
    D58_dyna = 10^logD_dyna,
    D58_iso = 10^logDref
  ) %>%
  mutate(rel_D58 = D58_dyna/D58_iso,
         relC0 = logC0_dyna - logC0) 

## Plot Ramp 10

this_d <- d_dyna %>%
  filter(Ramp == 10)

logN0 <- this_d %>% filter(time == 0) %>% summarize(l = mean(logN)) %>% pull(l)

times <- seq(0, max(this_d$time), length = 1000)

env_conditions <- tibble(time = c(0, (58-30)/10, max(d$time)),
                         temp = c(30, 58, 58))

pred_geeraerd_r10 <- global_models %>%
  map(coef) %>%
  map(.,
      ~ list(
        list(par = "logD",
             model = "Bigelow",
             temp = c("xref" = Tref, z = .[["z"]]),
             ref = .[["logDref"]]
        )
      )
  ) %>%
  imap(.,
       ~ predict_inactivation(times,
                              list(model = "Geeraerd_noTail", logN0 = logN0,
                                   C0 = my_COs[[.y]]
                                   # C0 = my_C0
                              ),
                              environment = "dynamic",
                              .x,
                              env_conditions)
  )

pred_geeraerd_r10 %>%
  map(., ~.$simulation) %>%
  imap_dfr(., ~ mutate(.x, strain_name = .y)) %>%
  ggplot() +
  geom_line(aes(x = time, y = logN, colour = strain_name)) +
  geom_point(aes(x = time, y = logN, colour = strain_name),
             data = this_d) +
  facet_wrap("strain_name") +
  coord_cartesian(ylim = c(0, 5.5))

dyna_fit_r10 <- this_d %>%
  split(.$strain_name) %>%
  map(., ~ select(., time, logN)) %>%
  map2(., my_guess,
       ~ fit_inactivation("dynamic", 
                          .x, 
                          "Geeraerd",
                          .y$guess, 
                          .y$known, 
                          secondary_models = sec_models_dyna,
                          env_conditions = env_conditions,
                          algorithm = "regression"
       )
  )

dyna_fit_r10 %>% map(plot) %>% plot_grid(plotlist = .)
dyna_fit_r10 %>% map(summary)
dyna_fit_r10 %>% map(coef)

#- Parameter comparison

global_models %>% 
  map(coef) %>%
  map(., ~ tibble(par = names(.), value = .)) %>%
  map(., ~ pivot_wider(., names_from = par, values_from = value)) %>%
  # imap_dfr(., ~ mutate(.x, strain = .y)) %>%
  # map(.,
  #     ~   mutate(.,
  #                b = log(10)/z,
  #                logD0 = log10(log(10)/exp(a)),
  #                logD58 = logD0 - (58)/z,
  #                D58 = 10^logD58
  #     )
  # ) %>%
  map2(., map(dyna_fit_r10, coef),
       ~ mutate(.x, logD_dyna = .y[["D_logref"]], logC0_dyna =.y[["logC0"]])
  ) %>%
  imap_dfr(., ~ mutate(.x, strain = .y)) %>%
  mutate(
    C0 = 10^logC0,
    C0_dyna = 10^logC0_dyna,
    D58_dyna = 10^logD_dyna,
    D58_iso = 10^logDref
  ) %>%
  mutate(rel_D58 = D58_dyna/D58_iso,
         relC0 = logC0_dyna - logC0) 

## Plot Ramp 14

this_d <- d_dyna %>%
  filter(Ramp == 14)

logN0 <- this_d %>% filter(time == 0) %>% summarize(l = mean(logN)) %>% pull(l)

times <- seq(0, max(this_d$time), length = 1000)

env_conditions <- tibble(time = c(0, (58-30)/14, max(d$time)),
                         temp = c(30, 58, 58))

pred_geeraerd_r14 <- global_models %>%
  map(coef) %>%
  map(.,
      ~ list(
        list(par = "logD",
             model = "Bigelow",
             temp = c("xref" = Tref, z = .[["z"]]),
             ref = .[["logDref"]]
        )
      )
  ) %>%
  imap(.,
       ~ predict_inactivation(times,
                              list(model = "Geeraerd_noTail", logN0 = logN0,
                                   C0 = my_COs[[.y]]
                                   # C0 = my_C0
                              ),
                              environment = "dynamic",
                              .x,
                              env_conditions)
  )

pred_geeraerd_r14 %>%
  map(., ~.$simulation) %>%
  imap_dfr(., ~ mutate(.x, strain_name = .y)) %>%
  ggplot() +
  geom_line(aes(x = time, y = logN, colour = strain_name)) +
  geom_point(aes(x = time, y = logN, colour = strain_name),
             data = this_d) +
  facet_wrap("strain_name") +
  coord_cartesian(ylim = c(0, 5.5))

dyna_fit_r14 <- this_d %>%
  split(.$strain_name) %>%
  map(., ~ select(., time, logN)) %>%
  map2(., my_guess,
       ~ fit_inactivation("dynamic", 
                          .x, 
                          "Geeraerd",
                          .y$guess, 
                          .y$known, 
                          secondary_models = sec_models_dyna,
                          env_conditions = env_conditions,
                          algorithm = "regression"
       )
  )

dyna_fit_r14 %>% map(plot) %>% plot_grid(plotlist = .)
dyna_fit_r14 %>% map(summary)
dyna_fit_r14 %>% map(coef)

#- Parameter comparison

global_models %>% 
  map(coef) %>%
  map(., ~ tibble(par = names(.), value = .)) %>%
  map(., ~ pivot_wider(., names_from = par, values_from = value)) %>%
  # imap_dfr(., ~ mutate(.x, strain = .y)) %>%
  # map(.,
  #     ~   mutate(.,
  #                b = log(10)/z,
  #                logD0 = log10(log(10)/exp(a)),
  #                logD58 = logD0 - (58)/z,
  #                D58 = 10^logD58
  #     )
  # ) %>%
  map2(., map(dyna_fit_r14, coef),
       ~ mutate(.x, logD_dyna = .y[["D_logref"]], logC0_dyna =.y[["logC0"]])
  ) %>%
  imap_dfr(., ~ mutate(.x, strain = .y)) %>%
  mutate(
    C0 = 10^logC0,
    C0_dyna = 10^logC0_dyna,
    D58_dyna = 10^logD_dyna,
    D58_iso = 10^logDref
  ) %>%
  mutate(rel_D58 = D58_dyna/D58_iso,
         relC0 = logC0_dyna - logC0) 

## Table 5

list(
  `2ºC/min` = dyna_fit,
  `5ºC/min` = dyna_fit_r5,
  `10ºC/min` = dyna_fit_r10,
  `14ºC/min` = dyna_fit_r14
) %>%
  map(.,
      ~ map(., ~ summary(.)$par)
  ) %>%
  map(.,
      ~ map(., ~ as_tibble(., rownames = "par"))
  ) %>%
  map(.,
      ~ imap_dfr(., ~ mutate(.x, strain = .y))
  ) %>%
  imap_dfr(., ~ mutate(.x, heating = .y)) %>%
  mutate(., 
         Estimate = round(Estimate, 2),
         `Std. Error` = round(`Std. Error`, 2)
  ) %>%
  unite(., value, Estimate, `Std. Error`) %>%
  select(par, strain, heating, value) %>%
  pivot_wider(names_from = par, values_from = value) %>%
  write_excel_csv2(., file = "Table_5.csv")

list(
  `2ºC/min` = dyna_fit,
  `5ºC/min` = dyna_fit_r5,
  `10ºC/min` = dyna_fit_r10,
  `14ºC/min` = dyna_fit_r14
) %>% 
  map(., ~ map(., ~ tibble(r = residuals(.)))) %>%
  map(., ~ imap_dfr(., ~ mutate(.x, strain = .y))) %>%
  map(., ~ mutate(., r2 = r^2)) %>%
  map(., ~ summarize(., ME = mean(r), MSE = mean(r2), .by = strain
  )
  ) %>%
  map(., ~ mutate(., RMSE = sqrt(MSE)))

## Parameter estimates

list(
  `2ºC/min` = dyna_fit,
  `5ºC/min` = dyna_fit_r5,
  `10ºC/min` = dyna_fit_r10,
  `14ºC/min` = dyna_fit_r14
) %>%
  map(.,
      ~ map(., ~ summary(.)$par)
  ) %>%
  map(.,
      ~ map(., ~ as_tibble(., rownames = "par"))
  ) %>%
  map(.,
      ~ imap_dfr(., ~ mutate(.x, strain = .y))
  ) %>%
  imap_dfr(., ~ mutate(.x, heating = .y)) %>%
  # mutate(., 
  #        Estimate = round(Estimate, 2),
  #        `Std. Error` = round(`Std. Error`, 2)
  # ) %>%
  # unite(., value, Estimate, `Std. Error`) %>%
  select(par, strain, heating, Estimate, `Std. Error`) %>%
  mutate(heating = factor(heating, levels = c("2ºC/min", "5ºC/min", "10ºC/min", "14ºC/min"))) %>%
  ggplot() +
  geom_col(aes(x = strain, y = Estimate, fill = heating), position = "dodge") +
  geom_errorbar(aes(x = strain, ymin = Estimate, ymax = Estimate + `Std. Error`, colour = heating), position = "dodge") +
  facet_wrap("par", scales = "free")

## Supp. Figure 3

p <- dyna_fit %>%
  map(plot) %>%
  map2(., pred_geeraerd,
       ~ .x + 
         geom_line(aes(x = time , y = logN), data = .y$simulation, linetype = 2) + 
         coord_cartesian(ylim = c(0, 5.5)) +
         labs(x = "Time (min)", y = "Concentration (log CFU/mL)")
  ) %>%
  plot_grid(plotlist = ., labels = "AUTO")

p
ggsave(p, filename = "supp_Figure_3.png", width = 12, height = 9, bg = "white")

## Supp. Figure 4

p <- dyna_fit_r5 %>%
  map(plot) %>%
  map2(., pred_geeraerd_r5,
       ~ .x + 
         geom_line(aes(x = time , y = logN), data = .y$simulation, linetype = 2) + 
         coord_cartesian(ylim = c(0, 5.5)) +
         labs(x = "Time (min)", y = "Concentration (log CFU/mL)")
  ) %>%
  plot_grid(plotlist = ., labels = "AUTO")

p
ggsave(p, filename = "supp_Figure_4.png", width = 12, height = 9, bg = "white")

## Supp. Figure 5

p <- dyna_fit_r10 %>%
  map(plot) %>%
  map2(., pred_geeraerd_r10,
       ~ .x + 
         geom_line(aes(x = time , y = logN), data = .y$simulation, linetype = 2) + 
         coord_cartesian(ylim = c(0, 5.5)) +
         labs(x = "Time (min)", y = "Concentration (log CFU/mL)")
  ) %>%
  plot_grid(plotlist = ., labels = "AUTO")

p
ggsave(p, filename = "supp_Figure_5.png", width = 12, height = 9, bg = "white")

## Supp. Figure 6

p <- dyna_fit_r14 %>%
  map(plot) %>%
  map2(., pred_geeraerd_r14,
       ~ .x + 
         geom_line(aes(x = time , y = logN), data = .y$simulation, linetype = 2) + 
         coord_cartesian(ylim = c(0, 5.5)) +
         labs(x = "Time (min)", y = "Concentration (log CFU/mL)")
  ) %>%
  plot_grid(plotlist = ., labels = "AUTO")

p
ggsave(p, filename = "supp_Figure_6.png", width = 12, height = 9, bg = "white")

## Table 4

pred_geeraerd %>% 
  map(., ~ .$simulation) %>%
  map(., ~ select(., time, logN)) %>%
  map(as.data.frame) %>%
  map(.,
      ~ modCost(model = .,
                obs = d_dyna %>% filter(Ramp == 2) %>% select(time, logN) %>% as.data.frame()
      )$res
  ) %>%
  map(., ~ tibble(r = .$res)) %>%
  map(., ~ summarize(., ME = mean(r)))

pred_geeraerd_r5 %>% 
  map(., ~ .$simulation) %>%
  map(., ~ select(., time, logN)) %>%
  map(as.data.frame) %>%
  map(.,
      ~ modCost(model = .,
                obs = d_dyna %>% filter(Ramp == 5) %>% select(time, logN) %>% as.data.frame()
      )$res
  ) %>%
  map(., ~ tibble(r = .$res)) %>%
  map(., ~ summarize(., ME = mean(r)))

pred_geeraerd_r10 %>% 
  map(., ~ .$simulation) %>%
  map(., ~ select(., time, logN)) %>%
  map(as.data.frame) %>%
  map(.,
      ~ modCost(model = .,
                obs = d_dyna %>% filter(Ramp == 10) %>% select(time, logN) %>% as.data.frame()
      )$res
  ) %>%
  map(., ~ tibble(r = .$res)) %>%
  map(., ~ summarize(., ME = mean(r)))

pred_geeraerd_r14 %>% 
  map(., ~ .$simulation) %>%
  map(., ~ select(., time, logN)) %>%
  map(as.data.frame) %>%
  map(.,
      ~ modCost(model = .,
                obs = d_dyna %>% filter(Ramp == 14) %>% select(time, logN) %>% as.data.frame()
      )$res
  ) %>%
  map(., ~ tibble(r = .$res)) %>%
  map(., ~ summarize(., ME = mean(r)))

## Figure 3A

pars_expo <- global_models %>%
  map(., ~ summary(.)$par) %>%
  map(., ~ as_tibble(., rownames = "par")) %>%
  imap_dfr(., ~ mutate(.x, strain = .y)) %>%
  mutate(heating = "Iso - exponential")

pars_stat <- global_models_stat %>%
  map(., ~ summary(.)$par) %>%
  map(., ~ as_tibble(., rownames = "par")) %>%
  imap_dfr(., ~ mutate(.x, strain = .y)) %>%
  mutate(heating = "Iso - stationary")


pars_dynamic <- list(
  `2ºC/min` = dyna_fit,
  `5ºC/min` = dyna_fit_r5,
  `10ºC/min` = dyna_fit_r10,
  `14ºC/min` = dyna_fit_r14
) %>%
  map(.,
      ~ map(., ~ summary(.)$par)
  ) %>%
  map(.,
      ~ map(., ~ as_tibble(., rownames = "par"))
  ) %>%
  map(.,
      ~ imap_dfr(., ~ mutate(.x, strain = .y))
  ) %>%
  imap_dfr(., ~ mutate(.x, heating = .y)) %>%
  # mutate(., 
  #        Estimate = round(Estimate, 2),
  #        `Std. Error` = round(`Std. Error`, 2)
  # ) %>%
  # unite(., value, Estimate, `Std. Error`) %>%
  select(par, strain, heating, Estimate, `Std. Error`)

p1 <- bind_rows(pars_expo, pars_stat, pars_dynamic) %>%
  filter(par %in% c("D_logref", "logDref")) %>% 
  mutate(heating = factor(heating,
                          levels = c("Iso - exponential",
                                     "Iso - stationary",
                                     "2ºC/min",
                                     "5ºC/min",
                                     "10ºC/min",
                                     "14ºC/min")
                          )) %>%
  ggplot(aes(x = strain, y = Estimate, fill = heating)) +
  geom_col(position = "dodge") +
  geom_errorbar(aes(ymin = Estimate, ymax = Estimate + `Std. Error`,
                    colour = heating),
                position = "dodge") +
  labs(x = "", y = "log D at 58ºC (log min)") +
  ggthemes::theme_few(base_size = 14) +
  scale_color_npg() +
  scale_fill_npg() +
  theme(legend.position = "none",
        axis.text = element_text(angle = 30, hjust = 1)
        )
## Figure 3B

p2 <- bind_rows(pars_expo, pars_stat, pars_dynamic) %>%
  filter(par %in% c("logC0", "D_logref", "logDref")) %>% 
  select(strain, heating, par, Estimate) %>%
  mutate(par = ifelse(par == "D_logref", "logDref", par)) %>%
  pivot_wider(names_from = par, values_from = Estimate) %>%
  mutate(C0 = 10^logC0, 
         Dref = 10^logDref,
         kref = log(10)/Dref,
         SL = log(C0 + 1)/kref
         ) %>%
  mutate(heating = factor(heating,
                          levels = c("Iso - exponential",
                                     "Iso - stationary",
                                     "2ºC/min",
                                     "5ºC/min",
                                     "10ºC/min",
                                     "14ºC/min")
  )) %>%
  ggplot(aes(x = strain, y = SL, fill = heating)) +
  geom_col(position = "dodge") +
  labs(x = "", y = "Shoulder length at 58ºC (min)") +
  ggthemes::theme_few(base_size = 14) +
  scale_color_npg() +
  scale_fill_npg() +
  theme(legend.title = element_blank(),
        legend.position = "right",
        axis.text = element_text(angle = 30, hjust = 1)
        )

p <- plot_grid(p1, p2, labels = "AUTO",
               rel_widths = c(.4, .6))
p
ggsave(p, filename = "Figure_3.png", width = 12, height = 4)




















