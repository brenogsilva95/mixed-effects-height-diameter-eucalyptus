## Packages
{library(ggplot2);library(dplyr);library(lme4);library(nlme);
  library(lmerTest);library(cAIC4);library(splines);library(gamlss);
  library(gridExtra);library(hnp);library(varTestnlme);library(car);
  library(Rcpp);library(MASS);library(Matrix);library(lmer);library(tidyverse);
  library(merDeriv);library(qqplotr); library(DHARMa)}

## Reading the dataset 
Eucalyptus_dat <- read.table("df_thesis.csv", 
                             header = T, sep=";", dec = ",") %>% 
  rename(cbh = cap,
         age = idade,
         htotal = altura,
         activity = atividade,
         stand = talhao) %>% 
  mutate(dbh = cbh / (pi * 10),
         h = htotal / 10,
         dg_dbh = dg / dbh,
         inv_dbh_age = 1 / (dbh * age),
         inv_dbh = 1 / dbh) %>%
  filter(h > 0, cod1 != "Q") %>%
  mutate(cod1 = as.factor(cod1),
         rf = as.factor(rf),
         stand = as.factor(stand)) %>%
  dplyr::select(h, dbh, hdom, dg, dg_dbh, age, activity,
                inv_dbh_age, inv_dbh, stand, rf)
Eucalyptus_dat$Arvore_key <- 1:nrow(Eucalyptus_dat)
Eucalyptus_dat$Arvore_key = as.factor(paste(Eucalyptus_dat$rf, Eucalyptus_dat$talhao, Eucalyptus_dat$Arvore_key, SEP="_"))
## Exploratory analysis 
############################################################################
##################### Table: descriptive analysis ##########################
############################################################################
Table <- Eucalyptus_dat %>%
  group_by(Classe_Idade = cut(age, breaks = seq(1, 10, by = 2), right = FALSE, include.lowest = TRUE)) %>%
  summarise(
    n = n_distinct(Arvore_key),          
    d_mean = mean(dbh, na.rm = TRUE), 
    d_max = max(dbh, na.rm = TRUE),   
    d_min = min(dbh, na.rm = TRUE),   
    h_mean = mean(h, na.rm = TRUE),    
    h_max = max(h, na.rm = TRUE),     
    h_min = min(h, na.rm = TRUE)       
  )
print(Table)
############################################################################
##################### Table: descriptive analysis ##########################
############################################################################
classes_d <- c(1, 5, 10, 15, 20, 25, 30, 35, 40)
classes_h <- c(1, 5, 10, 15, 20, 25, 30, 35, 40)
Eucalyptus_dat <- Eucalyptus_dat %>%
  mutate(
    d_class = cut(dbh, breaks = classes_d, right = FALSE, include.lowest = TRUE),
    h_class = cut(h, breaks = classes_h, right = FALSE, include.lowest = TRUE)
  )
todas_classes <- expand.grid(
  d_class = levels(cut(Eucalyptus_dat$dbh, breaks = classes_d, right = FALSE, include.lowest = TRUE)),
  h_class = levels(cut(Eucalyptus_dat$h, breaks = classes_h, right = FALSE, include.lowest = TRUE))
)
tabela_frequencia <- Eucalyptus_dat %>%
  distinct(Arvore_key, d_class, h_class) %>%  
  count(d_class, h_class, name = "num_arvores") %>%
  right_join(todas_classes, by = c("d_class", "h_class")) %>%
  replace_na(list(num_arvores = 0)) %>%
  pivot_wider(names_from = h_class, values_from = num_arvores, values_fill = 0) %>%
  arrange(d_class)
tabela_frequencia <- tabela_frequencia %>%
  mutate(Total = rowSums(across(-d_class))) %>%
  bind_rows(
    summarise(., d_class = "Total", across(-d_class, ~ sum(.x, na.rm = TRUE)))
  )
tabela_frequencia
############################################################################
##################### Graphs: descriptive analysis #########################
############################################################################
Eucalyptus_dat <- Eucalyptus_dat %>%
  mutate(Classe_age = cut(age, breaks = seq(1, 10, by = 2), right = FALSE, include.lowest = TRUE))
boxplot_d <- ggplot(Eucalyptus_dat, aes(x = Classe_age, y = dbh, fill = Classe_age)) +
  geom_boxplot() +
  labs(x = "Age Classes (years)", y = "D (cm)") +
  theme(
    axis.title.x = element_text(size = 18, color = "black"),
    axis.title.y = element_text(size = 18, color = "black"),
    axis.text.x = element_text(size = 20, angle = 45, color = "black", hjust = 1),
    axis.text.y = element_text(size = 20, color = "black"),
    plot.title = element_text(size = 45),
    legend.position = "none"
  )
violin_d <- ggplot(Eucalyptus_dat, aes(x = Classe_age, y = dbh, fill = Classe_age)) +
  geom_violin() +
  geom_boxplot(width = 0.1, color = "black", alpha = 0.7) + 
  labs(x = "Age Classes (years)", y = "D (cm)") +
  theme(
    axis.title.x = element_text(size = 18, color = "black"),
    axis.title.y = element_text(size = 18, color = "black"),
    axis.text.x = element_text(size = 20, angle = 45, color = "black", hjust = 1),
    axis.text.y = element_text(size = 20, color = "black"),
    plot.title = element_text(size = 45),
    legend.position = "none"
  )
boxplot_h <- ggplot(Eucalyptus_dat, aes(x = Classe_age, y = h, fill = Classe_age)) +
  geom_boxplot() +
  labs(x = "Age Classes (years)", y = "Total height (m)") +
  theme(
    axis.title.x = element_text(size = 18, color = "black"),
    axis.title.y = element_text(size = 18, color = "black"),
    axis.text.x = element_text(size = 20, angle = 45, color = "black", hjust = 1),
    axis.text.y = element_text(size = 20, color = "black"),
    plot.title = element_text(size = 45),
    legend.position = "none"
  )
violin_h <- ggplot(Eucalyptus_dat, aes(x = Classe_age, y = h, fill = Classe_age)) +
  geom_violin() +
  geom_boxplot(width = 0.1, color = "black", alpha = 0.7) + 
  labs(x = "Age Classes (years)", y = "Total height (m)") +
  theme(
    axis.title.x = element_text(size = 18, color = "black"),
    axis.title.y = element_text(size = 18, color = "black"),
    axis.text.x = element_text(size = 20, angle = 45, color = "black", hjust = 1),
    axis.text.y = element_text(size = 20, color = "black"),
    plot.title = element_text(size = 45),
    legend.position = "none"
  )
grid.arrange(boxplot_d, violin_d, boxplot_h, violin_h, ncol = 2)
## Dispersion graphs
## h vs hdom
ggplot(Eucalyptus_dat,
       aes(x = hdom, 
           y = h)) +
  geom_point(size = 1.2)+
  labs(x = 'Hd (m)', y = 'H (m)') +
  geom_smooth(se=FALSE)+
  theme(legend.position = "none",
        axis.title = element_text(size = 40),
        axis.text.x = element_text(color = "black", hjust=1),
        axis.text.y = element_text(color = "black", hjust=1),
        axis.text = element_text(size = 40),
        plot.title = element_text(size = 40),
        strip.text.x = element_text(size = 12))
## log(h) vs log(HD) 
ggplot(Eucalyptus_dat,
       aes(x = log(hdom), 
           y = log(h))) +
  geom_point(size = 1.2)+
  labs(x = 'log(Hd)', y = 'log(H)') +
  geom_smooth(se=FALSE)+
  theme(legend.position = "none",
        axis.title = element_text(size = 40),
        axis.text.x = element_text(color = "black", hjust=1),
        axis.text.y = element_text(color = "black", hjust=1),
        axis.text = element_text(size = 40),
        plot.title = element_text(size = 40),
        strip.text.x = element_text(size = 12))
## h vs dg/dbh 
ggplot(Eucalyptus_dat,
       aes(x = dg_dbh, 
           y = h)) +
  geom_point(size = 1.2)+
  labs(x = bquote(frac(Dg, DBH)), y = 'h') +
  geom_smooth(se=FALSE)+
  theme(legend.position = "none",
        axis.title = element_text(size = 40),
        axis.text.x = element_text(color = "black", hjust=1),
        axis.text.y = element_text(color = "black", hjust=1),
        axis.text = element_text(size = 40),
        plot.title = element_text(size = 40),
        strip.text.x = element_text(size = 12))
## log(h) vs log(dg/dbh) 
ggplot(Eucalyptus_dat,
       aes(x = log(dg_dbh), 
           y = log(h))) +
  geom_point(size = 1.2)+
  labs(x = bquote(log(frac(Dg, D))), y = 'log(H)') +
  geom_smooth(se = FALSE) +
  theme(legend.position = "none",
        axis.title = element_text(size = 40),
        axis.text.x = element_text(color = "black", hjust=1),
        axis.text.y = element_text(color = "black", hjust=1),
        axis.text = element_text(size = 40),
        plot.title = element_text(size = 40),
        strip.text.x = element_text(size = 12))
## h vs inv(DBH*age) 
ggplot(Eucalyptus_dat,
       aes(x = inv_dbh_age, 
           y = h)) +
  geom_point(size = 1.2)+
  labs(x = bquote(frac(1, D ~ "*" ~ I)), y = 'h') +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.01)) +
  geom_smooth(se=FALSE) +
  theme(legend.position = "none",
        axis.title = element_text(size = 40),
        axis.text.x = element_text(color = "black", hjust=1),
        axis.text.y = element_text(color = "black", hjust=1),
        axis.text = element_text(size = 40),
        plot.title = element_text(size = 40),
        strip.text.x = element_text(size = 12))
## log(h) vs inv(DBH*age) 
ggplot(Eucalyptus_dat,
       aes(x = inv_dbh_age, 
           y = log(h))) +
  geom_point(size = 1.2)+
  geom_smooth(se=FALSE) +
  labs(x = bquote(frac(1, D ~ "*" ~ I)), y = 'log(H)') +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.01)) +
  theme(legend.position = "none",
        axis.title = element_text(size = 40),
        axis.text.x = element_text(color = "black", hjust=1),
        axis.text.y = element_text(color = "black", hjust=1),
        axis.text = element_text(size = 40),
        plot.title = element_text(size = 40),
        strip.text.x = element_text(size = 12))
## h vs DBH 
ggplot(Eucalyptus_dat,
       aes(x = dbh, 
           y = h)) +
  geom_point(size = 1.2)+
  geom_smooth(se = FALSE) +
  labs(x = "DBH", y = 'h') +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.01)) +
  theme(legend.position = "none",
        axis.title = element_text(size = 40),
        axis.text.x = element_text(color = "black", hjust=1),
        axis.text.y = element_text(color = "black", hjust=1),
        axis.text = element_text(size = 40),
        plot.title = element_text(size = 40),
        strip.text.x = element_text(size = 12))
## log(h) vs log(DBH) 
ggplot(Eucalyptus_dat,
       aes(x = log(dbh), 
           y = log(h))) +
  geom_point(size = 1.2)+
  geom_smooth(se = FALSE) +
  labs(x = "log(DBH)", y = 'log(h)') +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.01)) +
  theme(legend.position = "none",
        axis.title = element_text(size = 40),
        axis.text.x = element_text(color = "black", hjust=1),
        axis.text.y = element_text(color = "black", hjust=1),
        axis.text = element_text(size = 40),
        plot.title = element_text(size = 40),
        strip.text.x = element_text(size = 12))
## h vs inv(DBH)
ggplot(Eucalyptus_dat,
       aes(x = inv_dbh, 
           y = h)) +
  geom_point(size = 1.2)+
  geom_smooth(se = FALSE) +
  labs(x = bquote(frac(1, DBH)), y = 'h') +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.01)) +
  theme(legend.position = "none",
        axis.title = element_text(size = 40),
        axis.text.x = element_text(color = "black", hjust=1),
        axis.text.y = element_text(color = "black", hjust=1),
        axis.text = element_text(size = 40),
        plot.title = element_text(size = 40),
        strip.text.x = element_text(size = 12))
## log(h) vs inv(DBH)
ggplot(Eucalyptus_dat,
       aes(x = inv_dbh, 
           y = log(h))) +
  geom_point(size = 1.2)+
  geom_smooth(se = FALSE) +
  labs(x = bquote(frac(1, D)), y = 'log(H)') +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.01)) +
  theme(legend.position = "none",
        axis.title = element_text(size = 40),
        axis.text.x = element_text(color = "black", hjust=1),
        axis.text.y = element_text(color = "black", hjust=1),
        axis.text = element_text(size = 40),
        plot.title = element_text(size = 40),
        strip.text.x = element_text(size = 12))
## h vs I
ggplot(Eucalyptus_dat,
       aes(x = age, 
           y = h)) +
  geom_point(size = 1.2)+
  geom_smooth(se = FALSE) +
  labs(x = "Age", y = 'h') +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.01)) +
  theme(legend.position = "none",
        axis.title = element_text(size = 40),
        axis.text.x = element_text(color = "black", hjust=1),
        axis.text.y = element_text(color = "black", hjust=1),
        axis.text = element_text(size = 40),
        plot.title = element_text(size = 40),
        strip.text.x = element_text(size = 12))
## log(h) vs log(I)
ggplot(Eucalyptus_dat,
       aes(x = log(age), 
           y = log(h))) +
  geom_point(size = 1.2)+
  geom_smooth(se = FALSE) +
  labs(x = "log(I)", y = 'log(H)') +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.01)) +
  theme(legend.position = "none",
        axis.title = element_text(size = 40),
        axis.text.x = element_text(color = "black", hjust=1),
        axis.text.y = element_text(color = "black", hjust=1),
        axis.text = element_text(size = 40),
        plot.title = element_text(size = 40),
        strip.text.x = element_text(size = 12))
## log(h) vs 1/I
ggplot(Eucalyptus_dat,
       aes(x = 1/age, 
           y = log(h))) +
  geom_point(size = 1.2)+
  geom_smooth(se = FALSE) +
  labs(x = bquote(frac(1, Age)), y = 'log(h)') +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.01)) +
  theme(legend.position = "none",
        axis.title = element_text(size = 40),
        axis.text.x = element_text(color = "black", hjust=1),
        axis.text.y = element_text(color = "black", hjust=1),
        axis.text = element_text(size = 40),
        plot.title = element_text(size = 40),
        strip.text.x = element_text(size = 12))
## DBH vs I
ggplot(Eucalyptus_dat,
       aes(x = age, 
           y = dbh)) +
  geom_point(size = 1.2)+
  geom_smooth(se = FALSE) +
  labs(x = "Age", y = 'DBH') +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.01)) +
  theme(legend.position = "none",
        axis.title = element_text(size = 40),
        axis.text.x = element_text(color = "black", hjust=1),
        axis.text.y = element_text(color = "black", hjust=1),
        axis.text = element_text(size = 40),
        plot.title = element_text(size = 40),
        strip.text.x = element_text(size = 12))
## log(DBH) vs log(I)
ggplot(Eucalyptus_dat,
       aes(x = log(age), 
           y = log(dbh))) +
  geom_point(size = 1.2)+
  geom_smooth(se = FALSE) +
  labs(x = "log(I)", y = 'log(D)') +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.01)) +
  theme(legend.position = "none",
        axis.title = element_text(size = 40),
        axis.text.x = element_text(color = "black", hjust=1),
        axis.text.y = element_text(color = "black", hjust=1),
        axis.text = element_text(size = 40),
        plot.title = element_text(size = 40),
        strip.text.x = element_text(size = 12))
## log(DBH) vs log(I)
ggplot(Eucalyptus_dat,
       aes(x = 1/age, 
           y = 1/dbh)) +
  geom_point(size = 1.2)+
  geom_smooth(se = FALSE) +
  labs(x = bquote(frac(1, Age)), y = bquote(frac(1, DBH))) +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.01)) +
  theme(legend.position = "none",
        axis.title = element_text(size = 40),
        axis.text.x = element_text(color = "black", hjust=1),
        axis.text.y = element_text(color = "black", hjust=1),
        axis.text = element_text(size = 40),
        plot.title = element_text(size = 40),
        strip.text.x = element_text(size = 12))

ggplot(Eucalyptus_dat,
       aes(x = dbh, 
           y = dg)) +
  geom_point(size = 1.2)+
  geom_smooth(se = FALSE) +
  labs(x = "DBH", y = "DG") +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.01)) +
  theme(legend.position = "none",
        axis.title = element_text(size = 40),
        axis.text.x = element_text(color = "black", hjust=1),
        axis.text.y = element_text(color = "black", hjust=1),
        axis.text = element_text(size = 40),
        plot.title = element_text(size = 40),
        strip.text.x = element_text(size = 12))


ggplot(Eucalyptus_dat,
       aes(x = log(dbh), 
           y = log(dg))) +
  geom_point(size = 1.2)+
  geom_smooth(se = FALSE) +
  labs(x = "log(DBH)", y = "log(DG)") +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.01)) +
  theme(legend.position = "none",
        axis.title = element_text(size = 40),
        axis.text.x = element_text(color = "black", hjust=1),
        axis.text.y = element_text(color = "black", hjust=1),
        axis.text = element_text(size = 40),
        plot.title = element_text(size = 40),
        strip.text.x = element_text(size = 12))

Table <- Eucalyptus_dat %>%
  mutate(
    Classe_Idade = cut(age, breaks = seq(1, 10, by = 2), right = FALSE, include.lowest = TRUE),
    Classe_Idade_Label = paste0("Age Class (years)\n", levels(Classe_Idade))[Classe_Idade]
  )

# Criação do gráfico com ggplot2
ggplot(Table, aes(x = dbh, y = h)) +
  geom_point() +  # Adiciona os pontos
  facet_wrap(~ Classe_Idade_Label, ncol = 4) +  
  labs(
    x = "D (cm)", 
    y = "H (m)", 
    title = ""
  ) +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 40),
    axis.text.x = element_text(color = "black", hjust = 1),
    axis.text.y = element_text(color = "black", hjust = 1),
    axis.text = element_text(size = 40),
    plot.title = element_text(size = 40),
    strip.text = element_text(size = 20, face = "bold")  
  )

############################################################################
############# Descriptive analysis: number of stands per FR ################
############################################################################
talhoes_por_RF <- Eucalyptus_dat %>%
  group_by(rf) %>%
  summarize(n_talhoes = n_distinct(stand)) %>%
  mutate(rf = paste0("R", row_number()))
print(talhoes_por_RF, n = 30)
table(talhoes_por_RF$n_talhoes)
talhoes_por_RF <- data.frame(
  rf = paste0("R", 1:30),
  n_talhoes = talhoes_por_RF$n_talhoes
)
talhoes_por_RF$rf <- factor(talhoes_por_RF$rf, levels = paste0("R", 1:30))
ggplot(talhoes_por_RF, aes(x = rf, y = n_talhoes)) +
  geom_bar(stat = "identity") +
  coord_cartesian(ylim = c(0,5))+
  geom_text(aes(label = n_talhoes), vjust = -0.5, size = 8) + 
  scale_fill_manual(values = 'gray') +
  labs(x = "Forest Regions", y = "Number of Stands") +
  theme(
    axis.title.x = element_text(size = 30, color = "black"),
    axis.title.y = element_text(size = 30, color = "black"),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 30, color = "black"),
    plot.title = element_text(size = 45),
    legend.position = "none"
  )

result <- talhoes_por_RF %>%
  summarize(
    min_talhao = min(n_talhoes),
    max_talhao = max(n_talhoes),
    mean_talhao = mean(n_talhoes),
    var_talhao = var(n_talhoes)
  ) 
print(result)
sum(talhoes_por_RF$n_talhoes)
############################################################################
arvores_por_talhoes <- Eucalyptus_dat %>%
  group_by(stand) %>%
  summarize(n_arvores = n_distinct(Arvore_key)) 
print(arvores_por_talhoes, n = 45)
sum(arvores_por_talhoes$n_arvores)
table(arvores_por_talhoes$n_arvores)
arvores_por_talhoes <- data.frame(
  S = paste0("S", 1:45),
  n_arvores = arvores_por_talhoes$n_arvores
)
arvores_por_talhoes$S<- factor(arvores_por_talhoes$S, levels = paste0("S", 1:45))
ggplot(arvores_por_talhoes, aes(x = S, y = n_arvores)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = n_arvores), vjust = -0.5, size = 5) + 
  scale_fill_manual(values = 'gray') +
  labs(x = "Stands", y = "Number of Trees") +
  theme(
    axis.title.x = element_text(size = 30, color = "black"),
    axis.title.y = element_text(size = 30, color = "black"),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 30, color = "black"),
    plot.title = element_text(size = 45),
    legend.position = "none"
  )

result_arv <- arvores_por_talhoes %>%
  summarize(
    min_talhao = min(n_arvores),
    max_talhao = max(n_arvores),
    mean_talhao = mean(n_arvores),
    var_talhao = var(n_arvores)
  ) 
print(result_arv)
## log(h) vs (log(Dg) + log(Dg/DBH)+ 1/DBH)
## extract residuals - residuals1

## log(1/(DBH*I)) vs (log(Dg) + log(Dg/DBH)+ 1/DBH)
## extract residuals - residuals2
## graphic of res1 vs res2

#####
## log(h) vs (log(Dg) + log(Dg/DBH)+ log(1/(DBH*I)))
## extract residuals - residuals1

## 1/DBH vs (log(Dg) + log(Dg/DBH)+ log(1/(DBH*I)))
## extract residuals - residuals2
## graphic of res1 vs res2
# Case 1:
Mod1 <- lm(log(h)~log(dg)+log(dg_dbh)+inv_dbh,
           data = Eucalyptus_dat) 
resid_1 <- residuals(Mod1)
# Case 2:
Mod2 <- lm(log(inv_dbh_age)~log(dg)+log(dg_dbh)+inv_dbh,
           data = Eucalyptus_dat) 
resid_2 <- residuals(Mod2)
# Residuals plot
resid_dataframe1 <- data.frame(resid_1 = resid_1, resid_2 = resid_2)
ggplot(resid_dataframe1, aes(x = resid_2, y = resid_1)) +
  geom_point() +
  geom_smooth(se=FALSE)+
  labs(x = "Mod2 Residuals",
       y = "Mod1 Residuals") +
  theme_bw() +
  theme(legend.position = "none",
        axis.title = element_text(size = 25),
        axis.text.x = element_text(color = "black", hjust=1),
        axis.text.y = element_text(color = "black", hjust=1),
        axis.text = element_text(size = 25),
        plot.title = element_text(size = 25),
        strip.text.x = element_text(size = 12))
# Case 3:
Mod3 <- lm(log(h)~log(dg)+log(dg_dbh)+log(inv_dbh_age),
           data = Eucalyptus_dat) 
resid_3 <- residuals(Mod3)
# Case 4:
Mod4 <- lm(inv_dbh~log(dg)+log(dg_dbh)+log(inv_dbh_age),
           data = Eucalyptus_dat)
resid_4 <- residuals(Mod4)
# Residuals plot
resid_dataframe2 <- data.frame(resid_3 = resid_3, resid_4 = resid_4)
ggplot(resid_dataframe2, aes(x = resid_4, y = resid_3)) +
  geom_point() +
  geom_smooth(se=FALSE)+
  labs(x = "Mod4 Residuals",
       y = "Mod3 Residuals") +
  theme_bw() +
  theme(legend.position = "none",
        axis.title = element_text(size = 25),
        axis.text.x = element_text(color = "black", hjust=1),
        axis.text.y = element_text(color = "black", hjust=1),
        axis.text = element_text(size = 25),
        plot.title = element_text(size = 25),
        strip.text.x = element_text(size = 12))

## Variance Component Tests 

## H0: no random effects 
## H1: sigma^2_stands 
## i.e. sigma^2_stands = 0 vs. sigma^2_stands > 0 
mod.1 <- lm(log(h)~log(hdom)+log(dg_dbh)+
                        (inv_dbh_age)+(inv_dbh),
                      Eucalyptus_dat) # without random effects
mod.2_st <- lmer(log(h)~log(hdom)+log(dg_dbh)+
                          (inv_dbh_age)+(inv_dbh) + (1|stand),
                        Eucalyptus_dat, REML = F) # random effects per stands
ranova(mod.2_st) 
2*(653.04527004179-612.83302429934) #= 80.42449
1-0.5*pchisq(80.4244914849,0)-0.5*pchisq(80.4244914849,1) # p-value: 0
## varTestnlme package
t1 = varCompTest(mod.2_st, mod.1);summary(t1)

## H0: sigma^2_stands 
## H1: sigma^2_stands & sigma^2_forestregion 
##  i.e. sigma^2_forestregion = 0 vs. sigma^2_forestregion > 0 
mod.2_st <- lmer(log(h)~log(hdom)+log(dg_dbh)+
                          (inv_dbh_age)+(inv_dbh) + (1|stand),
                        Eucalyptus_dat, REML = F)
mod.4_both <- lmer(log(h)~log(hdom)+log(dg_dbh)+
                        (inv_dbh_age)+(inv_dbh) + (1|stand) + (1|rf),
                      Eucalyptus_dat, REML = F) # random effects per stands & forestry region
ranova(mod.4_both) 
2*(655.11614773686-653.04527004179)
1-0.5*pchisq(4.1417553901399,0)-0.5*pchisq(4.1417553901399,1) # p-value: 0.020918881319204
## varTestnlme package
t2 = varCompTest(mod.4_both, mod.2_st);summary(t2)
# Variance components testing in mixed effects models
# Error in extractStruct.merMod(m1, m0, randm0) : 
#   the package does not currently support more than one level of random effects

## H0: no random effects 
## H1: sigma^2_rf 
## i.e. sigma^2_rf = 0 vs. sigma^2_rfs > 0 
mod.1 <- lm(log(h)~log(hdom)+log(dg_dbh)+
                        (inv_dbh_age)+(inv_dbh),
                      Eucalyptus_dat)
mod.3_fr <- lmer(log(h)~log(hdom)+log(dg_dbh)+
                      (inv_dbh_age)+(inv_dbh) + (1|rf),
                    Eucalyptus_dat, REML = F)
ranova(mod.3_fr) 
2*(646.87-612.83302429934) #= 68.07395
1-0.5*pchisq(68.07395,0)-0.5*pchisq(68.07395,1) # p-value: 0

## H0: sigma^2_forestregion
## H1: sigma^2_stands & sigma^2_forestregion 
##  i.e. sigma^2_forestregion = 0 vs. sigma^2_forestregion > 0 
mod.3_fr <- lmer(log(h)~log(hdom)+log(dg_dbh)+
                      (inv_dbh_age)+(inv_dbh) + (1|rf),
                    Eucalyptus_dat, REML = F)
mod.4_both <- lmer(log(h)~log(hdom)+log(dg_dbh)+
                        (inv_dbh_age)+(inv_dbh) + (1|stand) + (1|rf),
                      Eucalyptus_dat, REML = F)
ranova(mod.4_both) 
2*(655.11614773686-646.87)
1-0.5*pchisq(16.4923,0)-0.5*pchisq(16.4923,1) # p-value: 0

## H0: no random effects
## H1: sigma^2_stands & sigma^2_forestregion 
## i.e. sigma^2_forestregion & sigma^2_forestregion = 0 vs. > 0
mod.1 <- lm(log(h)~log(hdom)+log(dg_dbh)+
                        (inv_dbh_age)+(inv_dbh),
                      Eucalyptus_dat)
mod.4_both <- lmer(log(h)~log(hdom)+log(dg_dbh)+
                        (inv_dbh_age)+(inv_dbh) + (1|stand) + (1|rf),
                      Eucalyptus_dat, REML = F)
2*(655.11614773686-612.83302429934) #= 84.566246875
-1225.6660485987 - (-1310.2322954737) #= 84.566246875
1-0.25*pchisq(84.566246875,0)-0.5*pchisq(84.566246875,1)-0.25*pchisq(84.566246875,2) # p-value: 0

## Initial Diagnostic Analysis 
mod.1_ <- lm(log(h)~log(hdom)+log(dg_dbh)+
                       (inv_dbh_age)+(inv_dbh),
                     Eucalyptus_dat)
mod.2_st_ <- lmer(log(h)~log(hdom)+log(dg_dbh)+
                     (inv_dbh_age)+(inv_dbh) + (1|stand),
                   Eucalyptus_dat, REML = T)
mod.3_fr_ <- lmer(log(h)~log(hdom)+log(dg_dbh)+
                 (inv_dbh_age)+(inv_dbh) + (1|rf),
               Eucalyptus_dat, REML = T)
mod.4_both_ <- lmer(log(h)~log(hdom)+log(dg_dbh)+
                   (inv_dbh_age)+(inv_dbh) + (1|stand) + (1|rf),
                 Eucalyptus_dat, REML = T)

## AIC
AIC(mod.1_)
AIC(mod.2_st_)
AIC(mod.3_fr_)
AIC(mod.4_both_)
## cAIC
cAIC(mod.2_st_)
cAIC(mod.3_fr_)
cAIC(mod.4_both_)
## -2log(L)
-2*logLik(mod.1_)
-2*logLik(mod.2_st_)
-2*logLik(mod.3_fr_)
-2*logLik(mod.4_both_)

mod.1_with_lm <- lm(log(h)~log(hdom)+log(dg_dbh)+
                                    (inv_dbh_age)+(inv_dbh),
                                  Eucalyptus_dat)

mod.2_rf.st.tree_with_REML <- lme(log(h)~log(hdom)+log(dg_dbh)+
                                       (inv_dbh_age)+(inv_dbh), random = ~1|rf/stand/Arvore_key,
                                     Eucalyptus_dat, method="REML")
mod.3_st.tree_with_REML <- lme(log(h)~log(hdom)+log(dg_dbh)+
                                    (inv_dbh_age)+(inv_dbh), random = ~1|stand/Arvore_key,
                                  Eucalyptus_dat, method="REML")
mod.4_rf.tree_with_REML <- lme(log(h)~log(hdom)+log(dg_dbh)+
                                    (inv_dbh_age)+(inv_dbh), random = ~1|rf/Arvore_key,
                                  Eucalyptus_dat, method="REML")
mod.5_rf_with_REML <- lme(log(h)~log(hdom)+log(dg_dbh)+
                               (inv_dbh_age)+(inv_dbh), random = ~1|rf,
                             Eucalyptus_dat, method="REML")
mod.6_stand_with_REML <- lme(log(h)~log(hdom)+log(dg_dbh)+
                                  (inv_dbh_age)+(inv_dbh), random = ~1|stand,
                                Eucalyptus_dat, method="REML")
mod.7_tree_with_REML <- lme(log(h)~log(hdom)+log(dg_dbh)+
                                 (inv_dbh_age)+(inv_dbh), random = ~1|Arvore_key,
                               Eucalyptus_dat, method="REML")

###########################################################################
##############################  Test #######################################
############################################################################
########################## without vs. stands ##############################
############################################################################
logLik(mod.1_with_lm) 
logLik(mod.6_stand_with_REML) 
2*(644.3587-(612.833 )) 
1-0.5*pchisq(63.0514,0)-0.5*pchisq(63.0514,1)
############################################################################
##############################  Test #######################################
############################################################################
########################## stands vs. stands/tree ##########################
############################################################################
logLik(mod.6_stand_with_REML) 
logLik(mod.3_st.tree_with_REML) 
2*(645.42-(644.3587)) 
1-0.5*pchisq(1.9226,0)-0.5*pchisq(1.9226,1) 
############################################################################
##############################  Test #######################################
############################################################################
################### stands/tree vs. rf/stands/tree #########################
############################################################################
logLik(mod.3_st.tree_with_REML) 
logLik(mod.2_rf.st.tree_with_REML) 
2*(636.6843-(645.42)) 
1-0.5*pchisq(-17.4714,0)-0.5*pchisq(-17.4714,1) 
############################################################################
##############################  Test #######################################
############################################################################
########################## without vs. rf ##################################
############################################################################
logLik(mod.1_with_lm) 
logLik(mod.5_rf_with_REML) 
2*(637.8785-(612.833 )) 
1-0.5*pchisq(50.091,0)-0.5*pchisq(50.091,1) 
############################################################################
##############################  Test #######################################
############################################################################
########################## rf vs. rf/tree ##################################
############################################################################
logLik(mod.5_rf_with_REML) 
logLik(mod.4_rf.tree_with_REML) 
2*(638.90-(637.8785)) 
1-0.5*pchisq(2.043,0)-0.5*pchisq(2.043,1) 
############################################################################
##############################  Test #######################################
############################################################################
################### rf/tree vs. rf/stands/tree #############################
############################################################################
logLik(mod.4_rf.tree_with_REML) 
logLik(mod.2_rf.st.tree_with_REML) 
2*(636.6843 -(638.90)) 
1-0.5*pchisq(-4.4314,0)-0.5*pchisq(-4.4314,1) 
############################################################################
##############################  Test #######################################
############################################################################
########################## without vs. tree ################################
############################################################################
logLik(mod.1_with_lm) 
logLik(mod.7_tree_with_REML) 
2*(622.20-(612.833 )) 
1-0.5*pchisq(18.37,0)-0.5*pchisq(18.37,1) 

## AIC
AIC(mod.1_with_REML)
AIC(mod.2_rf.st.tree_with_REML)
AIC(mod.3_st.tree_with_REML)
AIC(mod.4_rf.tree_with_REML)
AIC(mod.5_rf_with_REML)
AIC(mod.6_stand_with_REML)
AIC(mod.7_tree_with_REML)
## cAIC
cAIC(mod.2_rf.st.tree_with_REML)
cAIC(mod.3_st.tree_with_REML)
cAIC(mod.4_rf.tree_with_REML)
cAIC(mod.5_rf_with_REML)
cAIC(mod.6_stand_with_REML)
cAIC(mod.7_tree_with_REML)


logLik(mod.1_with_REML)
logLik(mod.2_rf.st.tree_with_REML)
logLik(mod.3_st.tree_with_REML)
logLik(mod.4_rf.tree_with_REML)
logLik(mod.5_rf_with_REML)
logLik(mod.6_stand_with_REML)
logLik(mod.7_tree_with_REML)


VarCorr(mod.2_rf.st.tree_with_REML)
VarCorr(mod.3_st.tree_with_REML)
VarCorr(mod.4_rf.tree_with_REML)
VarCorr(mod.5_rf_with_REML)
VarCorr(mod.6_stand_with_REML)
VarCorr(mod.7_tree_with_REML)

## -2log(L)
## Local influence analysis 
# This code was taken from: http://www.ime.usp.br/~jmsinger/lmmdiagnostics.zip
#
# The theoretical procedure is based on: 
# Lesaffre, E. and G. Verbeke. 1998. Local influence in linear mixed models. 
# Biometrics, pages 570-582. doi: 10.2307/3109764.
#----
# Verbeke, G., and G. Molenberghs. 2000. Linear mixed models for longitudinal 
# data. Springer Science & Business Media, 568 p. 

#
# Initially checking whether the model readjusted via lme is the same as the 
# one via lmer.
# The fit model will be used in the modifiedlesafre_verbeke function process
#

fit <- lme(log(h)~log(hdom)+log(dg_dbh)+
             (inv_dbh_age)+(inv_dbh),
           random= list(activity = pdIdent(form = ~ 0 + as.factor(stand)), rf = ~1),
           control = lmeControl(opt = "optim"),
           data=Eucalyptus_dat) # this

fitLMER <- lmer(log(h)~log(hdom)+log(dg_dbh)+
                  (inv_dbh_age)+(inv_dbh) + (1|stand) + (1|rf),
                Eucalyptus_dat, REML = T)
#are equal?
all.equal(REMLcrit(fitLMER), c(-2*logLik(fit))) 
# fixed effects are equal? 
all.equal(fixef(fit), fixef(fitLMER))
# Ok, are equal!!!!
#
# Function: Modified Lesaffre-Verbeke index versus unit indices
#
modifiedlesafre_verbeke = function(dados, limit, plotid=NULL) {
  #
  # This function obtains the square root of a matrix                                                         ##
  #
  sqrt.matrix <- function(mat) {              
    mat <- as.matrix(mat)  
    singular_dec <- svd(mat)
    U <- singular_dec$u
    V <- singular_dec$v
    D <- diag(singular_dec$d)
    sqrtmatrix <- U %*% sqrt(D) %*% t(V)
  }
  #
  # This function extracts various objects of the function lme                                                ##
  #
  extract.lmeDesign2 <- function(m){
    start.level = 1
    data <- getData(m)
    grps <- nlme::getGroups(m)
    n <- length(grps)
    X <- list()
    grp.dims <- m$dims$ncol
    Zt <- model.matrix(m$modelStruct$reStruct, data)
    cov <- as.matrix(m$modelStruct$reStruct)
    i.col <- 1
    n.levels <- length(m$groups)
    Z <- matrix(0, n, 0)
    if (start.level <= n.levels) {
      for (i in 1:(n.levels - start.level + 1)) {
        if (length(levels(m$groups[[n.levels - i + 1]])) != 1)
        {
          X[[1]] <- model.matrix(~m$groups[[n.levels - i +
                                              1]] - 1, 
                                 contrasts.arg = c("contr.treatment",
                                                   "contr.treatment"))
        }
        else X[[1]] <- matrix(1, n, 1)
        X[[2]] <- as.matrix(Zt[, i.col:(i.col + grp.dims[i] -
                                          1)])
        i.col <- i.col + grp.dims[i]
        Z <- cbind(mgcv::tensor.prod.model.matrix(X),Z)
      }
      Vr <- matrix(0, ncol(Z), ncol(Z))
      start <- 1
      for (i in 1:(n.levels - start.level + 1)) {
        k <- n.levels - i + 1
        for (j in 1:m$dims$ngrps[i]) {
          stop <- start + ncol(cov[[k]]) - 1
          Vr[ncol(Z) + 1 - (stop:start),ncol(Z) + 1 - (stop:start)] <- cov[[k]]
          start <- stop + 1
        }
      }
    }
    X <- if (class(m$call$fixed) == "name" &&  !is.null(m$data$X)) {
      m$data$X
    } else   {
      model.matrix(formula(eval(m$call$fixed)),data)
    }
    y <- as.vector(matrix(m$residuals, ncol = NCOL(m$residuals))[,NCOL(m$residuals)] + 
                     matrix(m$fitted, ncol = NCOL(m$fitted))[,NCOL(m$fitted)])
    return(list(
      Vr = Vr,                                                                 
      X = X,
      Z = Z,
      sigmasq = m$sigma ^ 2,
      lambda = unique(diag(Vr)),
      y = y,
      k = n.levels
    )
    )
  }
  fit <- lme(log(h)~log(hdom)+log(dg_dbh)+
               (inv_dbh_age)+(inv_dbh),random = ~1|stand,
             dados, method="REML")
  data.fit <- extract.lmeDesign2(fit)
  data <-    getData(fit)
  y <- data.fit$y
  X <- data.fit$X
  N <- length(y)                                                               
  id <-  sort(as.numeric(getGroups(fit, level = 1)), index.return = TRUE)$x     
  subject <- as.numeric(unique(id))
  n <- length(as.numeric(names(table(id))))                                    
  vecni <- (table(id))                                                         
  p <- ncol(X)                                                                 
  n.levels <- length(fit$groups)                                              
  start.level <- 1
  Cgrps <- nlme::getGroups(fit, level = start.level)                           
  CCind <- levels((Cgrps))                                                     
  sigma2 <- fit$sigma^2
  obs <- numeric()
  
  for (i in 1:n)
  {
    obs <- append(obs,1:vecni[i])                                                
  }
  #
  # Construction of the Z and Gam matrices                                                                   ##
  #
  if (n.levels > 1) { 
    lZi <- list()
    lgi <- list()
    numrow <- numeric()
    
    mgroups <- fit$groups      
    for (n in 1:length(CCind)) {
      dgi <- data.frame(as.matrix(mgroups[mgroups == CCind[n], ]))
      nrowzi <- dim(dgi)[1]
      ncolzi <- 0
      # Number of repetitions of the variance components to construct the Gi matrix   
      girep <- as.numeric(length(levels(dgi[,1])))
      for (k in 2:n.levels) {
        girep <- c(girep,as.numeric(length(levels(dgi[,k]))))
      }
      # Number of columns of the Zi matrix   
      for (k in 1:n.levels) {
        ncolzi <- ncolzi + as.numeric(length(levels(dgi[,k])))
      }
      # Numbers of one's by columns of the Zi matrix
      auxi <- as.vector(table(dgi[,1]))
      for (i in 2:n.levels) {
        auxi <- c(auxi,as.vector(table(dgi[,i])))
      }
      # Matrix Zi
      l <- 1
      Zi <- matrix(0,nrowzi,ncolzi)
      # Inserting elements in Zi
      for (j in 1:ncolzi) {
        Zi[l:(l + auxi[j] - 1),j] <- rep(1,auxi[j]) 
        l <- l + auxi[j]
        if (l == (nrowzi + 1)) l <- 1
      }
      
      lZi[[n]] <- Zi
      numrow[n] <- dim(Zi)[1]
      comp.var <- as.matrix(fit1$modelStruct$reStruct)
      auxg <- rep(as.numeric(comp.var[1])*sigma2,girep[1])
      for (i in 2:length(girep)) {
        auxg <- c(auxg,rep(as.numeric(comp.var[i])*sigma2,girep[i]))
      }
      lgi[[n]] <- diag(auxg)
    }
    q <- dim(lgi[[1]])[1]                     # Dimensions of Gi matrices
    for (h in 2:length(CCind)) {
      q <- c(q,dim(lgi[[h]])[1])
    }
    Z <- lZi[[1]]
    for (k in 2:length(CCind)) {
      Z <- bdiag(Z,(lZi[[k]]))
    }
    Z <- as.matrix(Z)
    nrowZi <- lZi[[1]]                        # Dmensions of Zi matrices
    for (h in 2:length(CCind)) {
      nrowZi <- c(nrowZi,dim(lZi[[h]])[1])
    }
    
    Gam <- lgi[[1]]
    for (k in 2:length(CCind)) {
      Gam <- bdiag(Gam,(lgi[[k]]))
    }
    Gam <- as.matrix(Gam)
  }else{
    mataux <- model.matrix(fit$modelStruct$reStruct,data)
    mataux <- as.data.frame(cbind(mataux,id))
    lZi <- list()
    lgi <- list()
    
    for (i in (as.numeric(unique(id)))) { 
      lZi[[i]] <- as.matrix((subset(split(mataux,id == i,
                                          drop = T)$`TRUE`,select = -id)))          
      lgi[[i]] <- getVarCov(fit,type = "random.effects")
    }
    Z <- as.matrix(bdiag(lZi))
    g <- getVarCov(fit,type = "random.effects")
    q <- dim(g)[1]                                                           # Total number of random effects
    Gam <- as.matrix(kronecker(diag(length(as.numeric(unique(id)))),g))
  }
  #
  ## Estimate of the covariance matrix of conditional errors (homoskedastic conditional independence model)   ##
  #
  if (n.levels > 1) {   
    if (!inherits(fit, "lme")) 
      stop("object does not appear to be of class lme")
    grps <- nlme::getGroups(fit)
    n <- length(grps)                                                                      # Number of observations
    n.levels <- length(fit$groups)                                                         # Number of levels
    if (is.null(fit$modelStruct$corStruct)) 
      n.corlevels <- 0
    else n.corlevels <- length(all.vars(nlme::getGroupsFormula(fit$modelStruct$corStruct))) # Levels of the repeated measures
    if (n.levels < n.corlevels) {
      getGroupsFormula(fit$modelStruct$corStruct)
      vnames <- all.vars(nlme::getGroupsFormula(fit$modelStruct$corStruct))
      lab <- paste(eval(parse(text = vnames[1]), envir = fit$data))
      if (length(vnames) > 1) 
        for (i in 2:length(vnames)) {
          lab <- paste(lab, "/", eval(parse(text = vnames[i]), 
                                      envir = fit$data), sep = "")
        }
      grps <- factor(lab)
    }
    if (n.levels >= start.level || n.corlevels >= start.level) {
      if (n.levels >= start.level) 
        Cgrps <- nlme::getGroups(fit, level = start.level)                           # Level = 1
      else Cgrps <- grps
      Cind <- sort(as.numeric(Cgrps), index.return = TRUE)$ix                        # Indices of the observations
      rCind <- 1:n 
      rCind[Cind] <- 1:n
      Clevel <- levels(Cgrps)                                                        # Levels of the first nesting level
      n.cg <- length(Clevel)                                                         
      size.cg <- array(0, n.cg)
      for (i in 1:n.cg) size.cg[i] <- sum(Cgrps == Clevel[i])                        # Number of the observations by subject
    }
    else {
      n.cg <- 1
      Cind <- 1:n
    }
    if (is.null(fit$modelStruct$varStruct)) 
      w <- rep(fit$sigma, n)
    else {
      w <- 1/nlme::varWeights(fit$modelStruct$varStruct)
      group.name <- names(fit$groups)
      order.txt <- paste("ind<-order(data[[\"", group.name[1], 
                         "\"]]", sep = "")
      if (length(fit$groups) > 1) 
        for (i in 2:length(fit$groups)) order.txt <- paste(order.txt, 
                                                           ",data[[\"", group.name[i], "\"]]", sep = "")
      order.txt <- paste(order.txt, ")")
      eval(parse(text = order.txt))
      w[ind] <- w
      w <- w * fit$sigma
    }
    w <- w[Cind]
    if (is.null(fit$modelStruct$corStruct)) 
      lR <- array(1, n)
    else {
      c.m <- nlme::corMatrix(fit$modelStruct$corStruct)
      if (!is.list(c.m)) {
        lR <- c.m
        lR <- lR[Cind, ]
        lR <- lR[, Cind]
      }
      else {
        lR <- list()
        ind <- list()
        for (i in 1:n.cg) {
          lR[[i]] <- matrix(0, size.cg[i], size.cg[i])
          ind[[i]] <- 1:size.cg[i]
        }
        Roff <- cumsum(c(1, size.cg))
        gr.name <- names(c.m)
        n.g <- length(c.m)
        j0 <- rep(1, n.cg)
        ii <- 1:n
        for (i in 1:n.g) {
          Clev <- unique(Cgrps[grps == gr.name[i]])
          if (length(Clev) > 1) 
            stop("inner groupings not nested in outer!!")
          k <- (1:n.cg)[Clevel == Clev]
          j1 <- j0[k] + nrow(c.m[[i]]) - 1
          lR[[k]][j0[k]:j1, j0[k]:j1] <- c.m[[i]]
          ind1 <- ii[grps == gr.name[i]]
          ind2 <- rCind[ind1]
          ind[[k]][j0[k]:j1] <- ind2 - Roff[k] + 1
          j0[k] <- j1 + 1
        }
        for (k in 1:n.cg) {
          lR[[k]][ind[[k]], ] <- lR[[k]]
          lR[[k]][, ind[[k]]] <- lR[[k]]
        }
      }
    }
    if (is.list(lR)) {
      for (i in 1:n.cg) {
        wi <- w[Roff[i]:(Roff[i] + size.cg[i] - 1)]
        lR[[i]] <- as.vector(wi) * t(as.vector(wi) * lR[[i]]) # Matrix lR
      }
    }
    else if (is.matrix(lR)) {
      lR <- as.vector(w) * t(as.vector(w) * lR)
    }
    else {
      lR <- w^2 * lR
    }
    if (is.list(lR)) {
      R <- lR[[1]]
      for (k in 2:n.cg) {
        R <- bdiag(R,lR[[k]])
      }
      R <- as.matrix(R)
    }
    else{
      R <- diag(lR)
    }
  }else{
    R <- getVarCov(fit,type = "conditional",individual = 1)[[1]]
    for (i in 2:length(as.numeric(unique(id)))) {
      R <- as.matrix(bdiag(R,getVarCov(fit,
                                       type = "conditional",individual = i)[[1]] ) )
    }
  }
  #
  ## Construction of covariance matrix of Y           
  #
  V <- (Z %*% Gam %*% t(Z)) + R
  iV <- solve(V)                                                
  varbeta <- solve((t(X) %*% iV %*% X))
  #
  ## EBLUE                                                                                        ##
  #
  eblue <- as.vector(fixef(fit))
  #
  ## Residual analysis                                                                                      ##
  #  
  predm <- X %*% eblue                       # Predicted values for expected response
  resm <- (y - predm)                        # Marginal residuals
  #
  ## Variance of marginal residuals                                                                         ##
  # 
  var.resm <- V - X %*% solve(t(X) %*% iV %*% X) %*% t(X) 
  #
  ## Modified Lesaffre-Verbeke index                                                                        ##
  #
  lesverb <- rep(0,length(CCind)) 
  auxni <- as.vector(vecni)
  for (t in 1:length(CCind)) { 
    li <- sum(vecni[1:t - 1]) + 1
    ls <- sum(vecni[1:t])
    if (vecni[t] == 1) {
      
      auxr2 <- solve(sqrt(var.resm[li:ls,li:ls]))
      Ri <- (auxr2) %*% resm[li:ls]
      auxt <- diag(vecni[t]) - Ri %*% t(Ri)
      lesverb[t] <- sqrt(sum(diag(auxt %*% t(auxt))))
    }
    else
    {  
      
      auxr2 <- solve(sqrt.matrix(var.resm[li:ls,li:ls]))
      Ri <- auxr2 %*% resm[li:ls]
      auxt <- diag(vecni[t]) - Ri %*% t(Ri)
      lesverb[t] <- sqrt(sum(diag(auxt %*% t(auxt))))
    } 
  }
  
  lesverb <- lesverb/((as.numeric((table(id)))))
  LSlesverb <- as.numeric(quantile(lesverb,prob = 0.75) + 1.5*(quantile(lesverb, prob = 0.75) - quantile(lesverb,prob = 0.25)))
  LS2lesverb <- 2*mean(lesverb)
  # 
  ## This function constructs the diagnostic plots                                                         ##
  # 
  
  plotg = function(plotid){
    cat("\n To select the graphic use plotid \n
1 - Modified Lesaffre-Verbeke index versus unit indices
\n")
    cat("\n Graph plotting", plotid)
    
    if (plotid == 1)
    {
      par(mfrow = c(1,1),mar = c(11, 5, 1, 2))
      plot(lesverb,ylab = expression(paste(C*""[i])),
           xlab = "Individual number", cex = 1.2, cex.lab = 1.8, cex.axis = 1.3, 
           pch = 20, ylim = c(0,2*max(abs(range(lesverb)))))
      abline(h = LSlesverb,lty = 2)                            # 3rd quartile + 1.5*interquartile interval
      #     abline(h = LS2lesverb,lty = 3)                           # 2*mean(lesverbp)
      index = which(lesverb > LSlesverb)
      #      index1<-subject[index]
      if (length(index) > 0)
      {
        text(index, lesverb[index], index, adj = c(1,-.5), cex = 1.0, font = 2)
      }
    }
    
  }
  if (is.null(plotid)) {
    cat("\n To choose plot, select plotid \n
1 - Modified Lesaffre-Verbeke index versus unit indices
    \n")
    return(1);  
  }
  
  #
  # Generation of diagnostic plots                                                                   ##
  #
  for (g in plotid) {
    plotg(g)
    cat("\n Press ENTER to continue...")
    readline()
  }
  
  useful.results <- list(
    lesaffreverbeke.measure = cbind(Subject = as.numeric(unique(id)),LV.m = lesverb))
}
#
# plot: Modified Lesaffre-Verbeke index versus unit indices
#
x11()
modifiedlesafre_verbeke(dados = Eucalyptus_dat, 
                        limit = c(-5,5), plotid = c(1))
#
# The Modified Lesaffre-Verbeke index (C_i) shows that the outlier observation is 
# contained in stand 35 (Q4A). 
# After a detailed analysis, it was identified that this is observation 509, which
# refers to a tree 5.4 meters, with a DBH of 8.37 cm and classified as bifurcated, 
# i.e. it is understood that there was no evolution in height due to its 
# characteristic. 

# Removing the outlier observation
Eucalyptus_dat_New = Eucalyptus_dat[-509,]

# ## Variance Component Tests 
# 
# ## H0: no random effects 
# ## H1: sigma^2_stands 
# ## i.e. sigma^2_stands = 0 vs. sigma^2_stands > 0 
# mod.1_without <- lm(log(h)~log(hdom)+log(dg_dbh)+
#                             (inv_dbh_age)+(inv_dbh),
#                           Eucalyptus_dat_New)
# mod.2_st_without <- lmer(log(h)~log(hdom)+log(dg_dbh)+
#                               (inv_dbh_age)+(inv_dbh) + (1|stand),
#                             Eucalyptus_dat_New, REML = F)
# ranova(mod.2_st_without) 
# 2*(765.65-702.45) #= 126.4
# 1-0.5*pchisq(126.4,0)-0.5*pchisq(126.4,1) # p-value: 0
# ## varTestnlme package
# t1 = varCompTest(mod.2_st_without, mod.1_without);summary(t1)
# 
# ## H0: sigma^2_stands 
# ## H1: sigma^2_stands & sigma^2_forestregion 
# ##  i.e. sigma^2_forestregion = 0 vs. sigma^2_forestregion > 0 
# mod.2_st_without <- lmer(log(h)~log(hdom)+log(dg_dbh)+
#                               (inv_dbh_age)+(inv_dbh) + (1|stand),
#                             Eucalyptus_dat_New, REML = F)
# mod.4_both_without <- lmer(log(h)~log(hdom)+log(dg_dbh)+
#                             (inv_dbh_age)+(inv_dbh) + (1|stand) + (1|rf),
#                           Eucalyptus_dat_New, REML = F)
# ranova(mod.4_both_without) 
# 2*(766.73-765.65) #= 2.16
# 1-0.5*pchisq(2.16,0)-0.5*pchisq(2.16,1) # p-value: 0.07082235
# ## varTestnlme package
# t2 = varCompTest(mod.4_both_without, mod.2_st_without);summary(t2)
# # Variance components testing in mixed effects models
# # Error in extractStruct.merMod(m1, m0, randm0) : 
# # the package does not currently support more than one level of random effects
# 
# ## H0: no random effects 
# ## H1: sigma^2_rf 
# ## i.e. sigma^2_rf = 0 vs. sigma^2_rfs > 0 
# mod.1_without <- lm(log(h)~log(hdom)+log(dg_dbh)+
#                             (inv_dbh_age)+(inv_dbh),
#                           Eucalyptus_dat_New)
# mod.3_fr_without <- lmer(log(h)~log(hdom)+log(dg_dbh)+
#                           (inv_dbh_age)+(inv_dbh) + (1|rf),
#                         Eucalyptus_dat_New, REML = F)
# ranova(mod.3_fr_without) 
# 2*(745.73-702.45)
# 1-0.5*pchisq(86.56,0)-0.5*pchisq(86.56,1) # p-value: 0
# 
# ## H0: sigma^2_forestregion
# ## H1: sigma^2_stands & sigma^2_forestregion 
# ##  i.e. sigma^2_forestregion = 0 vs. sigma^2_forestregion > 0 
# mod.3_fr_without <- lmer(log(h)~log(hdom)+log(dg_dbh)+
#                           (inv_dbh_age)+(inv_dbh) + (1|rf),
#                         Eucalyptus_dat_New, REML = F)
# mod.4_both_without <- lmer(log(h)~log(hdom)+log(dg_dbh)+
#                             (inv_dbh_age)+(inv_dbh) + (1|stand) + (1|rf),
#                           Eucalyptus_dat_New, REML = F)
# ranova(mod.4_both_without) 
# 2*(766.73-745.73) #= 42
# 1-0.5*pchisq(42,0)-0.5*pchisq(42,1) # p-value: 0
# 
# ## H0: no random effects
# ## H1: sigma^2_stands & sigma^2_forestregion 
# ## i.e. sigma^2_forestregion & sigma^2_forestregion = 0 vs. > 0
# mod.1_without <- lm(log(h)~log(hdom)+log(dg_dbh)+
#                       (inv_dbh_age)+(inv_dbh),
#                     Eucalyptus_dat_New)
# mod.4_both_without <- lmer(log(h)~log(hdom)+log(dg_dbh)+
#                              (inv_dbh_age)+(inv_dbh) + (1|stand) + (1|rf),
#                            Eucalyptus_dat_New, REML = F)
# 2*(766.73-702.45) #= 128.56
# 1-0.25*pchisq(128.56,0)-0.5*pchisq(128.56,1)-0.25*pchisq(128.56,2) # p-value: 0
# 
# # Wald tests for fixed effects. 
# anova(mod.2_st_without, ddf = "Kenward-Roger") 

## fit with REML = T

mod.1_without_lm <- lm(log(h)~log(hdom)+log(dg_dbh)+
                      (inv_dbh_age)+(inv_dbh),
                    Eucalyptus_dat_New)

mod.2_rf.st.tree_without_REML <- lme(log(h)~log(hdom)+log(dg_dbh)+
                                  (inv_dbh_age)+(inv_dbh), random = ~1|rf/stand/Arvore_key,
                                Eucalyptus_dat_New, method="REML")
mod.3_st.tree_without_REML <- lme(log(h)~log(hdom)+log(dg_dbh)+
                                  (inv_dbh_age)+(inv_dbh), random = ~1|stand/Arvore_key,
                                Eucalyptus_dat_New, method="REML")
mod.4_rf.tree_without_REML <- lme(log(h)~log(hdom)+log(dg_dbh)+
                                  (inv_dbh_age)+(inv_dbh), random = ~1|rf/Arvore_key,
                                Eucalyptus_dat_New, method="REML")
mod.5_rf_without_REML <- lme(log(h)~log(hdom)+log(dg_dbh)+
                                 (inv_dbh_age)+(inv_dbh), random = ~1|rf,
                               Eucalyptus_dat_New, method="REML")
mod.6_stand_without_REML <- lme(log(h)~log(hdom)+log(dg_dbh)+
                                 (inv_dbh_age)+(inv_dbh), random = ~1|stand,
                               Eucalyptus_dat_New, method="REML")
mod.7_tree_without_REML <- lme(log(h)~log(hdom)+log(dg_dbh)+
                                 (inv_dbh_age)+(inv_dbh), random = ~1|Arvore_key,
                               Eucalyptus_dat_New, method="REML")
###########################################################################
##############################  Test #######################################
############################################################################
########################## without vs. stands ##############################
############################################################################
logLik(mod.1_without_lm) 
logLik(mod.6_stand_without_REML) 
2*(756.4761-(702.4498)) 
1-0.5*pchisq(108.0526,0)-0.5*pchisq(108.0526,1)
############################################################################
##############################  Test #######################################
############################################################################
########################## stands vs. stands/tree ##########################
############################################################################
logLik(mod.6_stand_without_REML) 
logLik(mod.3_st.tree_without_REML) 
2*(757.34-(756.4761)) 
1-0.5*pchisq(1.7278,0)-0.5*pchisq(1.7278,1) 
############################################################################
##############################  Test #######################################
############################################################################
################### stands/tree vs. rf/stands/tree #########################
############################################################################
logLik(mod.3_st.tree_without_REML) 
logLik(mod.2_rf.st.tree_without_REML) 
2*(738.1034-(757.34)) 
1-0.5*pchisq(-38.4732,0)-0.5*pchisq(-38.4732,1) 
############################################################################
##############################  Test #######################################
############################################################################
########################## without vs. rf ##################################
############################################################################
logLik(mod.1_without_lm) 
logLik(mod.5_rf_without_REML) 
2*(736.0902-(702.4498)) 
1-0.5*pchisq(67.2808,0)-0.5*pchisq(67.2808,1) 
############################################################################
##############################  Test #######################################
############################################################################
########################## rf vs. rf/tree ##################################
############################################################################
logLik(mod.5_rf_without_REML) 
logLik(mod.4_rf.tree_without_REML) 
2*(737.03-(736.0902)) 
1-0.5*pchisq(1.8796,0)-0.5*pchisq(1.8796,1) 
############################################################################
##############################  Test #######################################
############################################################################
################### rf/tree vs. rf/stands/tree #############################
############################################################################
logLik(mod.4_rf.tree_without_REML) 
logLik(mod.2_rf.st.tree_without_REML) 
2*(738.1034-(737.03)) 
1-0.5*pchisq(2.1468,0)-0.5*pchisq(2.1468,1) 
############################################################################
##############################  Test #######################################
############################################################################
########################## without vs. tree ################################
############################################################################
logLik(mod.1_without_lm) 
logLik(mod.7_tree_without_REML) 
2*(711.14-(702.4498)) 
1-0.5*pchisq(17.3804,0)-0.5*pchisq(17.3804,1) 




## AIC
AIC(mod.1_without_REML)
AIC(mod.2_rf.st.tree_without_REML)
AIC(mod.3_st.tree_without_REML)
AIC(mod.4_rf.tree_without_REML)
AIC(mod.5_rf_without_REML)
AIC(mod.6_stand_without_REML)
AIC(mod.7_tree_without_REML)
## cAIC
cAIC(mod.2_rf.st.tree_without_REML)
cAIC(mod.3_st.tree_without_REML)
cAIC(mod.4_rf.tree_without_REML)
cAIC(mod.5_rf_without_REML)
cAIC(mod.6_stand_without_REML)
cAIC(mod.7_tree_without_REML)
## -2log(L)
logLik(mod.1_without_REML)
logLik(mod.2_rf.st.tree_without_REML)
logLik(mod.3_st.tree_without_REML)
logLik(mod.4_rf.tree_without_REML)
logLik(mod.5_rf_without_REML)
logLik(mod.6_stand_without_REML)
logLik(mod.7_tree_without_REML)

VarCorr(mod.2_rf.st.tree_without_REML)
VarCorr(mod.3_st.tree_without_REML)
VarCorr(mod.4_rf.tree_without_REML)
VarCorr(mod.5_rf_without_REML)
VarCorr(mod.6_stand_without_REML)
VarCorr(mod.7_tree_without_REML)
# Least confounded residuals
source("functions_eucal_mixedmodels.R")
#----------------------------------------------------
mod.2_rf.st.tree_with_REML <- lme(log(h)~log(hdom)+log(dg_dbh)+
                                       (inv_dbh_age)+(inv_dbh), random = ~1|rf/stand/Arvore_key,
                                     Eucalyptus_dat, method="REML")
#----------------------------------------------------
envelope_mod.2_rf.st.tree_with_REML = res_lcr(fit = mod.2_rf.st.tree_with_REML, 
                            graph = "envelope")
hist_mod.2_rf.st.tree_with_REML = res_lcr(fit = mod.2_rf.st.tree_with_REML, 
                           graph = "hist")
mahalanobis_mod.2_rf.st.tree_with_REML = res_lcr(fit = mod.2_rf.st.tree_with_REML, 
                            graph = "mahalanobis distance")
dispersion_mod.2_rf.st.tree_with_REML = res_lcr(fit = mod.2_rf.st.tree_with_REML, 
                            graph = "dispersion")
#----------------------------------------------------
mod.2_st_nlme_without <- lme(log(h)~log(hdom)+log(dg_dbh)+
                               (inv_dbh_age)+(inv_dbh), random = ~1|rf/stand/Arvore_key,
                             Eucalyptus_dat_New, method="REML")
#----------------------------------------------------
envelope_mod.2_st_nlme_without = res_lcr(fit = mod.2_st_nlme_without, 
                            graph = "envelope")
hist_mod.2_st_nlme_without = res_lcr(fit = mod.2_st_nlme_without, 
                            graph = "hist")
mahalanobis_mod.2_st_nlme_without = res_lcr(fit = mod.2_st_nlme_without, 
                            graph = "mahalanobis distance")
dispersion_mod.2_st_nlme_without = res_lcr(fit = mod.2_st_nlme_without, 
                            graph = "dispersion")

grid.arrange(envelope_mod.2_st_nlme, hist_mod.2_st_nlme, 
             mahalanobis_mod.2_st_nlme, dispersion_mod.2_st_nlme,
             envelope_mod.2_st_nlme_without, hist_mod.2_st_nlme_without, 
             mahalanobis_mod.2_st_nlme_without, dispersion_mod.2_st_nlme_without, 
             ncol = 4)
#----------------------------------------------------
mod.3_fr_nlme_with <- lme(log(h)~log(hdom)+log(dg_dbh)+
                            (inv_dbh_age)+(inv_dbh), random = ~1|rf,
                          Eucalyptus_dat, method="REML")
#----------------------------------------------------
envelope_mod.3_fr_nlme = res_lcr(fit = mod.3_fr_nlme_with, 
                                 graph = "envelope")
hist_mod.3_fr_nlme = res_lcr(fit = mod.3_fr_nlme_with, 
                             graph = "hist")
mahalanobis_mod.3_fr_nlme = res_lcr(fit = mod.3_fr_nlme_with, 
                                    graph = "mahalanobis distance")
dispersion_mod.3_fr_nlme = res_lcr(fit = mod.3_fr_nlme_with, 
                                   graph = "dispersion")
#----------------------------------------------------
mod.3_fr_nlme_without <- lme(log(h)~log(hdom)+log(dg_dbh)+
                               (inv_dbh_age)+(inv_dbh), random = ~1|rf,
                             Eucalyptus_dat_New, method="REML")
#----------------------------------------------------
envelope_mod.3_fr_nlme_without = res_lcr(fit = mod.3_fr_nlme_without, 
                                         graph = "envelope")
hist_mod.3_fr_nlme_without = res_lcr(fit = mod.3_fr_nlme_without, 
                                     graph = "hist")
mahalanobis_mod.3_fr_nlme_without = res_lcr(fit = mod.3_fr_nlme_without, 
                                            graph = "mahalanobis distance")
dispersion_mod.3_fr_nlme_without = res_lcr(fit = mod.3_fr_nlme_without, 
                                           graph = "dispersion")

grid.arrange(envelope_mod.3_fr_nlme, hist_mod.3_fr_nlme, 
             mahalanobis_mod.3_fr_nlme, dispersion_mod.3_fr_nlme,
             envelope_mod.3_fr_nlme_without, hist_mod.3_fr_nlme_without, 
             mahalanobis_mod.3_fr_nlme_without, dispersion_mod.3_fr_nlme_without, 
             ncol = 4)

#----------------------------------------------------
mod.3_tree_nlme_with <- lme(log(h)~log(hdom)+log(dg_dbh)+
                            (inv_dbh_age)+(inv_dbh), random = ~1|Arvore_key,
                          Eucalyptus_dat, method="REML")
#----------------------------------------------------
envelope_mod.3_tree_nlme = res_lcr(fit = mod.3_tree_nlme_with, 
                                 graph = "envelope")
hist_mod.3_tree_nlme = res_lcr(fit = mod.3_tree_nlme_with, 
                             graph = "hist")
mahalanobis_mod.3_tree_nlme = res_lcr(fit = mod.3_tree_nlme_with, 
                                    graph = "mahalanobis distance")
dispersion_mod.3_tree_nlme = res_lcr(fit = mod.3_tree_nlme_with, 
                                   graph = "dispersion")
#----------------------------------------------------
mod.3_tree_nlme_without <- lme(log(h)~log(hdom)+log(dg_dbh)+
                               (inv_dbh_age)+(inv_dbh), random = ~1|Arvore_key,
                             Eucalyptus_dat_New, method="REML")
#----------------------------------------------------
envelope_mod.3_tree_nlme_without = res_lcr(fit = mod.3_tree_nlme_without, 
                                         graph = "envelope")
hist_mod.3_tree_nlme_without = res_lcr(fit = mod.3_tree_nlme_without, 
                                     graph = "hist")
mahalanobis_mod.3_tree_nlme_without = res_lcr(fit = mod.3_tree_nlme_without, 
                                            graph = "mahalanobis distance")
dispersion_mod.3_tree_nlme_without = res_lcr(fit = mod.3_tree_nlme_without, 
                                           graph = "dispersion")
#------------------- relative changes on the parameters ------------------- 
# Fixed Effects 
# Extraction of fixed effects parameter estimates
without_outlier_model <- lmer(log(h)~log(hdom)+log(dg_dbh)+
                                (inv_dbh_age)+(inv_dbh) + (1|stand) + (1|rf),
                              Eucalyptus_dat_New, REML = T)
with_outlier_model <- lmer(log(h)~log(hdom)+log(dg_dbh)+
                             (inv_dbh_age)+(inv_dbh) + (1|stand) + (1|rf),
                           Eucalyptus_dat, REML = T)
theta_L <- fixed.effects(without_outlier_model); theta_L
theta <- fixed.effects(with_outlier_model); theta
# Function to calculate the percentage difference
calcular_diferenca_percentual <- function(theta_L, theta) {
  diff_percentual <- ((theta_L - theta) / theta) * 100
  return(diff_percentual)
}
# Call the function to calculate the percentage difference for each parameter
(diferenca_percentual <- calcular_diferenca_percentual(theta_L, theta))

# Random Effects 
# Extract the random effects parameters from the models 
# without and with the outlier obs

comp_var_theta_L_t <- VarCorr(without_outlier_model)[[1]][1]; comp_var_theta_L_t
comp_var_theta_L_rf <- VarCorr(without_outlier_model)[[2]][1]; comp_var_theta_L_rf
comp_var_theta_L_error <- (summary(without_outlier_model)$sigma)^2; comp_var_theta_L_error 

comp_var_theta_t <- VarCorr(with_outlier_model)[[1]][1]; comp_var_theta_t
comp_var_theta_rf <- VarCorr(with_outlier_model)[[2]][1]; comp_var_theta_rf
comp_var_theta_error <- (summary(with_outlier_model)$sigma)^2; comp_var_theta_error 

comp_var_theta_L <- c(comp_var_theta_L_t, comp_var_theta_L_rf, comp_var_theta_L_error)
comp_var_theta <- c(comp_var_theta_t, comp_var_theta_rf, comp_var_theta_error)
# Initialize a vector to store percentage differences
diferencas_percentuais_comp_var <- numeric(length(comp_var_theta_L))
# Loop to calculate the percentage difference for each variance component
for (i in seq_along(comp_var_theta_L)) {
  diff_percentual <- ((comp_var_theta_L[i] - comp_var_theta[i]) / comp_var_theta[i]) * 100
  diferencas_percentuais_comp_var[i] <- diff_percentual
}
# Display the percentage differences of variance components
print(diferencas_percentuais_comp_var)
#-------------------------------------------------------------------------- 