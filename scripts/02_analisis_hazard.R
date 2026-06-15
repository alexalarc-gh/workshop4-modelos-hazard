# ============================================================
# Workshop 4 - Modelos Hazard
# Analisis de supervivencia aplicado a confiabilidad operacional
# Base sintetica de estaciones volcanicas en Chile
# - carga y prepara la base
# - calcula descriptivos
# - ajusta Kaplan-Meier y Log-Rank
# - ajusta modelo Cox
# - valida supuesto de riesgos proporcionales
# - exporta tablas y figuras limpias para presentacion RMarkdown
# ============================================================

# ------------------------------------------------------------
# 1. Librerias
# ------------------------------------------------------------

library(tidyverse)
library(readr)
library(janitor)
library(here)
library(survival)
library(survminer)
library(broom)
library(gt)
library(scales)

# ------------------------------------------------------------
# 2. Directorios de salida
# ------------------------------------------------------------

dir.create(here("outputs"), showWarnings = FALSE)
dir.create(here("outputs", "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("outputs", "tables"), recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 3. Carga de datos
# ------------------------------------------------------------

hazard_data_raw <- read_csv(
  here("data", "synthetic_hazard_stations_chile.csv"),
  show_col_types = FALSE
)

# ------------------------------------------------------------
# 4. Preparacion de variables
# ------------------------------------------------------------

hazard_data <- hazard_data_raw %>%
  clean_names() %>%
  mutate(
    event = as.integer(event),
    time_days = as.numeric(time_days),
    age_years = as.numeric(age_years),
    battery_margin_days = as.numeric(battery_margin_days),
    winter_operation = as.integer(winter_operation),

    energy_system = factor(
      energy_system,
      levels = c("robust", "standard", "limited"),
      labels = c("Robusto", "Estándar", "Limitado")
    ),

    snow_ice_exposure = factor(
      snow_ice_exposure,
      levels = c("low", "medium", "high", "very_high"),
      labels = c("Baja", "Media", "Alta", "Muy alta")
    ),

    comm_type = factor(
      comm_type,
      levels = c("satellite", "radio", "cellular"),
      labels = c("Satelital", "Radio", "Celular")
    ),

    access_difficulty = factor(
      access_difficulty,
      levels = c("low", "medium", "high"),
      labels = c("Baja", "Media", "Alta")
    ),

    maintenance_history = factor(
      maintenance_history,
      levels = c("preventive_ok", "preventive_delayed", "corrective_only"),
      labels = c("Preventiva vigente", "Preventiva retrasada", "Solo correctiva")
    ),

    winter_operation = factor(
      winter_operation,
      levels = c(0, 1),
      labels = c("No", "Sí")
    )
  )

# Chequeo minimo de consistencia
stopifnot(all(hazard_data$event %in% c(0, 1)))
stopifnot(all(hazard_data$time_days > 0))

# ------------------------------------------------------------
# 5. Resumen descriptivo general
# ------------------------------------------------------------

summary_general <- tibble(
  indicador = c(
    "Estaciones sintéticas",
    "Volcanes representados",
    "Eventos/fallas observadas",
    "Observaciones censuradas",
    "Proporción de eventos",
    "Tiempo mínimo observado (días)",
    "Tiempo mediano observado (días)",
    "Tiempo máximo observado (días)"
  ),
  valor = c(
    nrow(hazard_data),
    n_distinct(hazard_data$volcano_name),
    sum(hazard_data$event == 1),
    sum(hazard_data$event == 0),
    round(mean(hazard_data$event == 1), 3),
    min(hazard_data$time_days),
    median(hazard_data$time_days),
    max(hazard_data$time_days)
  )
)

write_csv(summary_general, here("outputs", "tables", "01_summary_general.csv"))

# Eventos por sistema energetico
events_by_energy <- hazard_data %>%
  count(energy_system, event) %>%
  group_by(energy_system) %>%
  mutate(prop = round(n / sum(n), 3)) %>%
  ungroup()

write_csv(events_by_energy, here("outputs", "tables", "02_events_by_energy.csv"))

# Eventos por exposicion a nieve/hielo
events_by_snow <- hazard_data %>%
  count(snow_ice_exposure, event) %>%
  group_by(snow_ice_exposure) %>%
  mutate(prop = round(n / sum(n), 3)) %>%
  ungroup()

write_csv(events_by_snow, here("outputs", "tables", "03_events_by_snow.csv"))

# Tabla compacta para describir datos en la presentacion
data_overview_for_slide <- tibble(
  item = c("Estaciones", "Volcanes", "Fallas", "Censuras", "Horizonte máximo", "Tasa de falla observada"),
  value = c(
    nrow(hazard_data),
    n_distinct(hazard_data$volcano_name),
    sum(hazard_data$event == 1),
    sum(hazard_data$event == 0),
    paste0(max(hazard_data$time_days), " días"),
    percent(mean(hazard_data$event == 1), accuracy = 0.1)
  )
)

write_csv(data_overview_for_slide, here("outputs", "tables", "00_data_overview_for_slide.csv"))

# ------------------------------------------------------------
# 6. Objeto de supervivencia
# ------------------------------------------------------------

surv_obj <- Surv(
  time = hazard_data$time_days,
  event = hazard_data$event
)

# ------------------------------------------------------------
# 7. Kaplan-Meier global, figura limpia para presentacion
# ------------------------------------------------------------

km_global <- survfit(surv_obj ~ 1, data = hazard_data)

p_km_global <- ggsurvplot(
  km_global,
  data = hazard_data,
  conf.int = TRUE,
  risk.table = FALSE,
  censor = TRUE,
  legend = "none",
  xlab = "Tiempo de observación (días)",
  ylab = "Probabilidad de supervivencia operacional",
  title = "Curva Kaplan-Meier global",
  subtitle = "Base sintética de estaciones autónomas de monitoreo volcánico",
  ggtheme = theme_minimal(base_size = 13)
)

p_km_global$plot <- p_km_global$plot +
  coord_cartesian(xlim = c(0, 760), ylim = c(0, 1.02)) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 12),
    panel.grid.minor = element_blank()
  )

ggsave(
  filename = here("outputs", "figures", "01_kaplan_meier_global.png"),
  plot = p_km_global$plot,
  width = 10,
  height = 5.6,
  dpi = 300
)

# ------------------------------------------------------------
# 8. Kaplan-Meier por sistema energetico, figura limpia
# ------------------------------------------------------------

km_energy <- survfit(surv_obj ~ energy_system, data = hazard_data)

p_km_energy <- ggsurvplot(
  km_energy,
  data = hazard_data,
  conf.int = TRUE,
  risk.table = FALSE,
  censor = TRUE,
  pval = TRUE,
  pval.coord = c(25, 0.18),
  legend.title = "Sistema energético",
  legend.labs = c("Robusto", "Estándar", "Limitado"),
  xlab = "Tiempo de observación (días)",
  ylab = "Probabilidad de supervivencia operacional",
  title = "Supervivencia operacional según sistema energético",
  subtitle = "Comparación Kaplan-Meier y prueba Log-Rank",
  ggtheme = theme_minimal(base_size = 13)
)

p_km_energy$plot <- p_km_energy$plot +
  coord_cartesian(xlim = c(0, 760), ylim = c(0, 1.02)) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 12),
    panel.grid.minor = element_blank()
  )

ggsave(
  filename = here("outputs", "figures", "02_kaplan_meier_energy_system.png"),
  plot = p_km_energy$plot,
  width = 10.8,
  height = 5.9,
  dpi = 300
)

# ------------------------------------------------------------
# 9. Kaplan-Meier por exposicion a nieve/hielo, figura limpia
# ------------------------------------------------------------

km_snow <- survfit(surv_obj ~ snow_ice_exposure, data = hazard_data)

p_km_snow <- ggsurvplot(
  km_snow,
  data = hazard_data,
  conf.int = FALSE,
  risk.table = FALSE,
  censor = TRUE,
  pval = TRUE,
  pval.coord = c(25, 0.18),
  legend.title = "Exposición a nieve/hielo",
  legend.labs = c("Baja", "Media", "Alta", "Muy alta"),
  xlab = "Tiempo de observación (días)",
  ylab = "Probabilidad de supervivencia operacional",
  title = "Supervivencia operacional según exposición a nieve/hielo",
  subtitle = "Mayor exposición representa condiciones ambientales más severas",
  ggtheme = theme_minimal(base_size = 13)
)

p_km_snow$plot <- p_km_snow$plot +
  coord_cartesian(xlim = c(0, 760), ylim = c(0, 1.02)) +
  guides(color = guide_legend(nrow = 1), fill = guide_legend(nrow = 1)) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 12),
    panel.grid.minor = element_blank()
  )

ggsave(
  filename = here("outputs", "figures", "03_kaplan_meier_snow_ice.png"),
  plot = p_km_snow$plot,
  width = 11.5,
  height = 5.9,
  dpi = 300
)

# ------------------------------------------------------------
# 10. Pruebas Log-Rank
# ------------------------------------------------------------

logrank_energy <- survdiff(surv_obj ~ energy_system, data = hazard_data)
logrank_snow <- survdiff(surv_obj ~ snow_ice_exposure, data = hazard_data)

logrank_results <- tibble(
  comparacion = c("Sistema energético", "Exposición nieve/hielo"),
  chisq = c(logrank_energy$chisq, logrank_snow$chisq),
  gl = c(
    length(logrank_energy$n) - 1,
    length(logrank_snow$n) - 1
  ),
  p_value = c(
    1 - pchisq(logrank_energy$chisq, length(logrank_energy$n) - 1),
    1 - pchisq(logrank_snow$chisq, length(logrank_snow$n) - 1)
  )
) %>%
  mutate(
    chisq = round(chisq, 3),
    p_value = round(p_value, 5),
    interpretacion = if_else(
      p_value < 0.05,
      "Diferencias estadísticamente relevantes entre curvas",
      "No se observan diferencias estadísticamente relevantes al 5%"
    )
  )

write_csv(logrank_results, here("outputs", "tables", "04_logrank_results.csv"))

# ------------------------------------------------------------
# 11. Modelo Cox de riesgos proporcionales
# ------------------------------------------------------------

cox_model <- coxph(
  Surv(time_days, event) ~
    energy_system +
    snow_ice_exposure +
    age_years +
    battery_margin_days +
    winter_operation,
  data = hazard_data,
  x = TRUE,
  y = TRUE
)

cox_summary <- summary(cox_model)

# Tabla de Hazard Ratios
cox_hr_table <- tidy(
  cox_model,
  exponentiate = TRUE,
  conf.int = TRUE
) %>%
  mutate(
    estimate = round(estimate, 3),
    conf.low = round(conf.low, 3),
    conf.high = round(conf.high, 3),
    p.value = round(p.value, 5),
    interpretacion = case_when(
      estimate > 1 ~ "Aumenta el hazard relativo",
      estimate < 1 ~ "Reduce el hazard relativo",
      TRUE ~ "Sin cambio relativo"
    )
  ) %>%
  rename(
    variable = term,
    hazard_ratio = estimate,
    ic95_inf = conf.low,
    ic95_sup = conf.high,
    p_value = p.value
  )

write_csv(cox_hr_table, here("outputs", "tables", "05_cox_hazard_ratios.csv"))

# ------------------------------------------------------------
# 12. Forest plot de Hazard Ratios, figura limpia
# ------------------------------------------------------------

forest_data <- cox_hr_table %>%
  mutate(
    variable = recode(
      variable,
      "energy_systemEstándar" = "Energía: estándar vs robusto",
      "energy_systemLimitado" = "Energía: limitado vs robusto",
      "snow_ice_exposureMedia" = "Nieve/hielo: media vs baja",
      "snow_ice_exposureAlta" = "Nieve/hielo: alta vs baja",
      "snow_ice_exposureMuy alta" = "Nieve/hielo: muy alta vs baja",
      "age_years" = "Antigüedad del sistema",
      "battery_margin_days" = "Margen de batería",
      "winter_operationSí" = "Operación en invierno crítico"
    ),
    variable = fct_reorder(variable, hazard_ratio)
  )

p_forest <- ggplot(
  forest_data,
  aes(x = hazard_ratio, y = variable)
) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.7) +
  geom_segment(
    aes(x = ic95_inf, xend = ic95_sup, y = variable, yend = variable),
    linewidth = 0.75
  ) +
  geom_point(size = 3.2) +
  scale_x_log10(
    breaks = c(0.3, 0.5, 1, 2, 3, 5, 10, 30),
    labels = c("0.3", "0.5", "1", "2", "3", "5", "10", "30")
  ) +
  coord_cartesian(xlim = c(0.25, 30)) +
  labs(
    title = "Modelo Cox: Hazard Ratios estimados",
    subtitle = "HR > 1 indica mayor riesgo instantáneo de falla; HR < 1 indica menor riesgo relativo",
    x = "Hazard Ratio, escala logarítmica",
    y = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 12),
    panel.grid.minor = element_blank()
  )

ggsave(
  filename = here("outputs", "figures", "04_forest_plot_cox.png"),
  plot = p_forest,
  width = 10.8,
  height = 5.9,
  dpi = 300
)

# ------------------------------------------------------------
# 13. Metricas del modelo
# ------------------------------------------------------------

model_metrics <- tibble(
  metrica = c(
    "AIC",
    "Concordance index C",
    "SE Concordance",
    "Likelihood ratio test p-value",
    "Número de observaciones",
    "Número de eventos"
  ),
  valor = c(
    round(AIC(cox_model), 3),
    round(as.numeric(cox_summary$concordance[1]), 3),
    round(as.numeric(cox_summary$concordance[2]), 3),
    round(as.numeric(cox_summary$logtest["pvalue"]), 5),
    nrow(hazard_data),
    sum(hazard_data$event == 1)
  )
)

write_csv(model_metrics, here("outputs", "tables", "06_model_metrics.csv"))

# Tabla compacta de metricas para presentacion
model_metrics_for_slide <- model_metrics %>%
  filter(metrica %in% c("AIC", "Concordance index C", "Likelihood ratio test p-value")) %>%
  mutate(
    lectura = case_when(
      metrica == "AIC" ~ "Criterio de información para comparar modelos",
      metrica == "Concordance index C" ~ "Capacidad discriminante moderada-alta",
      metrica == "Likelihood ratio test p-value" ~ "El modelo global es estadísticamente relevante",
      TRUE ~ ""
    )
  )

write_csv(model_metrics_for_slide, here("outputs", "tables", "06b_model_metrics_for_slide.csv"))

# ------------------------------------------------------------
# 14. Validacion del supuesto de riesgos proporcionales
# ------------------------------------------------------------

ph_test <- cox.zph(cox_model)

ph_table <- as.data.frame(ph_test$table) %>%
  rownames_to_column("variable") %>%
  as_tibble() %>%
  mutate(
    across(where(is.numeric), ~ round(.x, 5)),
    interpretacion = if_else(
      p < 0.05,
      "Posible violación del supuesto PH",
      "Sin evidencia de violación PH"
    )
  )

write_csv(ph_table, here("outputs", "tables", "07_ph_assumption_test.csv"))

# Grafico diagnostico de Schoenfeld como respaldo tecnico
png(
  filename = here("outputs", "figures", "05_schoenfeld_diagnostics.png"),
  width = 1800,
  height = 1200,
  res = 180
)

par(mfrow = c(2, 3))
plot(ph_test)

dev.off()

# ------------------------------------------------------------
# 15. Tablas compactas para presentacion
# ------------------------------------------------------------

cox_slide_table <- cox_hr_table %>%
  select(variable, hazard_ratio, ic95_inf, ic95_sup, p_value) %>%
  mutate(
    variable = recode(
      variable,
      "energy_systemEstándar" = "Energía: estándar vs robusto",
      "energy_systemLimitado" = "Energía: limitado vs robusto",
      "snow_ice_exposureMedia" = "Nieve/hielo: media vs baja",
      "snow_ice_exposureAlta" = "Nieve/hielo: alta vs baja",
      "snow_ice_exposureMuy alta" = "Nieve/hielo: muy alta vs baja",
      "age_years" = "Antigüedad del sistema",
      "battery_margin_days" = "Margen de batería",
      "winter_operationSí" = "Operación en invierno crítico"
    ),
    hr_ic95 = paste0(
      hazard_ratio,
      " [", ic95_inf, "; ", ic95_sup, "]"
    )
  ) %>%
  select(variable, hr_ic95, p_value)

write_csv(cox_slide_table, here("outputs", "tables", "08_cox_table_for_slide.csv"))

cox_table_for_presentation <- cox_hr_table %>%
  transmute(
    variable = recode(
      variable,
      "energy_systemEstándar" = "Energía: estándar vs robusto",
      "energy_systemLimitado" = "Energía: limitado vs robusto",
      "snow_ice_exposureMedia" = "Nieve/hielo: media vs baja",
      "snow_ice_exposureAlta" = "Nieve/hielo: alta vs baja",
      "snow_ice_exposureMuy alta" = "Nieve/hielo: muy alta vs baja",
      "age_years" = "Antigüedad del sistema",
      "battery_margin_days" = "Margen de batería",
      "winter_operationSí" = "Operación en invierno crítico"
    ),
    HR = round(hazard_ratio, 2),
    IC95 = paste0("[", round(ic95_inf, 2), "; ", round(ic95_sup, 2), "]"),
    p_value = round(p_value, 4),
    lectura = case_when(
      p_value < 0.05 & HR > 1 ~ "Aumenta significativamente el hazard",
      p_value < 0.05 & HR < 1 ~ "Reduce significativamente el hazard",
      p_value >= 0.05 & HR > 1 ~ "Tendencia a mayor hazard, no significativa",
      p_value >= 0.05 & HR < 1 ~ "Tendencia protectora, no significativa",
      TRUE ~ "Sin cambio relevante"
    )
  ) %>%
  arrange(p_value)

write_csv(cox_table_for_presentation, here("outputs", "tables", "09_cox_table_for_presentation.csv"))

# ------------------------------------------------------------
# 16. Mensaje final
# ------------------------------------------------------------

message("============================================================")
message("Analisis Hazard completado correctamente.")
message("Figuras corregidas y exportadas en: outputs/figures")
message("Tablas exportadas en: outputs/tables")
message("Archivos de figuras generados:")
print(list.files(here("outputs", "figures"), full.names = FALSE))
message("Archivos de tablas generados:")
print(list.files(here("outputs", "tables"), full.names = FALSE))
message("============================================================")
