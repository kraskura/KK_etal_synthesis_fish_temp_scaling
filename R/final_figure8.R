# Author: Krista Kraskura 
# Date: april 6 
# theme functions and color assignment functions generated with help from Claude code AI; data not shared only the conceptual functionality of the code. 

library(ggplot2)
library(dplyr)
library(cowplot)
library(colorspace)

base_font <- "Helvetica"

# -------------------------------------------------------
cols.amr<-c(colorBlindness::SteppedSequential5Steps[c(21, 22, 23, 24,  25)]) 
cols.rmr<-c(colorBlindness::SteppedSequential5Steps[c(11, 12, 13,14,  15)]) # "#6B990F" "#A3CC51" "#E5FFB2"
cols.fas<-c(colorBlindness::SteppedSequential5Steps[c(1,2, 3,4, 5)])  # "#990F0F" "#CC5151" "#FFB2B2"
cols.as<-c(colorBlindness::SteppedSequential5Steps[c(16,17,  18, 19, 20)]) 



# read in data: ----------
dat_main<-read.csv(here("./Data_exports/models/scaling_parameters.csv")) # main model
dat_inset1<-read.csv(here("./Data_exports/juvenilemodels/scaling_parameters.csv")) # main model
dat_inset2<-read.csv(here("./Data_exports/adultmodels/scaling_parameters.csv")) # main model

dat_main<- dat_main %>% 
  filter(Temperature == 10 |
           Temperature == 20 |
           Temperature == 30)
dat_inset1<- dat_inset1 %>% 
  filter(Temperature == 10 |
           Temperature == 20 |
           Temperature == 30)
dat_inset2<- dat_inset2 %>% 
  filter(Temperature == 10 |
           Temperature == 20 |
           Temperature == 30)

# --- Color palettes -----------------------------------------
grey_ramp <- colorRampPalette(c("white", "#111111"))(30)

perf_hue_as <- cols.as[3]   # color for this performance panel (category 2)
perf_hue_fas <- cols.fas[3]   # color for this performance panel (category 2)
perf_hue_amr <- cols.amr[3]   # color for this performance panel (category 2)
perf_hue_rmr <- cols.rmr[3]   # color for this performance panel (category 2)

make_perf_ramp <- function(hue, n = 30) {
  light <- {
    rgb_vals <- col2rgb(hue) / 255
    rgb((1 - rgb_vals[1]) * 0.5 + rgb_vals[1],
        (1 - rgb_vals[2]) * 0.5 + rgb_vals[2],
        (1 - rgb_vals[3]) * 0.5 + rgb_vals[3])
  }
  dark <- {
    rgb_vals <- col2rgb(hue) / 255
    rgb(rgb_vals[1] * 0.5,
        rgb_vals[2] * 0.5,
        rgb_vals[3] * 0.5)
  }
  colorRampPalette(c(light, hue, dark))(n)
}


assign_colors <- function(dat, perf_hue_as, perf_hue_fas, 
                               perf_hue_amr, perf_hue_rmr,
                          n_temp) {
  
  ramp_as  <- make_perf_ramp(perf_hue_as,  n = n_temp)
  ramp_fas <- make_perf_ramp(perf_hue_fas, n = n_temp)
  ramp_amr <- make_perf_ramp(perf_hue_amr, n = n_temp)
  ramp_rmr <- make_perf_ramp(perf_hue_rmr, n = n_temp)
  
  dat %>%
    mutate(
      line_color = case_when(
        temp_categ == "ecol_relev" ~ grey_ramp[Temperature],
        temp_categ == "warm" & performance == "AS"  ~ ramp_as[Temperature],
        temp_categ == "warm" & performance == "FAS" ~ ramp_fas[Temperature],
        temp_categ == "warm" & performance == "MMR" ~ ramp_amr[Temperature],
        temp_categ == "warm" & performance == "RMR" ~ ramp_rmr[Temperature]
      )
    )
}



# --- Apply ---
dat_main   <- assign_colors(dat_main,
                            perf_hue_as = perf_hue_as,
                            perf_hue_fas = perf_hue_fas,
                            perf_hue_rmr = perf_hue_rmr,
                            perf_hue_amr = perf_hue_amr,
                            n_temp = 30)
dat_inset1 <- assign_colors(dat_inset1,
                            perf_hue_as = perf_hue_as,
                            perf_hue_fas = perf_hue_fas,
                            perf_hue_rmr = perf_hue_rmr,
                            perf_hue_amr = perf_hue_amr,
                            n_temp = 30)
dat_inset2 <- assign_colors(dat_inset2,
                            perf_hue_as = perf_hue_as,
                            perf_hue_fas = perf_hue_fas,
                            perf_hue_rmr = perf_hue_rmr,
                            perf_hue_amr = perf_hue_amr,
                            n_temp = 30)



# --- X bounds and segment endpoints --------------------------------
x_min <- -5
x_max <-  7

dat_main <- dat_main %>%
  mutate(
    x_start = x_min,
    x_end   = x_max,
    y_start = Intercept + Slope * x_min,
    y_end   = Intercept + Slope * x_max,
    row_id  = row_number()
  )

dat_inset1 <- dat_inset1 %>%
  mutate(
    x_start = x_min,
    x_end   = x_max,
    y_start = Intercept + Slope * x_min,
    y_end   = Intercept + Slope * x_max,
    row_id  = row_number()
  )

dat_inset2 <- dat_inset2 %>%
  mutate(
    x_start = x_min,
    x_end   = x_max,
    y_start = Intercept + Slope * x_min,
    y_end   = Intercept + Slope * x_max,
    row_id  = row_number()
  )

# Shared y limits across all three plots
y_lim <- range(c(dat_main$y_start,   dat_main$y_end,
                  dat_inset1$y_start, dat_inset1$y_end,
                  dat_inset2$y_start, dat_inset2$y_end))

x_lim <- c(x_min, x_max)

# --- Themes ------------------------------------------
panel_theme <- theme_classic(base_size = 11, base_family = base_font) +
  theme(
    axis.text.y      = element_text(size = 11, colour = "black", family = base_font),
    axis.text.x      = element_text(size = 11, colour = "black", family = base_font),
    axis.line.y      = element_line(colour = "black", linewidth = 0.5),
    axis.line.x      = element_line(colour = "black", linewidth = 0.5),
    axis.ticks.y     = element_line(linewidth = 0.5),
    axis.ticks.x     = element_line(linewidth = 0),
    axis.title.y     = element_text(size = 11, family = base_font),
    axis.title.x     = element_text(size = 11, family = base_font),
    panel.border     = element_rect(linetype = "solid", fill = NA, colour = "black"),
    plot.title       = element_blank(),
    legend.position  = "none",
    plot.background  = element_rect(fill = "white", color = NA)
  )

inset_theme <- theme_classic(base_size = 8, base_family = base_font) +
  theme(
    legend.position  = "none",
    axis.title       = element_blank(),
    axis.text        = element_blank(),
    axis.ticks       = element_blank(),
    axis.line.y      = element_line(colour = "black", linewidth = 0.4),
    axis.line.x      = element_line(colour = "black", linewidth = 0.4),
    panel.border     = element_rect(linetype = "solid", fill = NA, colour = "black"),
    plot.title       = element_text(size = 8, face = "bold", hjust = 0.5,
                                    family = base_font, margin = margin(b = 1)),
    plot.background  = element_rect(fill = "transparent", color = NA,
                                    linewidth = 0.4),
    plot.margin      = margin(3, 3, 3, 3),
    aspect.ratio     = 1
  )


# ******************************************************************
# --- Main plot AS OPTIMAL ---------------------------------------------
performance <- "AS"
temp_category <- "ecol_relev"
plot_main <- ggplot(dat_main[dat_main$performance == performance & 
                            dat_main$temp_categ == temp_category,],
                    aes(x     = x_start,
                        xend  = x_end,
                        y     = y_start,
                        yend  = y_end,
                        color = line_color,
                        group = row_id)) +
  geom_segment(alpha =1, linewidth = 1.1) +
  scale_color_identity() +
  scale_x_continuous(limits = x_lim) +
  coord_cartesian(ylim = y_lim) +
  geom_text(data        = data.frame(x = x_max, y = max(y_lim),
                                     label = "AS: Optimal Temperatures"),
            aes(x = x, y = y, label = label),
            hjust        = 1.05,
            vjust        = 1.5,
            size         = 3.5,
            family       = base_font,
            fontface     = "bold",
            inherit.aes  = FALSE) +
  labs(x = expression(italic(ln)~Mass~(g)), y = expression(italic(ln)~MR~(mgO[2]~h^-1))) +
  panel_theme

# --- Inset 1 --- juveniles
plot_inset1 <- ggplot(dat_inset1[dat_inset1$performance == performance & 
                              dat_inset1$temp_categ == temp_category,],
                      aes(x     = x_start,
                          xend  = x_end,
                          y     = y_start,
                          yend  = y_end,
                          color = line_color,
                          group = row_id)) +
  geom_segment(alpha =1, linewidth = 0.6) +
  scale_color_identity() +
  scale_x_continuous(limits = x_lim) +
  coord_cartesian(ylim = y_lim) +
  labs(title = "Juveniles") +
  inset_theme

# --- Inset 2 --- adults
plot_inset2 <- ggplot(dat_inset2[dat_inset2$performance == performance & 
                              dat_inset2$temp_categ == temp_category,],
                      aes(x     = x_start,
                          xend  = x_end,
                          y     = y_start,
                          yend  = y_end,
                          color = line_color,
                          group = row_id)) +
  geom_segment(alpha = 1, linewidth = 0.6) +
  scale_color_identity() +
  scale_x_continuous(limits = x_lim) +
  coord_cartesian(ylim = y_lim) +
  labs(title = "Adults") +
  inset_theme


# --- Compose ---
composed_as_er <- ggdraw(plot_main) +
  draw_plot(plot_inset1,
            x = 0.45,
            y = 0.17,
            width = 0.32,
            height = 0.32) +
  draw_plot(plot_inset2,
            x = 0.69,
            y = 0.17,
            width = 0.32,
            height = 0.32)

# composed_as_er


# ******************************************************************
# ******************************************************************
# --- Main plot AS WAMR ---------------------------------------------
performance <- "AS"
temp_category <- "warm"
plot_main <- ggplot(dat_main[dat_main$performance == performance & 
                            dat_main$temp_categ == temp_category,],
                    aes(x     = x_start,
                        xend  = x_end,
                        y     = y_start,
                        yend  = y_end,
                        color = line_color,
                        group = row_id)) +
  geom_segment(alpha =1, linewidth = 1.1) +
  scale_color_identity() +
  scale_x_continuous(limits = x_lim) +
  coord_cartesian(ylim = y_lim) +
  geom_text(data        = data.frame(x = x_max, y = max(y_lim),
                                     label = "AS: Warm Temperatures"),
            aes(x = x, y = y, label = label),
            hjust        = 1.05,
            vjust        = 1.5,
            size         = 3.5,
            family       = base_font,
            fontface     = "bold",
            inherit.aes  = FALSE) +
  labs(x = expression(italic(ln)~Mass~(g)), y = expression(italic(ln)~MR~(mgO[2]~h^-1))) +
  panel_theme

# --- Inset 1 --- juveniles
plot_inset1 <- ggplot(dat_inset1[dat_inset1$performance == performance & 
                              dat_inset1$temp_categ == temp_category,],
                      aes(x     = x_start,
                          xend  = x_end,
                          y     = y_start,
                          yend  = y_end,
                          color = line_color,
                          group = row_id)) +
  geom_segment(alpha =1, linewidth = 0.6) +
  scale_color_identity() +
  scale_x_continuous(limits = x_lim) +
  coord_cartesian(ylim = y_lim) +
  labs(title = "Juveniles") +
  inset_theme

# --- Inset 2 --- adults
plot_inset2 <- ggplot(dat_inset2[dat_inset2$performance == performance & 
                              dat_inset2$temp_categ == temp_category,],
                      aes(x     = x_start,
                          xend  = x_end,
                          y     = y_start,
                          yend  = y_end,
                          color = line_color,
                          group = row_id)) +
  geom_segment(alpha = 1, linewidth = 0.6) +
  scale_color_identity() +
  scale_x_continuous(limits = x_lim) +
  coord_cartesian(ylim = y_lim) +
  labs(title = "Adults") +
  inset_theme


# --- Compose ---
composed_as_w <- ggdraw(plot_main) +
  draw_plot(plot_inset1,
            x = 0.45,
            y = 0.17,
            width = 0.32,
            height = 0.32) +
  draw_plot(plot_inset2,
            x = 0.69,
            y = 0.17,
            width = 0.32,
            height = 0.32)

# composed_as_w


# ******************************************************************
# ******************************************************************
# --- Main plot FAS WAMR ---------------------------------------------
performance <- "FAS"
temp_category <- "warm"
plot_main <- ggplot(dat_main[dat_main$performance == performance & 
                            dat_main$temp_categ == temp_category,],
                    aes(x     = x_start,
                        xend  = x_end,
                        y     = y_start,
                        yend  = y_end,
                        color = line_color,
                        group = row_id)) +
  geom_segment(alpha =1, linewidth = 1.1) +
  scale_color_identity() +
  scale_x_continuous(limits = x_lim) +
  coord_cartesian(ylim = y_lim) +
  geom_text(data        = data.frame(x = x_max, y = max(y_lim),
                                     label = "FAS: Warm Temperatures"),
            aes(x = x, y = y, label = label),
            hjust        = 1.05,
            vjust        = 1.5,
            size         = 3.5,
            family       = base_font,
            fontface     = "bold",
            inherit.aes  = FALSE) +
  labs(x = expression(italic(ln)~Mass~(g)), y = expression(italic(ln)~MR~(mgO[2]~h^-1))) +
  panel_theme

# --- Inset 1 --- juveniles
plot_inset1 <- ggplot(dat_inset1[dat_inset1$performance == performance & 
                              dat_inset1$temp_categ == temp_category,],
                      aes(x     = x_start,
                          xend  = x_end,
                          y     = y_start,
                          yend  = y_end,
                          color = line_color,
                          group = row_id)) +
  geom_segment(alpha =1, linewidth = 0.6) +
  scale_color_identity() +
  scale_x_continuous(limits = x_lim) +
  coord_cartesian(ylim = y_lim) +
  labs(title = "Juveniles") +
  inset_theme

# --- Inset 2 --- adults
plot_inset2 <- ggplot(dat_inset2[dat_inset2$performance == performance & 
                              dat_inset2$temp_categ == temp_category,],
                      aes(x     = x_start,
                          xend  = x_end,
                          y     = y_start,
                          yend  = y_end,
                          color = line_color,
                          group = row_id)) +
  geom_segment(alpha = 1, linewidth = 0.6) +
  scale_color_identity() +
  scale_x_continuous(limits = x_lim) +
  coord_cartesian(ylim = y_lim) +
  labs(title = "Adults") +
  inset_theme


# --- Compose ---
composed_fas_w <- ggdraw(plot_main) +
  draw_plot(plot_inset1,
            x = 0.45,
            y = 0.17,
            width = 0.32,
            height = 0.32) +
  draw_plot(plot_inset2,
            x = 0.69,
            y = 0.17,
            width = 0.32,
            height = 0.32)

# composed_fas_w

# ******************************************************************
# ******************************************************************
# --- Main plot FAS OPTIMAL ---------------------------------------------
performance <- "FAS"
temp_category <- "ecol_relev"
plot_main <- ggplot(dat_main[dat_main$performance == performance & 
                            dat_main$temp_categ == temp_category,],
                    aes(x     = x_start,
                        xend  = x_end,
                        y     = y_start,
                        yend  = y_end,
                        color = line_color,
                        group = row_id)) +
  geom_segment(alpha =1, linewidth = 1.1) +
  scale_color_identity() +
  scale_x_continuous(limits = x_lim) +
  coord_cartesian(ylim = y_lim) +
  geom_text(data        = data.frame(x = x_max, y = max(y_lim),
                                     label = "FAS: Optimal Temperatures"),
            aes(x = x, y = y, label = label),
            hjust        = 1.05,
            vjust        = 1.5,
            size         = 3.5,
            family       = base_font,
            fontface     = "bold",
            inherit.aes  = FALSE) +
  labs(x = expression(italic(ln)~Mass~(g)), y = expression(italic(ln)~MR~(mgO[2]~h^-1))) +
  panel_theme

# --- Inset 1 --- juveniles
plot_inset1 <- ggplot(dat_inset1[dat_inset1$performance == performance & 
                              dat_inset1$temp_categ == temp_category,],
                      aes(x     = x_start,
                          xend  = x_end,
                          y     = y_start,
                          yend  = y_end,
                          color = line_color,
                          group = row_id)) +
  geom_segment(alpha =1, linewidth = 0.6) +
  scale_color_identity() +
  scale_x_continuous(limits = x_lim) +
  coord_cartesian(ylim = y_lim) +
  labs(title = "Juveniles") +
  inset_theme

# --- Inset 2 --- adults
plot_inset2 <- ggplot(dat_inset2[dat_inset2$performance == performance & 
                              dat_inset2$temp_categ == temp_category,],
                      aes(x     = x_start,
                          xend  = x_end,
                          y     = y_start,
                          yend  = y_end,
                          color = line_color,
                          group = row_id)) +
  geom_segment(alpha = 1, linewidth = 0.6) +
  scale_color_identity() +
  scale_x_continuous(limits = x_lim) +
  coord_cartesian(ylim = y_lim) +
  labs(title = "Adults") +
  inset_theme


# --- Compose ---
composed_fas_er <- ggdraw(plot_main) +
  draw_plot(plot_inset1,
            x = 0.45,
            y = 0.17,
            width = 0.32,
            height = 0.32) +
  draw_plot(plot_inset2,
            x = 0.69,
            y = 0.17,
            width = 0.32,
            height = 0.32)

# composed_fas_er



# ******************************************************************
# ******************************************************************
# --- Main plot MMR OPTIMAL ---------------------------------------------
performance <- "MMR"
temp_category <- "ecol_relev"
plot_main <- ggplot(dat_main[dat_main$performance == performance & 
                            dat_main$temp_categ == temp_category,],
                    aes(x     = x_start,
                        xend  = x_end,
                        y     = y_start,
                        yend  = y_end,
                        color = line_color,
                        group = row_id)) +
  geom_segment(alpha =1, linewidth = 1.1) +
  scale_color_identity() +
  scale_x_continuous(limits = x_lim) +
  coord_cartesian(ylim = y_lim) +
  geom_text(data        = data.frame(x = x_max, y = max(y_lim),
                                     label = "MMR: Optimal Temperatures"),
            aes(x = x, y = y, label = label),
            hjust        = 1.05,
            vjust        = 1.5,
            size         = 3.5,
            family       = base_font,
            fontface     = "bold",
            inherit.aes  = FALSE) +
  labs(x = expression(italic(ln)~Mass~(g)), y = expression(italic(ln)~MR~(mgO[2]~h^-1))) +
  panel_theme

# --- Inset 1 --- juveniles
plot_inset1 <- ggplot(dat_inset1[dat_inset1$performance == performance & 
                              dat_inset1$temp_categ == temp_category,],
                      aes(x     = x_start,
                          xend  = x_end,
                          y     = y_start,
                          yend  = y_end,
                          color = line_color,
                          group = row_id)) +
  geom_segment(alpha =1, linewidth = 0.6) +
  scale_color_identity() +
  scale_x_continuous(limits = x_lim) +
  coord_cartesian(ylim = y_lim) +
  labs(title = "Juveniles") +
  inset_theme

# --- Inset 2 --- adults
plot_inset2 <- ggplot(dat_inset2[dat_inset2$performance == performance & 
                              dat_inset2$temp_categ == temp_category,],
                      aes(x     = x_start,
                          xend  = x_end,
                          y     = y_start,
                          yend  = y_end,
                          color = line_color,
                          group = row_id)) +
  geom_segment(alpha = 1, linewidth = 0.6) +
  scale_color_identity() +
  scale_x_continuous(limits = x_lim) +
  coord_cartesian(ylim = y_lim) +
  labs(title = "Adults") +
  inset_theme


# --- Compose ---
composed_amr_er <- ggdraw(plot_main) +
  draw_plot(plot_inset1,
            x = 0.45,
            y = 0.17,
            width = 0.32,
            height = 0.32) +
  draw_plot(plot_inset2,
            x = 0.69,
            y = 0.17,
            width = 0.32,
            height = 0.32)

# composed_amr_er


# ******************************************************************
# ******************************************************************
# --- Main plot MMR WARM ---------------------------------------------
performance <- "MMR"
temp_category <- "warm"
plot_main <- ggplot(dat_main[dat_main$performance == performance & 
                            dat_main$temp_categ == temp_category,],
                    aes(x     = x_start,
                        xend  = x_end,
                        y     = y_start,
                        yend  = y_end,
                        color = line_color,
                        group = row_id)) +
  geom_segment(alpha =1, linewidth = 1.1) +
  scale_color_identity() +
  scale_x_continuous(limits = x_lim) +
  coord_cartesian(ylim = y_lim) +
  geom_text(data        = data.frame(x = x_max, y = max(y_lim),
                                     label = "MMR: Warm Temperatures"),
            aes(x = x, y = y, label = label),
            hjust        = 1.05,
            vjust        = 1.5,
            size         = 3.5,
            family       = base_font,
            fontface     = "bold",
            inherit.aes  = FALSE) +
  labs(x = expression(italic(ln)~Mass~(g)), y = expression(italic(ln)~MR~(mgO[2]~h^-1))) +  panel_theme

# --- Inset 1 --- juveniles
plot_inset1 <- ggplot(dat_inset1[dat_inset1$performance == performance & 
                              dat_inset1$temp_categ == temp_category,],
                      aes(x     = x_start,
                          xend  = x_end,
                          y     = y_start,
                          yend  = y_end,
                          color = line_color,
                          group = row_id)) +
  geom_segment(alpha =1, linewidth = 0.6) +
  scale_color_identity() +
  scale_x_continuous(limits = x_lim) +
  coord_cartesian(ylim = y_lim) +
  labs(title = "Juveniles") +
  inset_theme

# --- Inset 2 --- adults
plot_inset2 <- ggplot(dat_inset2[dat_inset2$performance == performance & 
                              dat_inset2$temp_categ == temp_category,],
                      aes(x     = x_start,
                          xend  = x_end,
                          y     = y_start,
                          yend  = y_end,
                          color = line_color,
                          group = row_id)) +
  geom_segment(alpha = 1, linewidth = 0.6) +
  scale_color_identity() +
  scale_x_continuous(limits = x_lim) +
  coord_cartesian(ylim = y_lim) +
  labs(title = "Adults") +
  inset_theme


# --- Compose ---
composed_amr_w <- ggdraw(plot_main) +
  draw_plot(plot_inset1,
            x = 0.45,
            y = 0.17,
            width = 0.32,
            height = 0.32) +
  draw_plot(plot_inset2,
            x = 0.69,
            y = 0.17,
            width = 0.32,
            height = 0.32)

# composed_amr_w


# ******************************************************************
# ******************************************************************
# --- Main plot RMR EOPTIMAL ---------------------------------------------
performance <- "RMR"
temp_category <- "ecol_relev"
plot_main <- ggplot(dat_main[dat_main$performance == performance & 
                            dat_main$temp_categ == temp_category,],
                    aes(x     = x_start,
                        xend  = x_end,
                        y     = y_start,
                        yend  = y_end,
                        color = line_color,
                        group = row_id)) +
  geom_segment(alpha =1, linewidth = 1.1) +
  scale_color_identity() +
  scale_x_continuous(limits = x_lim) +
  coord_cartesian(ylim = y_lim) +
  geom_text(data        = data.frame(x = x_max, y = max(y_lim),
                                     label = "RMR: Optimal Temperatures"),
            aes(x = x, y = y, label = label),
            hjust        = 1.05,
            vjust        = 1.5,
            size         = 3.5,
            family       = base_font,
            fontface     = "bold",
            inherit.aes  = FALSE) +
  labs(x = expression(italic(ln)~Mass~(g)), y = expression(italic(ln)~MR~(mgO[2]~h^-1))) +  panel_theme

# --- Inset 1 --- juveniles
plot_inset1 <- ggplot(dat_inset1[dat_inset1$performance == performance & 
                              dat_inset1$temp_categ == temp_category,],
                      aes(x     = x_start,
                          xend  = x_end,
                          y     = y_start,
                          yend  = y_end,
                          color = line_color,
                          group = row_id)) +
  geom_segment(alpha =1, linewidth = 0.6) +
  scale_color_identity() +
  scale_x_continuous(limits = x_lim) +
  coord_cartesian(ylim = y_lim) +
  labs(title = "Juveniles") +
  inset_theme

# --- Inset 2 --- adults
plot_inset2 <- ggplot(dat_inset2[dat_inset2$performance == performance & 
                              dat_inset2$temp_categ == temp_category,],
                      aes(x     = x_start,
                          xend  = x_end,
                          y     = y_start,
                          yend  = y_end,
                          color = line_color,
                          group = row_id)) +
  geom_segment(alpha = 1, linewidth = 0.6) +
  scale_color_identity() +
  scale_x_continuous(limits = x_lim) +
  coord_cartesian(ylim = y_lim) +
  labs(title = "Adults") +
  inset_theme


# --- Compose ---
composed_rmr_er <- ggdraw(plot_main) +
  draw_plot(plot_inset1,
            x = 0.45,
            y = 0.17,
            width = 0.32,
            height = 0.32) +
  draw_plot(plot_inset2,
            x = 0.69,
            y = 0.17,
            width = 0.32,
            height = 0.32)

# composed_rmr_er


# ******************************************************************
# ******************************************************************
# --- Main plot RMR WAMR ---------------------------------------------
performance <- "RMR"
temp_category <- "warm"
plot_main <- ggplot(dat_main[dat_main$performance == performance & 
                            dat_main$temp_categ == temp_category,],
                    aes(x     = x_start,
                        xend  = x_end,
                        y     = y_start,
                        yend  = y_end,
                        color = line_color,
                        group = row_id)) +
  geom_segment(alpha =1, linewidth = 1.1) +
  scale_color_identity() +
  scale_x_continuous(limits = x_lim) +
  coord_cartesian(ylim = y_lim) +
  geom_text(data        = data.frame(x = x_max, y = max(y_lim),
                                     label = "RMR: Warm Temperatures"),
            aes(x = x, y = y, label = label),
            hjust        = 1.05,
            vjust        = 1.5,
            size         = 3.5,
            family       = base_font,
            fontface     = "bold",
            inherit.aes  = FALSE) +
  labs(x = expression(italic(ln)~Mass~(g)), y = expression(italic(ln)~MR~(mgO[2]~h^-1))) +
  panel_theme

# --- Inset 1 --- juveniles
plot_inset1 <- ggplot(dat_inset1[dat_inset1$performance == performance & 
                              dat_inset1$temp_categ == temp_category,],
                      aes(x     = x_start,
                          xend  = x_end,
                          y     = y_start,
                          yend  = y_end,
                          color = line_color,
                          group = row_id)) +
  geom_segment(alpha =1, linewidth = 0.6) +
  scale_color_identity() +
  scale_x_continuous(limits = x_lim) +
  coord_cartesian(ylim = y_lim) +
  labs(title = "Juveniles") +
  inset_theme

# --- Inset 2 --- adults
plot_inset2 <- ggplot(dat_inset2[dat_inset2$performance == performance & 
                              dat_inset2$temp_categ == temp_category,],
                      aes(x     = x_start,
                          xend  = x_end,
                          y     = y_start,
                          yend  = y_end,
                          color = line_color,
                          group = row_id)) +
  geom_segment(alpha = 1, linewidth = 0.6) +
  scale_color_identity() +
  scale_x_continuous(limits = x_lim) +
  coord_cartesian(ylim = y_lim) +
  labs(title = "Adults") +
  inset_theme


# --- Compose ---
composed_rmr_w <- ggdraw(plot_main) +
  draw_plot(plot_inset1,
            x = 0.45,
            y = 0.17,
            width = 0.32,
            height = 0.32) +
  draw_plot(plot_inset2,
            x = 0.69,
            y = 0.17,
            width = 0.32,
            height = 0.32)

# composed_rmr_w

# SAVE ---------
composed<-cowplot::plot_grid(composed_amr_er,
                   composed_amr_w,
                   composed_rmr_er,
                   composed_rmr_w,
                   composed_as_er,
                   composed_as_w,
                   composed_fas_er,
                   composed_fas_w, 
                   nrow = 4, ncol = 2) 
ggsave(here("./Figures/Figure8.png"), composed,  width = 6.5, height = 12, dpi = 300)
