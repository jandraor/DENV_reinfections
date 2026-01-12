library(ggplot2)
library(scales)
library(patchwork)
library(viridis)

my_theme <- theme_classic() +
  theme(axis.line = element_line(color = "grey65"),
        axis.text = element_text(size = 12, colour = "grey50"),
        axis.title = element_text(size = 14, colour = "grey60"),
        plot.subtitle = element_text(colour = "grey55"))

theme_set(my_theme)

CPC_clr  <- "#6E291C"
KFCS_clr <- "#F08C61"
NMC_all  <- "#59979C"
NMC_PCR  <- "#B2CABD"

DENV_1_4    <- c("#9AC4C5", "#7CA3CD", "#7260A8", "#ECEDE0") # Palette 880
untyped     <- "#154B66"
subclinical <- c("grey70")

source("./R_scripts/plots_fig1.R")
source("./R_scripts/plots_fig2.R")
source("./R_scripts/plots_fig3.R")
