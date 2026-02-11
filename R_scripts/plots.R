library(ggplot2)
library(scales)
library(patchwork)
library(viridis)

my_theme <- theme_classic() +
  theme(axis.line = element_line(color = "grey65"),
        axis.text = element_text(size = 14, colour = "grey50"),
        axis.title = element_text(size = 14, colour = "grey60"),
        plot.subtitle = element_text(size = 14, colour = "grey55"))

theme_set(my_theme)

CPC_clr  <- "#6E291C"
KFCS_clr <- "#F08C61"
NMC_all  <- "#59979C"
NMC_PCR  <- "#B2CABD"

DENV_1_4    <- c("#E35A43", "#F5CE5A", "#7B92DC", "#ED47C9") # Palette 824
untyped     <- "#383838"
subclinical <- c("grey85")

source("./R_scripts/plots_fig1.R")
source("./R_scripts/plots_fig2.R")
source("./R_scripts/plots_fig3.R")
source("./R_scripts/plots_inference.R")
source("./R_scripts/plots_supplementary.R")
