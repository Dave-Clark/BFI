hydroPlot <- ggplot(hydroData[! is.na(value) & variable == "meanHydrophob"],
    aes(x = bfi, y = value, fill = geol, shape = .id)) +
  geom_point(size = 4, alpha = 0.7) +
  scale_shape_manual(values = c(21, 22)) +
  scale_fill_manual("Geology", values = c("white", "grey", "darkseagreen3")) +
  guides(fill = guide_legend(override.aes=list(shape=21))) +
  geom_ribbon(data = allHydroPreds[variable == "meanHydrophob"],
    aes(x = bfi, y = value, ymin = lowCI, ymax = uppCI),
    alpha = 0.4, fill = "grey", colour = NA) +
  geom_line(data = allHydroPreds[variable == "meanHydrophob"],
    aes(x = bfi, y = value, linetype = .id)) +
  labs(x = "Base flow index", y = "Average hydrophobicity", shape = "",
    linetype = "") +
  theme_bw() +
  theme(axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    strip.text.x = element_text(size = 14),
    plot.title = element_text(size = 20),
    panel.grid = element_blank())

chargePlot <- ggplot(hydroData[! is.na(value) & variable == "meanCharge"],
    aes(x = bfi, y = value, fill = geol, shape = .id)) +
  geom_point(size = 4, alpha = 0.7) +
  scale_shape_manual(values = c(21, 22)) +
  scale_fill_manual("Geology", values = c("white", "grey", "darkseagreen3")) +
  guides(fill = guide_legend(override.aes=list(shape=21))) +
  geom_ribbon(data = allHydroPreds[variable == "meanCharge"],
    aes(x = bfi, y = value, ymin = lowCI, ymax = uppCI),
    alpha = 0.4, fill = "grey", colour = NA) +
  geom_line(data = allHydroPreds[variable == "meanCharge"],
    aes(x = bfi, y = value, linetype = .id)) +
  labs(x = "Base flow index", y = "Average net charge", shape = "",
    linetype = "") +
  theme_bw() +
  theme(axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    strip.text.x = element_text(size = 14),
    plot.title = element_text(size = 20),
    panel.grid = element_blank())

caiPlot <- ggplot(hydroData[! is.na(value) & variable == "meanCAI"],
    aes(x = bfi, y = value, fill = geol, shape = .id)) +
  geom_point(size = 4, alpha = 0.7) +
  scale_shape_manual(values = c(21, 22)) +
  scale_fill_manual("Geology", values = c("white", "grey", "darkseagreen3")) +
  guides(fill = guide_legend(override.aes=list(shape=21))) +
  geom_ribbon(data = allHydroPreds[variable == "meanCAI"],
    aes(x = bfi, y = value, ymin = lowCI, ymax = uppCI),
    alpha = 0.4, fill = "grey", colour = NA) +
  geom_line(data = allHydroPreds[variable == "meanCAI"],
    aes(x = bfi, y = value, linetype = .id)) +
  labs(x = "Base flow index", y = "Average CAI", shape = "",
    linetype = "") +
  theme_bw() +
  theme(axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    strip.text.x = element_text(size = 14),
    plot.title = element_text(size = 20),
    panel.grid = element_blank())

gcPlot <- ggplot(hydroData[! is.na(value) & variable == "meanGC"],
    aes(x = bfi, y = value, fill = geol, shape = .id)) +
  geom_point(size = 4, alpha = 0.7) +
  scale_shape_manual(values = c(21, 22)) +
  scale_fill_manual("Geology", values = c("white", "grey", "darkseagreen3")) +
  guides(fill = guide_legend(override.aes=list(shape=21))) +
  geom_ribbon(data = allHydroPreds[variable == "meanGC"],
    aes(x = bfi, y = value, ymin = lowCI, ymax = uppCI),
    alpha = 0.4, fill = "grey", colour = NA) +
  geom_line(data = allHydroPreds[variable == "meanGC"],
    aes(x = bfi, y = value, linetype = .id)) +
  labs(x = "Base flow index", y = "Average GC content (%)", shape = "",
    linetype = "") +
  theme_bw() +
  theme(axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    strip.text.x = element_text(size = 14),
    plot.title = element_text(size = 20),
    panel.grid = element_blank())

fig <- caiPlot + hydroPlot +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A")

ggsave("/media/storage/Dropbox/cv/funding_applications/2021_NERC_NEOF_Pilot/Fig_1.png", fig, height = 4, width = 9, device = "png", dpi = 300)
