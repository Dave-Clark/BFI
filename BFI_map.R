library(rnrfa)
library(ggplot2)
library(data.table)
library(rnaturalearth)
library(rnaturalearthhires)
library(sf)
library(viridis)

x <- catalogue()
setDT(x)
bfiDat <- x[, .(latitude, longitude, bfihost)]

bfiDat <- st_as_sf(bfiDat, coords = c("longitude", "latitude"), crs = 4326)

uk <- ne_countries(scale = "large", returnclass="sf",
  country = "United Kingdom")

bfiMap <- ggplot(data = uk) +
  geom_sf() +
  geom_sf(data = bfiDat, aes(col = bfihost), size = 3) +
  scale_colour_viridis(name = "Base flow index") +
  theme(axis.text = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    panel.grid = element_blank()) +
  coord_sf(xlim = c(-10, 1.53066))

ggsave("/media/storage/BFI_paper/figures/BFI_map.pdf", bfiMap, height = 9,
  width = 7, device = "pdf")
