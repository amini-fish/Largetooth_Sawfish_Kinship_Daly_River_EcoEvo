# Packages
library(dplyr)
library(tidygraph)
library(ggraph)
library(igraph)
library(ggplot2)
library(stringr)

set.seed(42)

kin_pairs <- read.csv("C:/Users/samue/Desktop/Research/Honours/Chapter_2_Relatedness_paper/Kinference/kinference_kin.csv")

nodes_df <- read.csv("C:/Users/samue/Desktop/Research/Honours/10658_prelim_report/Pristis_metadata_SRA.csv")

nodes_df <- nodes_df %>%
  dplyr::filter(!year_caught == c("2021"))

table(nodes_df$year_caught)

kin_pairs

table(kin_pairs$rel, kin_pairs$pw_yr) ## sanity checked - no cross cohort FSP fuck yeah

glimpse(kin_pairs)

kin_pairs$rel


# --- 1) Keep only FSP / HSP edges; keep ALL individuals as nodes
edges_keep <- kin_pairs %>%
  dplyr::mutate(rel = str_to_upper(rel)) %>%
  dplyr::filter(!rel == ("UNRELATED")) %>%
  transmute(
    from = id_1,
    to   = id_2,
    kin_type  = factor(rel),         # keep if you want linetype by location
  )

edges_keep

# --- 2) Nodes (all individuals) ---
nodes_keep <- nodes_df %>%
  dplyr::distinct(id, .keep_all = TRUE) %>%
  dplyr::mutate(Year = as.factor(year_caught),
         sex  = factor(sex))

# --- 3) Graph with isolates kept ---
tg <- tbl_graph(nodes = nodes_keep, edges = edges_keep, directed = FALSE) %>%
  activate(nodes) %>% arrange(Year, id) %>%   # group around circle by Year
  activate(edges)

# --- 4) Layout for year labels (node positions only) ---
pos <- ggraph::create_layout(tg, layout = "linear", circular = TRUE)

year_labels <- pos %>%
  as.data.frame() %>%
  group_by(Year) %>%
  summarise(x = mean(x), y = mean(y), .groups = "drop") %>%
  mutate(r = sqrt(x^2 + y^2),
         x_lab = x * 1.12 / ifelse(r == 0, 1, r),
         y_lab = y * 1.12 / ifelse(r == 0, 1, r))



glimpse(tg)

# --- 5) Plot from the GRAPH (edge attrs preserved), positions manual ---
plot <- ggraph(tg, layout = "linear", circular = TRUE) +
  
  geom_edge_arc(aes(edge_colour = kin_type),
                lineend = "round",
                width = 0.8,
                show.legend = TRUE) +
  
  scale_edge_colour_manual(
    values = c("HALF-SIBS" = "skyblue", "FULL-SIBS" = "navyblue"), 
    name = "Kin Type"
  ) +
  
  geom_node_point(aes(shape = sex, 
                      fill = Site_ID), 
                  colour = "black", 
                  size = 4, 
                  stroke = 0.6) +
  
  scale_shape_manual(values = c("F" = 21, "M" = 24, "Unknown" = 23)) +
  
  scale_fill_grey(name = "Location", drop = FALSE) + 
  
  guides(
    fill = guide_legend(override.aes = list(shape = 21, size = 5, colour = "black", stroke = 0.6)),
    edge_colour = guide_legend(override.aes = list(edge_width = 2), order = 2),
    shape       = guide_legend(order = 3)) +
  
  coord_fixed(xlim = c(-1.1, 1.1), ylim = c(-1.1, 1.1)) +
  
  theme_void() +
  
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank())

print(plot)

devEMF::emf("C:/Users/samue/Desktop/CKMRsim_Network2.emf", width = 10, height = 10)  
print(plot)
dev.off()

library(ggnewscale)

## add year facets 

pos_df <- as.data.frame(pos) |> dplyr::mutate(ord = dplyr::row_number())
N <- nrow(pos_df)
R <- max(sqrt(pos_df$x^2 + pos_df$y^2), na.rm = TRUE)

# 2) One arc segment per Year (contiguous groups assumed)
year_strips <- pos_df |>
  group_by(Year) |>
  summarise(start_idx = min(ord) - 1, end_idx = max(ord), .groups = "drop") |>
  mutate(
    start = 2*pi*start_idx/N,      # radians
    end   = 2*pi*end_idx/N,
    mid   = (start + end)/2,
    r0    = R*1.06,                # inner radius of strip
    r1    = R*1.16,                # outer radius of strip
    lab_r = R*1.20,                # label radius
    x_lab = lab_r*cos(mid),
    y_lab = lab_r*sin(mid),
    ang   = (mid*180/pi) %% 360,
    hjust = ifelse(ang > 90 & ang < 270, 1, 0),
    ang_text = ifelse(ang > 90 & ang < 270, ang + 180, ang)
  )


kinf_circle_plot <- plot + ggnewscale::new_scale_fill() +
  ggforce::geom_arc_bar(
    data = year_strips,
    aes(x0 = 0, y0 = 0, r0 = r0, r = r1, start = start, end = end, fill = Year),
    inherit.aes = FALSE,
    colour = "black", 
    show.legend = TRUE) + 
  scale_fill_brewer(palette = "YlOrRd") +
  theme(legend.position = "none")

print(kinf_circle_plot)

devEMF::emf("C:/Users/samue/Desktop/kinference_Network.emf", width = 10, height = 10)  
print(plot)
dev.off()

joined_plot <- ggpubr::ggarrange(rel.hist, plot, 
                                 labels = c("A", "B"),
                                 legend = "bottom",
                                 nrow = 2,
                                 ncol = 1,
                                 common.legend = TRUE)





#geom_text(data = year_labels, aes(x = x_lab, y = y_lab, label = Year), 
#           fontface = "bold", size = 4)

dev.off()


## Kin only network: 

sibs <- read.csv("C:/Users/samue/Desktop/Research/Honours/Chapter_2_Relatedness_paper/Kinference/kinference_kin.csv")

meta <- read.csv("C:/Users/samue/Desktop/Research/Honours/10658_prelim_report/Pristis_metadata_SRA.csv")


names(sibs)

names(sibs)

dim(sibs)

kinNWdata <- sibs %>%
  dplyr::select(id_1, id_2, rel, pw_yr, pw_loc)

#This makes our data frame from which the pariwise network plot between select individuals will be drawn 

network <- igraph::graph_from_data_frame(d = kinNWdata, directed = TRUE); print(network) # works

df <- data.frame(id = igraph::V(network)$name)
df

vertices <- dplyr::left_join(df, meta, by = "id") %>%
  dplyr::select(id, sex, birth_year, Site_ID)

vertices <- as_tibble(vertices, what = "vertices")

vertices

network <- igraph::graph_from_data_frame(d = kinNWdata, directed = TRUE,
                                         vertices = vertices ) 

layout <- ggraph::create_layout(network, layout = 'igraph', 
                                circular = FALSE, algorithm = 'fr')

## Plot the network ##

kin_network <- ggraph::ggraph(network, layout = layout) + 
  
  ggraph::geom_edge_link(
    aes(
      edge_colour = factor(kinNWdata$rel), 
      edge_linetype = kinNWdata$pw_loc_word), 
    width = 2) +
  
  ggraph::scale_edge_colour_manual(values = c("navyblue", "skyblue"),
                                   name = "Kin Type",
                                   aesthetics = "edge_colour",
                                   na.value = "grey50") +
  
  ggraph::geom_node_point(aes(shape = sex, 
                              fill = as.factor(birth_year)),
                          colour = "black",
                          size = 5) +
  
  scale_shape_manual(values = c("F" = 21, "M" = 24, "Unknown" = 23)) +
  
  scale_fill_brewer(palette = "YlOrRd") + 
  
  labs(shape = "Sex") +
  
  ggplot2::theme_void() +
  
  ggplot2::theme(
    panel.grid = element_blank(), 
    axis.text = element_blank(), 
    axis.title = element_blank(), 
    axis.ticks = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 15),
    legend.position = "none",
    plot.margin = unit(rep(1,5), "cm"), 
    margin = margin(20, 0, 40, 10)) +
  
  guides(
    fill = guide_legend(override.aes = list(shape = 21, size = 5, colour = "black", stroke = 0.6)),
    
    edge_colour = guide_legend(override.aes = list(edge_width = 2), order = 2),
    
    shape = guide_legend(order = 3))

print(kin_network)

devEMF::emf("C:/Users/samue/Desktop/kinference_plot.emf", width = 12, height = 7)  
print(kinference_combined)
dev.off()


kinference_combined <- ggpubr::ggarrange(kinf_circle_plot,
                                         kin_network,
                          ncol = 2, 
                          nrow = 1, 
                          labels = c("A", "B"), 
                          legend = "none")

print(kinference_combined)

devEMF::emf("C:/Users/samue/Desktop/kinference_plot.emf", width = 16, height = 8)  
print(kinference_combined)
dev.off()

## Combine CKMR and Kinference results

all_combined <- ggarrange(ckmr_plot, 
                          ckmr_kin_network1,
                          kinf_circle_plot,
                          kin_network,
                                 ncol = 2, 
                                 nrow = 2, 
                                 labels = c("A", "B","C", "D"), 
                                 legend = "none")

print(all_combined)

library(ggplotify)
library(ggpubr)

p1 <- as.ggplot(ckmr_plot)
p2 <- as.ggplot(ckmr_kin_network1)
p3 <- as.ggplot(kinf_circle_plot)
p4 <- as.ggplot(kin_network)

all_combined <- ggarrange(
  p1, p2, p3, p4,
  ncol = 2, nrow = 2,
  font.label = list(size = 18),
  labels = c("A","B","C","D"),
  legend = "none"
)

all_combined

devEMF::emf("C:/Users/samue/Desktop/all_combined_plot.emf", width = 16, height = 16)  
print(all_combined)
dev.off()





ckmr_plot
ckmr_kin_network1
plot
kin_network







