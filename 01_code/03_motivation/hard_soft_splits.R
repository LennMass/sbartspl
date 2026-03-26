# ============================================================
#  Replication of Figure 3 in Linero (2022)
#  "SoftBart: Soft Bayesian Additive Regression Trees"
#  arXiv:2210.16375
#
#  Hard vs. Soft decision tree on [0,1]^2
#  - Same 4 colours in both panels
#  - Soft panel: per-leaf gradient white (phi~0) -> full colour (phi~1)
#
#  Required: ggplot2, dplyr, patchwork
# ============================================================

library(ggplot2)
library(dplyr)
library(patchwork)

# ── Logistic soft split ───────────────────────────────────────
psi <- function(x, C, tau) 1 / (1 + exp(-(x - C) / tau))

# ── Tree structure ────────────────────────────────────────────
# b1: split x1 at C1=0.7  (root)
# b2: split x2 at C2=0.4  (left  child of b1)
# b3: split x2 at C3=0.6  (right child of b1)
# Leaves:  mu1 (LL), mu2 (LR), mu3 (RL), mu4 (RR)
C1  <- 0.7
C2  <- 0.4
C3  <- 0.6
tau <- 0.04    # bandwidth — decrease for sharper, increase for smoother

# ── Leaf colours (same in BOTH panels) ───────────────────────
col <- c("1" = "#F8766D",   # mu_t1  red/salmon
				 "2" = "#7CAE00",   # mu_t2  green
				 "3" = "#00BFC4",   # mu_t3  teal
				 "4" = "#C77CFF")   # mu_t4  purple

# Blend hex colour with white; w=0 -> white, w=1 -> full colour
blend_white <- function(hex, w) {
	r0 <- col2rgb(hex)
	r  <- round((1 - w) * 255 + w * r0[1])
	g  <- round((1 - w) * 255 + w * r0[2])
	b  <- round((1 - w) * 255 + w * r0[3])
	rgb(r, g, b, maxColorValue = 255)
}

# ── Grid ──────────────────────────────────────────────────────
n <- 300
grid <- expand.grid(x1 = seq(0, 1, length.out = n),
										x2 = seq(0, 1, length.out = n))

# ── HARD leaf assignment ──────────────────────────────────────
hard <- grid %>%
	mutate(
		leaf = case_when(
			x1 <= C1 & x2 <= C2 ~ "1",
			x1 <= C1 & x2 >  C2 ~ "2",
			x1 >  C1 & x2 <= C3 ~ "3",
			TRUE                 ~ "4"
		),
		hex = col[leaf]        # same colour, full saturation everywhere
	)

# ── SOFT leaf membership probabilities ────────────────────────
soft <- grid %>%
	mutate(
		L1   = psi(x1, C1, tau),  R1 = 1 - L1,
		L2   = psi(x2, C2, tau),  R2 = 1 - L2,
		L3   = psi(x2, C3, tau),  R3 = 1 - L3,
		phi1 = L1 * L2,
		phi2 = L1 * R2,
		phi3 = R1 * L3,
		phi4 = R1 * R3,
		# dominant leaf
		leaf = case_when(
			phi1 >= pmax(phi2, phi3, phi4) ~ "1",
			phi2 >= pmax(phi1, phi3, phi4) ~ "2",
			phi3 >= pmax(phi1, phi2, phi4) ~ "3",
			TRUE                           ~ "4"
		),
		# phi of the dominant leaf
		phi_dom = case_when(
			leaf == "1" ~ phi1,
			leaf == "2" ~ phi2,
			leaf == "3" ~ phi3,
			TRUE        ~ phi4
		),
		# rescale phi_dom to a blend weight in [0,1]
		# phi_dom ~ 0.25 (near triple boundary) -> white
		# phi_dom ~ 1.0  (deep in leaf)          -> full colour
		w = pmin(1, pmax(0, (phi_dom - 0.25) / 0.75)),
		# per-point blended colour
		hex = mapply(function(lf, wt) blend_white(col[lf], wt),
								 leaf, w)
	)

# ── Shared elements ───────────────────────────────────────────
base_theme <- theme_bw(base_size = 13) +
	theme(panel.grid.minor = element_blank(),
				plot.title       = element_text(hjust = 0.5, face = "bold"),
				legend.position  = "bottom",
				axis.title       = element_text(size = 12))

leaf_legend <- c("1" = col["1"], "2" = col["2"],
								 "3" = col["3"], "4" = col["4"])

label_df <- data.frame(
	x   = c(0.35, 0.35, 0.85, 0.85),
	y   = c(0.18, 0.70, 0.28, 0.82),
	lbl = c("mu[t1]","mu[t2]","mu[t3]","mu[t4]")
)

# ── HARD plot ─────────────────────────────────────────────────
p_hard <- ggplot(hard, aes(x1, x2)) +
	geom_point(aes(colour = leaf), size = 0.3, shape = 15) +
	scale_colour_manual(
		values = leaf_legend,
		labels = c(expression(mu[t1]), expression(mu[t2]),
							 expression(mu[t3]), expression(mu[t4])),
		name   = NULL
	) +
	# branch lines
	geom_vline(xintercept = C1, linewidth = 0.7) +
	annotate("segment",
					 x=0, xend=C1, y=C2, yend=C2, linewidth=0.7) +
	annotate("segment",
					 x=C1, xend=1, y=C3, yend=C3, linewidth=0.7) +
	# leaf labels
	annotate("text", x=label_df$x, y=label_df$y,
					 label=label_df$lbl, parse=TRUE,
					 size=4.5, fontface="bold") +
	scale_x_continuous(breaks=c(0,.25,.5,.75,1),
										 limits=c(0,1), expand=c(.005,.005)) +
	scale_y_continuous(breaks=c(0,.25,.5,.75,1),
										 limits=c(0,1), expand=c(.005,.005)) +
	labs(x=expression(x[1]), y=expression(x[2]),
			 title="Hard decision tree") +
	base_theme +
	guides(colour = guide_legend(
		override.aes = list(size=4, shape=15), nrow=1))

# ── SOFT plot ─────────────────────────────────────────────────
# We use scale_colour_identity so each point gets its exact
# blended hex colour; legend overridden with the pure leaf colours
p_soft <- ggplot(soft, aes(x1, x2)) +
	geom_point(aes(colour = hex), size = 0.3, shape = 15) +
	scale_colour_identity(
		guide  = "legend",
		breaks = unname(col),                  # pure colours for legend
		labels = c(expression(mu[t1]), expression(mu[t2]),
							 expression(mu[t3]), expression(mu[t4])),
		name   = NULL
	) +
	# dashed boundary lines
	geom_vline(xintercept = C1, linewidth=0.7, linetype="dashed") +
	annotate("segment",
					 x=0, xend=C1, y=C2, yend=C2,
					 linewidth=0.7, linetype="dashed") +
	annotate("segment",
					 x=C1, xend=1, y=C3, yend=C3,
					 linewidth=0.7, linetype="dashed") +
	# leaf labels
	annotate("text", x=label_df$x, y=label_df$y,
					 label=label_df$lbl, parse=TRUE,
					 size=4.5, fontface="bold") +
	scale_x_continuous(breaks=c(0,.25,.5,.75,1),
										 limits=c(0,1), expand=c(.005,.005)) +
	scale_y_continuous(breaks=c(0,.25,.5,.75,1),
										 limits=c(0,1), expand=c(.005,.005)) +
	labs(x=expression(x[1]), y=expression(x[2]),
			 title="Soft decision tree") +
	base_theme +
	guides(colour = guide_legend(
		override.aes = list(size=4, shape=15,
												colour=unname(col)),
		nrow=1))

# ── Combine & save ────────────────────────────────────────────
fig3 <- p_hard + p_soft +
	plot_layout(guides = "collect") &
	theme(legend.position = "bottom")

ggsave("figure3_softbart.pdf",
			 plot=fig3, width=10, height=5.2, device="pdf")

ggsave("figure3_softbart.png",
			 plot=fig3, width=10, height=5.2, dpi=300)

message("Saved: figure3_softbart.pdf  and  figure3_softbart.png")