# Pre-amble ---------------------------------------------------------------

# Clearing Memory
rm(list = ls())

# Loading packages
require(magrittr)
require(tidyverse)
require(cowplot)


#### MANUAL CHANGE REQUIRED ####
# Load models
load("./data/spec_curve_models.RData")


#### MANUAL CHANGE REQUIRED ####
# Names of treat vars
# Important they are in the same order as the Models columns
treatvars <- c("treat_1", "treat_2")

col_location <- which(names(Models) %in% treatvars)

curves <- cbind(treatvars, col_location) %>% as.data.frame()
curves$col_location %<>% as.character() %>%  as.numeric()


#### MANUAL CHANGE REQUIRED ####
# Direction of curves
# 1 for ascending, -1 for descending
curves$direct <- c(1)

#### MANUAL CHANGE REQUIRED ####
# Order to print the options in
printOrder <- c("Model_Type", "Additional_Controls", 
                "Y_Var_Form", "Exclude_Rec", 
                "Exclude_Long_Treat")

# Saves number of models that are significant and the qunatiles
perSig <- c()
q05 <- c()
q95 <- c()

# Loop over different treat vars
# for(curve_index in 1:nrow(curves)){
  
  # Testing
  curve_index <- 1
  # 
  Models$estimate = Models[ , (curves[curve_index, 2])]
  Models$ub = Models[ , (curves[curve_index, 2] + 1)]
  Models$lb = Models[ , (curves[curve_index, 2] + 2)]
  
  # For testing or subsample
  l <- which(Models$estimate == 0)
  
  if(length(l) > 0){
    Models <- Models[-l,]
    print("Warning, some models may not have run!")
  }
  
  # Getting Significance
  Models$sig <- (sign(Models$ub) == sign(Models$lb)) %>% as.factor()
  
  # Ordering by estimate
  Models <- Models[order(Models$estimate * curves$direct[curve_index]),]
  Models$Order <- 1:nrow(Models)
  
  
  #### MANUAL CHANGE REQUIRED ####
  # Get preferred specification
  pref <- which(Models$Model_Type == "Triple Diff" &
                  Models$Additional_Controls == "No" &
                  Models$Y_Var_Form == "Log" &
                  Models$Exclude_Rec == "Yes" &
                  Models$Exclude_Long_Treat == "Yes")
  
  if(length(pref) == 0) print("Preferred Spec not found!")
  
  # Plotting ----------------------------------------------------------------
  # Curve (shouldn't need to change any of this)
  
  # Makes sure that all ticks aren't gray if all significant
  if(length(unique(Models$sig)) == 1){
    gstart = 0
  } else {
    gstart = .7
  }
  
  # Makes top curve
  curve <- ggplot(data = Models) +
    geom_point(
      mapping = aes(x = Order, y = estimate, color = sig)
    ) +
    scale_color_grey(start = gstart,
                     end = 0) +
    geom_linerange(
      mapping =  aes(x = Order, ymin = lb, ymax = ub), colour = "blue", size = .2, alpha = .8
    ) +  
    geom_vline(xintercept = pref, color = "red", linetype = "dashed") + 
    theme(legend.position = "none") +
    labs(x = "Regression Number", y = "Estimate")
  
  #### MANUAL CHANGE REQUIRED ####
  # Specs (need to edit this based on column names and order)
  Models %>% 
    gather(key, value, Model_Type:Exclude_Long_Treat) -> plotDat
  
  # Makes bottom plot
  specs <- ggplot(data = plotDat, 
                  aes(x = plotDat$Order,
                      y = plotDat$value,
                      color = plotDat$sig)) + 
    scale_color_grey(start = gstart,
                     end = 0) +
    geom_point(size = 4,
               shape = 124) +
    facet_grid(rows = vars(key), scales = "free", space = "free") + 
    theme(
      axis.line = element_line("black", size = .5),
      legend.position = "none",
      panel.spacing = unit(.75, "lines"),
      axis.text = element_text(colour = "black"),
      strip.text.x = element_blank(),
      strip.text.y = element_text(face = "bold",
                                  size = 8),
      strip.background.y = element_blank()) +
    labs(x = "", y = "")
  
  #### MANUAL CHANGE REQUIRED ####
  # Fixing height for vars (lines in gp$heights that have null in the name are the ones to change)
  gp <- ggplotGrob(specs)
  gp$heights[7] <- gp$heights[7] * 5
  gp$heights[9] <- gp$heights[9] * 5
  gp$heights[11] <- gp$heights[11] * 5
  gp$heights[13] <- gp$heights[13] * 7
  gp$heights[15] <- gp$heights[15] * 8
  
  plot_grid(curve,
            gp,
            labels = c(),
            align = "v",
            axis = "rbl",
            rel_heights = c(2,6),
            ncol = 1)
  
  
  #### MANUAL CHANGE REQUIRED ####
  # This is where it gets saved to
  savename = paste0("../../../../Apps/Overleaf/SBA Household Collateral Bunching/figures/spec_curve_", 
                    curves[curve_index, 1], ".png")
  
  
  #### MANUAL CHANGE REQUIRED ####
  # You'll need to play around with the height and width
  ggsave(filename = savename,
         width = 10,
         height = 10, 
         units = "in")
  
  # This is saving some info about the curves (like percent significant, etc.)
  nrow(Models)
  q05[curve_index] <- quantile(Models$estimate, .05)
  q95[curve_index] <- quantile(Models$estimate, .95)
  
  m <- which(Models$ub < 0) %>% length()
  n <- which(Models$lb > 0) %>% length()
  l <- max(m, n)
  perSig[curve_index] <- l/nrow(Models)
  
}


perSig
q05
q95


