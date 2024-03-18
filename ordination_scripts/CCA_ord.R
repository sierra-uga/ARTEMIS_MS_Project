CCA_ord <- function(top_CCA_ref, actual_CCA_ref, title, name_of_graph){
  arrowmat <- scores(top_CCA_ref, display = "species")
  
  # Add labels, make a data.frame
  arrowdf <- data.frame(labels = rownames(arrowmat), arrowmat)
  
  # Define the arrow aesthetic mapping
  arrow_map = aes(xend = CCA1, yend = CCA2, x = 0, y = 0, shape = NULL, color = NULL, 
                  label = labels)
  label_map = aes(x = 1.2 * CCA1, y = 1.2 * CCA2, shape = NULL, color = NULL, 
                  label = labels)
  # Make a new graphic
  arrowhead = arrow(length = unit(0.02, "npc"))
  name_of_graph = actual_CCA_ref + ggtitle(title) + geom_segment(arrow_map, size = 0.5, data = arrowdf, color = "gray", 
                                                                 arrow = arrowhead) + geom_text(label_map, size = 3, data = arrowdf) + remove_grid + color_point
  name_of_graph
}