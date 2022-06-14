library(scales)
library(ggplot2)
library(tidyverse)
library(reshape2)
library(argparser)

parse_arguments <- function() {
  p <- arg_parser('evaluation_parser')
  p <- add_argument(p, "--infile", help="csv inputfile with data to plot")
  p <- add_argument(p, "--outfile", help="name of pdf with plots")
  argv <- parse_args(p)
  return(argv)
}

get_melted_data <- function(data, remove.na = TRUE, normalize_value = TRUE, id.vars) {
  data.melted = pivot_longer(data, cols = !all_of(id.vars), names_to = "variable")
  if (remove.na) {
    data.melted <- filter(data.melted, !is.na(value))
  }
  if(normalize_value) {
    data.melted$value <- data.melted$value / 1000
  }
  return (data.melted)
}

get_plot_num_vertices_as_x <- function(data.melted, title, measurement, transformation, use_log_y = FALSE, column_for_faceting, y_axis_desc, print_min_max = FALSE) {
  data.melted <- data.melted %>% filter(variable == measurement)
  if(!missing(transformation)) {
    data.melted = transformation(data.melted)
  }

  data.melted <- data.melted %>% mutate(degree = (num_edges / num_vertices))
  data.melted$degree <- as.factor(data.melted$degree)
  data.melted$contender <- as.factor(data.melted$contender)
  
  data.group_by <- data.frame()
  if ("num_changed_edges" %in% colnames(data.melted)) {
    data.group_by <- group_by(data.melted, contender, num_vertices, num_edges, degree, num_changed_edges, input_graph)
  } else {
    data.group_by <- group_by(data.melted, contender, num_vertices, num_edges, degree, input_graph)
  }
  data.mean <- summarise(data.group_by, min = min(value), max = max(value), mean = mean(value))
  
  #use this if you want one legend containing all discriminators
    #g <- ggplot(data.mean, aes(x=num_vertices, y=mean, colour = interaction(degree, contender), linetype=interaction(degree, contender), shape=interaction(degree, contender), group = interaction(degree, contender)))
  g <- ggplot(data.mean, aes(x=num_vertices, y=mean, colour = degree, linetype=contender, shape=contender, group = interaction(degree, contender)))
  g <- g + geom_point()
  if(print_min_max) {
    g <- g + geom_pointrange(mapping = aes(ymin = min, ymax = max))
  }
  g <- g + geom_line()
  g <- g + theme_bw()
  g <- g + theme(axis.text.x = element_text(angle = 0), strip.background = element_rect(fill = "white", colour = "white"))
  if(!missing(column_for_faceting)) {
  g <- g + facet_grid(~ eval(as.name(column_for_faceting)))
  }
  g <- g + labs(title = title, y = y_axis_desc, x = "num vertices")
  legend.title = "contender, degree"
  contender.num = nlevels(data.mean$contender)
  degree.num = nlevels(data.mean$degree)

  #use this if you want one legend containing all discriminators
    #g <- g + scale_linetype_manual(legend.title,
    #                             values = rep(1:contender.num, each = degree.num))
    #g <- g + scale_shape_manual(legend.title,
    #                             values = 15 + rep(1:contender.num, each = degree.num))

    #g <- g + scale_color_manual(legend.title,
    #                          values = rep(1:degree.num, contender.num))
  if(use_log_y) {
    g <- g + scale_y_continuous(trans = log10_trans())
  }
    #breaks = trans_breaks("log2", function(x) 2^x),
  labels = trans_format("log2", math_format(2^.x))
  g = g + scale_x_continuous(trans = "log2", labels = trans_format("log2", math_format(2^.x)))
  return (g)
}

get_plot_degree_as_x <- function(data.melted, title, measurement, transformation, use_log_y = FALSE, column_for_faceting, y_axis_desc, print_min_max = FALSE) {
  data.melted <- data.melted %>% filter(variable == measurement)
  if(!missing(transformation)) {
    data.melted = transformation(data.melted)
  }

  data.melted <- data.melted %>% mutate(degree = (num_edges / num_vertices))
  data.melted$degree <- as.factor(data.melted$degree)
  data.melted$contender <- as.factor(data.melted$contender)
  data.melted$num_vertices <- as.factor(data.melted$num_vertices)
  data.group_by <- group_by(data.melted, contender, num_vertices, num_edges, degree, input_graph)
  data.mean <- summarise(data.group_by, min = min(value), max = max(value), mean = mean(value))
  
  g <- ggplot(data.mean, aes(x=degree, y=mean, colour = num_vertices, linetype=contender, shape=contender, group = interaction(num_vertices, contender)))
  g <- g + geom_point()
  if(print_min_max) {
    g <- g + geom_pointrange(mapping = aes(ymin = min, ymax = max))
  }
  g <- g + geom_line()
  g <- g + theme_bw()
  g <- g + theme(axis.text.x = element_text(angle = 0), strip.background = element_rect(fill = "white", colour = "white"))
  if(!missing(column_for_faceting)) {
  g <- g + facet_grid(~ eval(as.name(column_for_faceting)))
  }
  g <- g + labs(title = title, y = y_axis_desc, x = "avg. degree")
  legend.title = "contender, degree"
  if(use_log_y) {
    g <- g + scale_y_continuous(trans = log10_trans())
  }
    #breaks = trans_breaks("log2", function(x) 2^x),
  labels = trans_format("log2", math_format(2^.x))
  #g = g + scale_x_continuous(trans = "log2", labels = trans_format("log2", math_format(2^.x)))
  return (g)
}


get_plot_num_changed_edges_as_x_facetted_by_num_vertices <- function(data.melted, 
                                            title, 
                                            measurement, 
                                            transformation, 
                                            use_log_y = FALSE,
                                            y_axis_desc, 
                                            print_min_max = FALSE,
                                            exclude_yes_instances = FALSE) {
  data.melted <- data.melted %>% filter(variable == measurement)
  
  if (exclude_yes_instances) {
    data.melted <- data.melted %>% filter(num_changed_edges != 0)
  }
  
  log10_num_changed_edges <- data.melted[,"num_changed_edges"]
  log10_num_changed_edges <- log10_num_changed_edges[log10_num_changed_edges != 0]
  log10_num_changed_edges <- unique(log10(log10_num_changed_edges))
  # breaks <- 10^(log10_num_changed_edges)
  labels <- paste0("10^", log10_num_changed_edges)
  if (!exclude_yes_instances && 0 %in% data.melted$num_changed_edges) {
    # breaks <- c(0, breaks)
    labels <- c("0", labels)
  }
  
  
  if(!missing(transformation)) {
    data.melted = transformation(data.melted)
  }
  
  data.melted <- data.melted %>% mutate(degree = (num_edges / num_vertices))
  data.melted$degree <- as.factor(data.melted$degree)
  data.melted$contender <- as.factor(data.melted$contender)
  data.melted$num_changed_edges <- as.factor(data.melted$num_changed_edges)
  data.group_by <- group_by(data.melted, contender, num_vertices, num_edges, degree, num_changed_edges, input_graph)
  data.mean <- summarise(data.group_by, min = min(value), max = max(value), mean = mean(value))
  
  g <- ggplot(data.mean, aes(x=num_changed_edges, y=mean, colour = degree, linetype=contender, shape=contender, group = interaction(degree, contender)))
  g <- g + geom_point()
  if(print_min_max) {
    g <- g + geom_pointrange(mapping = aes(ymin = min, ymax = max))
  }
  g <- g + geom_line()
  g <- g + theme_bw()
  g <- g + theme(axis.text.x = element_text(angle = 0), strip.background = element_rect(fill = "white", colour = "white"))
  
  num_vertices.vals <- unique(data.mean$num_vertices)
  num_vertices.labels <- sapply(num_vertices.vals, function(val) paste0("n = ", toString(val)))
  names(num_vertices.labels) <- sapply(num_vertices.vals, function(val) toString(val))
  g <- g + facet_grid(~ num_vertices, labeller = labeller(num_vertices = num_vertices.labels))

  g <- g + labs(title = title, y = y_axis_desc, x = "num changed edges")
  legend.title = "contender, degree"
  contender.num = nlevels(data.mean$contender)
  degree.num = nlevels(data.mean$degree)
  
  if(use_log_y) {
    g <- g + scale_y_continuous(trans = log10_trans())
  }
  
  g = g + scale_x_discrete(labels=labels)
  return (g)
  
}

get_plot_num_changed_edges_as_x_with_fixed_num_vertices <- function(data.melted, 
                                                                    title, 
                                                                    measurement, 
                                                                    transformation, 
                                                                    fixed_num_vertices, 
                                                                    use_log_y = FALSE, 
                                                                    column_for_faceting, 
                                                                    y_axis_desc, 
                                                                    print_min_max = FALSE,
                                                                    exclude_yes_instances = FALSE) {
  data.melted <- data.melted %>% filter(variable == measurement)
  data.melted <- data.melted %>% filter(num_vertices == fixed_num_vertices)
  
  if (exclude_yes_instances) {
    data.melted <- data.melted %>% filter(num_changed_edges != 0)
  }
  
  log10_num_changed_edges <- data.melted[,"num_changed_edges"]
  log10_num_changed_edges <- log10_num_changed_edges[log10_num_changed_edges != 0]
  log10_num_changed_edges <- unique(log10(log10_num_changed_edges))
  # breaks <- 10^(log10_num_changed_edges)
  labels <- paste0("10^", log10_num_changed_edges)
  if (!exclude_yes_instances && 0 %in% data.melted$num_changed_edges) {
    # breaks <- c(0, breaks)
    labels <- c("0", labels)
  }
  
  
  
  title <- paste0(title, " (n = ", fixed_num_vertices, ")")
  if(!missing(transformation)) {
    data.melted = transformation(data.melted)
  }
  
  data.melted <- data.melted %>% mutate(degree = (num_edges / num_vertices))
  data.melted$degree <- as.factor(data.melted$degree)
  data.melted$contender <- as.factor(data.melted$contender)
  data.melted$num_changed_edges <- as.factor(data.melted$num_changed_edges)
  data.group_by <- group_by(data.melted, contender, num_vertices, num_edges, degree, num_changed_edges, input_graph)
  data.mean <- summarise(data.group_by, min = min(value), max = max(value), mean = mean(value))
  
  g <- ggplot(data.mean, aes(x=num_changed_edges, y=mean, colour = degree, linetype=contender, shape=contender, group = interaction(degree, contender)))
  g <- g + geom_point()
  if(print_min_max) {
    g <- g + geom_pointrange(mapping = aes(ymin = min, ymax = max))
  }
  g <- g + geom_line()
  g <- g + theme_bw()
  g <- g + theme(axis.text.x = element_text(angle = 0), strip.background = element_rect(fill = "white", colour = "white"))
  if(!missing(column_for_faceting)) {
    g <- g + facet_grid(~ eval(as.name(column_for_faceting)))
  }
  g <- g + labs(title = title, y = y_axis_desc, x = "num changed edges")
  legend.title = "contender, degree"
  contender.num = nlevels(data.mean$contender)
  degree.num = nlevels(data.mean$degree)
  
  if(use_log_y) {
    g <- g + scale_y_continuous(trans = log10_trans())
  }
  
  g = g + scale_x_discrete(labels=labels)
  return (g)
  
}
