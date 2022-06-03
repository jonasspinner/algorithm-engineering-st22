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
  data.group_by <- group_by(data.melted, contender, num_vertices, num_edges, degree, input_graph)
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
  g <- g + facet_grid(~ column_for_faceting)
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

get_plot_degree_as_y <- function(data.melted, title, measurement, transformation, use_log_y = FALSE, column_for_faceting, y_axis_desc, print_min_max = FALSE) {
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
  g <- g + facet_grid(~ column_for_faceting)
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
