library(scales)
library(ggplot2)
library(tidyverse)
library(reshape2)
library(argparser)

source("utils.r")

options(tibble.width = Inf)



argv = parse_arguments()
infile = argv$infile
outfile = argv$outfile

# read data
data = read_csv(infile, col_types = cols( p = col_factor()), guess_max = 2000)
# keys that identify a run
id.vars = c("contender", "iteration", "num_vertices", "num_edges", "input_graph");
# clean data
data = filter(data, result_correct == "yes")
data = filter(data, iteration != 0) # use a warm up run
data = select(data, -result_correct)
# transform from wide to long format
data.melted = get_melted_data(data, normalize_value = FALSE, id.vars = id.vars)



# start output in pdf
pdfname = argv$outfile
pdf(paste("", pdfname, ".pdf",sep=""), width=10, height=5)
g = get_plot_num_vertices_as_x(data.melted, title = "Running time", "execution_time_ns", use_log_y = TRUE, y_axis_desc = "time (ns)")
plot(g)

time_per_edge <- function(data) {
  data$value = data$value / data$num_edges
  return (data)
}

g = get_plot_num_vertices_as_x(data.melted, title = "Time per Edge", "execution_time_ns", transformation = time_per_edge, y_axis_desc = "time / m (ns)")
plot(g)
g = get_plot_degree_as_y(data.melted, title = "Time per Edge", "execution_time_ns", transformation = time_per_edge, y_axis_desc = "time / m (ns)")
plot(g)
g = get_plot_degree_as_y(data.melted, title = "Time per Edge with min max", "execution_time_ns", transformation = time_per_edge, y_axis_desc = "time / m (ns)", print_min_max = TRUE)
plot(g)
