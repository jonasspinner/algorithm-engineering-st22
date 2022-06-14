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
data = read_csv(infile, col_types = "cciiiiiic", guess_max = 2000)
# keys that identify a run
id.vars = c("contender", "iteration", "num_vertices", "num_edges", "num_changed_edges", "input_graph");
# clean data
data = filter(data, result_correct == "yes")
data = filter(data, iteration != 0) # use warm up runs
data = select(data, -result_correct)
# transform from wide to long format
data.melted = get_melted_data(data, normalize_value = FALSE, id.vars = id.vars)

# start output in pdf
pdfname = argv$outfile
ending <- ""
if (!endsWith(pdfname, ".pdf")) ending <- ".pdf"
pdf(paste0("", pdfname, ending), width=10, height=5)


# Print some example plots to the output pdf. These examples are meant as basic 
# plots for your results as well as an inspiration and a starting point for 
# writing your own plots. Feel free to use the methods in utils.r and have a 
# look at the code.

# Plot with num_changed_edges on the x-axis for all runs on a graph with 2^14 vertices. Facetted by graph degree.
g = get_plot_num_changed_edges_as_x_with_fixed_num_vertices(data.melted, 
                                                            title = "Running time", 
                                                            "execution_time_ns",  
                                                            use_log_y = TRUE, 
                                                            y_axis_desc = "time (ns)",
                                                            fixed_num_vertices = 2^14,
                                                            exclude_yes_instances = FALSE,
                                                            column_for_faceting = "degree"
                                                            )
plot(g)

# Same plot as above except not facetted but all lines in one graph.
g = get_plot_num_changed_edges_as_x_with_fixed_num_vertices(data.melted, 
                                                            title = "Running time", 
                                                            "execution_time_ns",  
                                                            use_log_y = TRUE, 
                                                            y_axis_desc = "time (ns)",
                                                            fixed_num_vertices = 2^14,
                                                            exclude_yes_instances = FALSE
)
plot(g)

# Plots with num_changed_edges on the x-axis. One plot for each occurring value of num_vertices.
g = get_plot_num_changed_edges_as_x_facetted_by_num_vertices(data.melted,
                                    title = "Running times",
                                    "execution_time_ns",
                                    use_log_y = TRUE,
                                    y_axis_desc = "time (ns)",
                                    exclude_yes_instances = FALSE)
plot(g)


# Plot with num_changed_edges on the x-axis for all runs on a graph with 2^14 vertices. Applies transformation time_per_edge to each y value (execution_time_ns).
time_per_edge <- function(data) {
  data$value = data$value / data$num_edges
  return (data)
}
g = get_plot_num_changed_edges_as_x_with_fixed_num_vertices(data.melted,
                                                            title = "Time per Edge",
                                                            "execution_time_ns",
                                                            transformation = time_per_edge,
                                                            use_log_y = TRUE,
                                                            y_axis_desc = "time / m (ns)",
                                                            fixed_num_vertices = 2^14,
                                                            exclude_yes_instances = TRUE)
plot(g)


# Plot with num_vertices on the x-axis facetted by num_changed edges. Applies transformation remove_yes_instances to data before plotting, removing all runs with num_changed_edges == 0.
remove_yes_instances <- function(data) {
  data <- subset(data, num_changed_edges != 0)
  return(data)
}
g = get_plot_num_vertices_as_x(data.melted,
                               title = "Running time (all 'no' instances)",
                               "execution_time_ns",
                               transformation = remove_yes_instances,
                               use_log_y = TRUE,
                               y_axis_desc = "time (ns)",
                               column_for_faceting = "num_changed_edges")
plot(g)


# Plot with num_vertices on the x-axis. Applies transformation only_yes_instances to data before plotting, removing all runs with num_changed_edges > 0.
only_yes_instances <- function(data) {
  data <- subset(data, num_changed_edges == 0)
  return(data)
}
g = get_plot_num_vertices_as_x(data.melted, 
                               title = "Running time ('yes' instances)", 
                               "execution_time_ns", 
                               transformation = only_yes_instances, 
                               use_log_y = TRUE, 
                               y_axis_desc = "time (ns)")
plot(g)


# Plot with num_vertices on the x-axis. Applies transformation only_100_changed_edges to data before plotting, removing all runs with num_changed_edges != 100.
only_100_changed_edges <- function(data) {
  data <- subset(data, num_changed_edges == 100)
  return(data)
}
g = get_plot_num_vertices_as_x(data.melted, 
                               title = "Running time (only 'no' instances with 100 changed edges)", 
                               "execution_time_ns",
                               transformation = only_100_changed_edges, 
                               use_log_y = TRUE, 
                               y_axis_desc = "time (ns)")
plot(g)

# Same as above except with average degrees on x-axis.
g = get_plot_degree_as_x(data.melted, 
                         title = "Running time (only 'no' instances with 100 changed edges)", 
                         "execution_time_ns", 
                         transformation = only_100_changed_edges, 
                         y_axis_desc = "time (ns)")
plot(g)
