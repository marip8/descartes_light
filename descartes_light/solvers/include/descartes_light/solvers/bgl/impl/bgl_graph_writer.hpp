#pragma once

#include <descartes_light/descartes_macros.h>
DESCARTES_IGNORE_WARNINGS_PUSH
#include <fstream>
DESCARTES_IGNORE_WARNINGS_POP

#include <descartes_light/solvers/bgl/boost_graph_types.h>

const static std::string FILLCOLOR_ATTR = "fillcolor";
const static std::string STYLE_ATTR = "style";
const static std::string LABEL_ATTR = "label";
const static std::string RANKDIR_ATTR = "rankdir";

namespace descartes_light
{
// Forward declare helper function from solver
template <typename FloatType>
VertexDesc<FloatType> getLowestCostLastVertex(const SubGraph<FloatType>&);

/**
 * @brief Adds properties to a sub-graph and its children for display in a .dot file
 */
template <typename FloatType>
inline void decorateSubGraph(SubGraph<FloatType>& g)
{
  // Set the top-level graph properties
  boost::get_property(g, boost::graph_name) = "G";  
  boost::get_property(g, boost::graph_vertex_attribute)[STYLE_ATTR] = "filled";
  boost::get_property(g, boost::graph_graph_attribute)[RANKDIR_ATTR] = "LR";

  // Add graph names to the children of the main graph
  {
    typename SubGraph<FloatType>::children_iterator start, end;
    boost::tie(start, end) = g.children();
    for (auto it = start; it != end; ++it)
    {
      std::stringstream ss;
      ss << "cluster_" << std::distance(start, it);
      boost::get_property(*it, boost::graph_name) = ss.str();
    }
  }

  // Add name and color properties to the graph vertices
  VertexIt<FloatType> start, end;
  boost::tie(start, end) = boost::vertices(g);
  for (auto it = start; it != end; ++it)
  {
    std::stringstream name;
    name << std::setprecision(4) << "v" << *it << ": " << g[*it].distance;
    boost::get(boost::vertex_attribute, g)[*it][LABEL_ATTR] = name.str();

    // Add colors
    auto& color_prop = boost::get(boost::vertex_attribute, g)[*it][FILLCOLOR_ATTR];
    switch (g[*it].color)
    {
      case boost::default_color_type::white_color:
        color_prop = "white";
        break;
      case boost::default_color_type::gray_color:
        color_prop = "gray";
        break;
      case boost::default_color_type::black_color:
        color_prop = "gray30";
        break;
      case boost::default_color_type::red_color:
        color_prop = "red";
        break;
      default:
        break;
    }
  }

  // Add properties to the edges of the subgraph
  auto weights = boost::get(boost::edge_weight, g);

  EdgeIt<FloatType> first, last;
  boost::tie(first, last) = boost::edges(g);
  for (auto it = first; it != last; ++it)
  {
    std::stringstream ss;
    ss << std::setprecision(4) << weights[*it];
    boost::get(boost::edge_attribute, g)[*it][LABEL_ATTR] = ss.str();
  }
}

template <typename FloatType>
void BGLLadderGraphSolver<FloatType>::writeGraph(const std::string& filename)
{
  std::ofstream file(filename);
  if (!file.good())
    throw std::runtime_error("Failed to open file '" + filename + "'");

  decorateSubGraph(graph_);
  boost::write_graphviz(file, graph_);
}

template <typename FloatType>
void BGLLadderGraphSolver<FloatType>::writeGraphWithPath(const std::string& filename)
{
  std::ofstream file(filename);
  if (!file.good())
    throw std::runtime_error("Failed to open file '" + filename + "'");

  // Get the path through the graph
  const std::vector<VertexDesc<FloatType>> path = reconstructPath(source_, getLowestCostLastVertex(graph_));

  // Decorate the graph with properties
  decorateSubGraph(graph_);

  // Colorize the path
  for (const VertexDesc<FloatType>& v : path)
  {
    boost::get(boost::vertex_attribute, graph_)[v][FILLCOLOR_ATTR] = "green";
  }

  boost::write_graphviz(file, graph_);
}

}  // namespace descartes_light
