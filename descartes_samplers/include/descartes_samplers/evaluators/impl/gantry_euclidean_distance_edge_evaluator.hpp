/*
 * Software License Agreement (Apache License)
 *
 * Copyright (c) 2016, Southwest Research Institute
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
#ifndef DESCARTES_SAMPLERS_EVALUATORS_GANTRY_EUCLIDEAN_DISTANCE_EDGE_EVALUATOR_HPP
#define DESCARTES_SAMPLERS_EVALUATORS_GANTRY_EUCLIDEAN_DISTANCE_EDGE_EVALUATOR_HPP

#include <descartes_light/descartes_macros.h>
DESCARTES_IGNORE_WARNINGS_PUSH
#include <cmath>
DESCARTES_IGNORE_WARNINGS_POP

#include <descartes_samplers/evaluators/gantry_euclidean_distance_edge_evaluator.h>

namespace descartes_light
{
template <typename FloatType>
GantryEuclideanDistanceEdgeEvaluator<FloatType>::GantryEuclideanDistanceEdgeEvaluator(int dof)
  : EdgeEvaluator<FloatType>(static_cast<std::size_t>(dof))
{
}

template <typename FloatType>
std::pair<bool, FloatType> GantryEuclideanDistanceEdgeEvaluator<FloatType>::considerEdge(const FloatType* start,
                                                                                         const FloatType* end)
{
  FloatType cost = 0.0;
  for (std::size_t i = 2; i < this->dof_; ++i)
    cost += std::pow(end[i] - start[i], FloatType(2));

  return std::make_pair(true, cost);
}

}  // namespace descartes_light

#endif  // DESCARTES_SAMPLERS_EVALUATORS_GANTRY_EUCLIDEAN_DISTANCE_EDGE_EVALUATOR_HPP
