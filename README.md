# License
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses/.

# Manuscript
The primary purpose in distributing this source code is to enable readers to reproduce the numerical results reported in the manuscript "Navigating the Landscape of Fear" [preprint forthcoming] by Marissa Gee, Nicolas Gonzalez-Granda, Sunay Joshi, Nagaprasad Rundrapatna, Anne Somalwar, Stephen P. Ellner, and Alexander Vladimirsky.


# Contributions
- All authors contributed to the formulation of the mathematical model.
- MG, SJ, NR, and AS developed the solution code and visualization. NGG also contributed to the visualization.
- MG, AS, and NGG contributed to data collection and processing.
- MG, SJ, AS, AV, and SPE developed and refined the numerical experiments.
- MG, NGG, SJ, and AS were responsible for early drafting of the manuscript and results.
- MG, AV, and SPE revised and developed the manuscript, and all authors contributed to final editing and approval.


# Instructions
## Requirements
The C++ Code requires one external library: [Boost](http://www.boost.org/), which is used for implementation of multidimensional arrays and heaps.

Currently the C++ code is run through a Makefile that assumes the libraries are installed in the `../lib/` directory.
If those libraries are installed elsewhere, you will have to modify the [Makefile](https://github.com/eikonal-equation/Navigating_LoF/blob/main/Makefile) to ensure the library is properly linked.

The code uses the C++17 standard, and can be compiled using both gcc and icpc.

The results are visualized via [MATLAB](https://www.mathworks.com/products/matlab.html).
The MATLAB code requires the file [arrow.m](https://www.mathworks.com/matlabcentral/fileexchange/278-arrow) available via the MATLAB File Exchange.

### Data
The data for our example based on foraging Samango monkeys can be found in the is primarily based on Figure 1 from ["Living in a landscape of fear: the impact of predation, resource availability and habitat structure on primate range use"](https://doi.org/10.1016/j.anbehav.2013.11.027) by Ben T. Coleman and Russell A. Hill. Our data was generated via the file [ImageRead.m](https://github.com/eikonal-equation/Navigating_LoF/blob/main/visualization/imread/ImageRead.m)  and can be found in the `visualization/imread` directory. 

## Running the Code
Assuming the libraries are appropriately linked, you should be able to compile the code and run the test cases from the manuscript.

The results of the code are stored in a directory labeled `output` that must be created. Additionally, paths are stored in directories `output/RR_PatchSelection`, `output/RR_FoodDepletion` `output/RR_PathDeformation`, and`output/RW_Paths` which must also be created in advance.

When running each example, the format of the inputs is given by

`./main [EXAMPLE NAME] [UTILITY FUNCTION] [STENCIL TYPE] [PATH TRACING]`

#### Possible example names:
`ExampleRiskReward`
`ExampleRiskRewardPT`
`ExampleRealWorld`
`Example9`

#### Possible utility functions
`L` (Linear)
`R` (Square root)
`M` (Sigmoid)
`U` (Wide sigmoid)

#### Stencil type:
`“”` (Default, 5 pt stencil)
`9` (9 point stencil)

Path tracing:
`“”` (Default, compute value functions then perform path tracing)
`P` (Do not compute value functions and read them from file instead, then follow prompts to simulate trajectories. Must have value function data files already generated)

For example, to generate the data needed for Figure 3, run the following four commands
`./main ExampleRiskReward L 9`
`./main ExampleRiskReward R 9`
`./main ExampleRiskReward M 9`
`./main ExampleRiskReward U 9`

Some commands will trigger prompts to specify quantities such as whether to include predator interactions, whether to include food depletion, and how many realized trajectories to generate. When the instructions for generating data below include information beyond the initial command, it is entered via these prompts.

### To generate the data:
All examples were run with the 9 point stencil option.
- Figures 2, 3 -  Run ExampleRiskReward for L, R, M, and U.
- Figure 4 - Run ExampleRiskRewardPT for L, R, and M. Specify one realization, and run each with and without food depletion.
- Figures 5, 6 - Run ExampleRealWorld for L, R, and M. Run with path tracing and specify starting mode 1, no food depletion, allow mode switches, no specified initial energy, 1000 realizations, and many starting locations via prompts.
- Figures 7, 8 - Run ExampleRealWorld for L, R, and M with path tracing. Specify starting mode 1, no food depletion, no mode switches, unspecified initial energy, no specified spotting times. Then run three more times with the same settings but with specified spotting times at 0, 0.5, and 1.0.
- Figures 9, S2 - S4 - Run ExampleMultistage for L, R, and M.
- Figure S1 - Run ExampleRiskRewardPT for L, R, and M. Specify 100 realizations, and run each with and without food depletion.

## Visualizing Output
Code for generating the figures in the manuscript can be found in the [visualization](https://github.com/eikonal-equation/Navigating_LoF/blob/main/visualization) directory. Each file first calls [initialization.m](https://github.com/eikonal-equation/Navigating_LoF/blob/main/visualization/initialization.m), which specifies the example and utility function to be visualized. The code assumes that all data was generated using the 9-point stencil.

### To generate each figure
- Figure 1 - UtilityFunctions.m
- Figure 2 - EnvironmentVisualizer.m, `example = "RiskReward"`
- Figure 3 - PatchSelection.m, `example = "RiskReward"`
- Figure 4 - PathDeformation.m with `depletionComparison = false`, `example = "RiskRewardPT"`
- Figure 5 - EnvironmentVisualizer.m, `example = "RealWorld"`
- Figure 6 - PathFrequency.m and Utilization.m, `example = "RealWorld"`
- Figures 7, 8 - RealWorldWindow.m, `example = "RealWorld"`
- Figure 9 - PathCloud.m, `example = "MultiStage"`
- Figure S1 - DepletionComparison.m, PathDeformation.m with `depletionComparison = true`, `example = "RiskRewardPT"`
- Figure S2 - EnvironmentVisualizer.m, `example = "Multistage"`
- Figure 12 - PathCloud.m, `example = "Multistage"`
