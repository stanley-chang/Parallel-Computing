#include <iostream>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <getopt.h>

#include "a-star-navigator.h"
#include "VideoOutput.h"
#include "Utility.h"
#include <omp.h>

void simulate_waves(ProblemData &problemData, int t) {
    auto &islandMap = problemData.islandMap;
    float (&secondLastWaveIntensity)[MAP_SIZE][MAP_SIZE] = *problemData.secondLastWaveIntensity;
    float (&lastWaveIntensity)[MAP_SIZE][MAP_SIZE] = *problemData.lastWaveIntensity;
    float (&currentWaveIntensity)[MAP_SIZE][MAP_SIZE] = *problemData.currentWaveIntensity;
#pragma omp parallel for num_threads(4) schedule(static, 256)
    for (int x = 1; x < MAP_SIZE - 1; ++x) {
        for (int y = 1; y < MAP_SIZE - 1; ++y) {
            currentWaveIntensity[x][y] = islandMap[x][y] >= LAND_THRESHOLD ? 0.0f : std::clamp(lastWaveIntensity[x][y] + (lastWaveIntensity[x][y] - secondLastWaveIntensity[x][y] + ATTACK_FACTOR*(lastWaveIntensity[x][y - 1]
                                                                                                                                                                                                   + lastWaveIntensity[x - 1][y]
                                                                                                                                                                                                   + lastWaveIntensity[x + 1][y]
                                                                                                                                                                                                   + lastWaveIntensity[x][y + 1]
                                                                                                                                                                                                   - 4 * lastWaveIntensity[x][y])) * std::clamp(
                    ENERGY_PRESERVATION_FACTOR * (LAND_THRESHOLD - 0.1f * islandMap[x][y]),0.0f, 1.0f),0.0f,1.0f);


            /*
             * TODO@Students: This is the first part of the main calculation. Here, we simulate the waves on the
             * stormy seas. Look at all TODOS for parallelization.
             */

//            // Simulate some waves
//
//            // The acceleration is the relative difference between the current point and the last.
//            float acceleration = lastWaveIntensity[x][y - 1]
//                                 + lastWaveIntensity[x - 1][y]
//                                 + lastWaveIntensity[x + 1][y]
//                                 + lastWaveIntensity[x][y + 1]
//                                 - 4 * lastWaveIntensity[x][y];
//
//            // The acceleration is multiplied with an attack value, specifying how fast the system can accelerate.
//            acceleration *= ATTACK_FACTOR;
//
//            // The last_velocity is calculated from the difference between the last intensity and the
//            // second to last intensity
//            float last_velocity = lastWaveIntensity[x][y] - secondLastWaveIntensity[x][y];
//
//            // energy preserved takes into account that storms lose energy to their environments over time. The
//            // ratio of energy preserved is higher on open water, lower close to the shore and 0 on land.
//            float energyPreserved = std::clamp(
//                    ENERGY_PRESERVATION_FACTOR * (LAND_THRESHOLD - 0.1f * islandMap[x][y]),
//                    0.0f, 1.0f);
//
//            // There aren't any waves on land.
//            if (islandMap[x][y] >= LAND_THRESHOLD) {
//                currentWaveIntensity[x][y] = 0.0f;
//            } else {
//                currentWaveIntensity[x][y] =
//                        std::clamp(lastWaveIntensity[x][y] + (last_velocity + acceleration) * energyPreserved,
//                                   0.0f,
//                                   1.0f);
//            }
        }
    }
}

/*
 * Since all pirates like navigating by the stars, Captain Jack's favorite pathfinding algorithm is called A*.
 * Unfortunately, sometimes you just have to make do with what you have. So here we use a search algorithm that searches
 * the entire domain every time step and calculates all possible ship positions.
 */
bool findPathWithExhaustiveSearch(ProblemData &problemData, int timestep,
                                  std::vector<Position2D> &pathOutput) {
    auto &start = problemData.shipOrigin;
    auto &portRoyal = problemData.portRoyal;
    auto &islandMap = problemData.islandMap;
    auto &currentWaveIntensity = *problemData.currentWaveIntensity;
    auto &lastWaveIntensity = *problemData.lastWaveIntensity;

    /*
     * TODO@Students: This is the second part of the main calculation. Here, we find the shortest path to get back to
     * Port Royal without our ship capsizing. Take a look at whether this should be parallelized.
     */



    // The Jolly Mon (Jack's ship) is not very versatile. It can only move along the four cardinal directions by one
    // square each and along their diagonals. Alternatively, it can just try to stay where it is.
    Position2D neighbors[] = {
            Position2D(-1, 0),
            Position2D(0, -1),
            Position2D(1, 0),
            Position2D(0, 1),
            Position2D(-1, -1),
            Position2D(-1, 1),
            Position2D(1, -1),
            Position2D(1, 1),
            Position2D(0, 0)
    };

//    std::cerr << "Searching for position: " << portRoyal.x << " " << portRoyal.y << std::endl;
    int numPossiblePositions = 0;

    bool (&currentShipPositions)[MAP_SIZE][MAP_SIZE] = *problemData.currentShipPositions;
    bool (&previousShipPositions)[MAP_SIZE][MAP_SIZE] = *problemData.previousShipPositions;

    // We could always have been at the start in the previous frame since we get to choose when we start our journey.
    previousShipPositions[start.x][start.y] = true;
    #pragma omp parallel for num_threads(4) schedule(static, 256)
        // Ensure that our new buffer is set to zero. We need to ensure this because we are reusing previously used buffers.
        for (int x = 0; x < MAP_SIZE; ++x) {
            for (int y = 0; y < MAP_SIZE; ++y) {
                currentShipPositions[x][y] = false;
            }
        }

    // Do the actual path finding.
    for (int x = 0; x < MAP_SIZE; ++x) {
        for (int y = 0; y < MAP_SIZE; ++y) {
            // If there is no possibility to reach this position then we don't need to process it.
            if (!previousShipPositions[x][y]) {
                continue;
            }
            Position2D previousPosition(x, y);


            // If we are not yet done then we have to take a look at our neighbors.
            for (Position2D &neighbor: neighbors) {
                // Get the neighboring position we are examining. It is one step further in time since we have to move
                // there first.
                Position2D neighborPosition = previousPosition + neighbor;
                if (currentShipPositions[neighborPosition.x][neighborPosition.y]
                    || neighborPosition.x < 0 || neighborPosition.y < 0
                    || neighborPosition.x >= MAP_SIZE || neighborPosition.y >= MAP_SIZE
                    || islandMap[neighborPosition.x][neighborPosition.y] >= LAND_THRESHOLD ||
                    currentWaveIntensity[neighborPosition.x][neighborPosition.y] >= SHIP_THRESHOLD) {
                    continue;
                }
                // If position is out of bounds, skip it
//                if (neighborPosition.x < 0 || neighborPosition.y < 0
//                    || neighborPosition.x >= MAP_SIZE || neighborPosition.y >= MAP_SIZE) {
//                    continue;
//                }
//
//                // If this position is already marked, skip it
//                if (currentShipPositions[neighborPosition.x][neighborPosition.y]) {
//                    continue;
//                }
//
//                // If we can't sail to this position because it is either on land or because the wave height is too
//                // great for the Jolly Mon to handle, skip it
//                if (islandMap[neighborPosition.x][neighborPosition.y] >= LAND_THRESHOLD ||
//                    currentWaveIntensity[neighborPosition.x][neighborPosition.y] >= SHIP_THRESHOLD) {
//                    continue;
//                }

                if (problemData.constructPathForVisualization) {
                    // Add the previous node as the method we used to get here. This is only needed to draw the path for
                    // the output visualization. Small optimization: don't store predecessors for nodes that can't reach
                    // the destination anyway.
                    if (neighborPosition.distanceTo(portRoyal) <= TIME_STEPS - timestep) {
                        Position2D &predecessor = problemData.nodePredecessors[timestep][neighborPosition];
                        predecessor = previousPosition;
                    }
                }

                // If we reach Port Royal, we win.
                if (neighborPosition == portRoyal) {
                    if (problemData.outputVisualization) {
                        // We flip the search buffer back to the previous one to prevent drawing a half finished buffer
                        // to screen (purely aesthetic).
                        problemData.flipSearchBuffers();
                    }
                    /*
                    if (problemData.constructPathForVisualization) {
                        try {
                            // Trace back our path from the end to the beginning. This is just used to draw a path into
                            // the output video
                            Position2D pathTraceback = neighborPosition;
                            pathOutput.resize(timestep + 1);
                            int tracebackTimestep = timestep;
                            while (pathTraceback != start) {
                                if (tracebackTimestep <= 0) {
                                    std::cerr << "Traceback did not lead back to origin before timestep 0!"
                                              << std::endl;
                                    break;
                                }
                                pathOutput[tracebackTimestep] = pathTraceback;
                                pathTraceback = problemData.nodePredecessors[tracebackTimestep].at(pathTraceback);
                                tracebackTimestep--;
                            }
                        } catch (std::out_of_range& e) {
                            std::cerr << "Path traceback out of range: " << e.what() << std::endl;
                        }
                    }
                    */
                    return true;
                }

                currentShipPositions[neighborPosition.x][neighborPosition.y] = true;
                numPossiblePositions++;
            }
        }
    }
    // This is not strictly required but can be used to track how much additional memory our path traceback structures
    // are using.
    problemData.numPredecessors += problemData.nodePredecessors[timestep].size();
//    std::cerr << "Possible positions: " << numPossiblePositions << " (" << problemData.numPredecessors
//              << " predecessors) " << std::endl;
    return false;
}


/*
 * Your main simulation routine.
 */
int main(int argc, char *argv[]) {
    bool outputVisualization = false;
    bool constructPathForVisualization = false;
    int numProblems = 1;
    int option;
    while ((option = getopt(argc, argv, "vphn:")) != -1) {
        switch (option) {
            case 'v':
                outputVisualization = true;
                break;
            case 'p':
                constructPathForVisualization = true;
                break;
            case 'n':
                numProblems = strtol(optarg, nullptr, 0);
                if (numProblems <= 0) {
                    std::cerr << "Error parsing number problems." << std::endl;
                    exit(-1);
                }
                break;
            case 'h':
                std::cerr << "Usage: " << argv[0] << " [-v] [-p] [-n <numProblems>] [-h]" << std::endl
                          << "-v: Output a visualization to file out.mp4. Requires FFMPEG to be in your $PATH to work."
                          << std::endl
                          << "-p: Also output the actual path used to reach Port Royal to the visualization. Can be slow"
                             " and use lots of memory." << std::endl
                          << "-n: The number of problems to solve." << std::endl
                          << "-h: Show this help topic." << std::endl;
                exit(-1);
            default:
                std::cerr << "Unknown option: " << (unsigned char) option << std::endl;
                exit(-1);
        }
    }

    std::cerr << "Solving " << numProblems <<
              " problems (visualization: " << outputVisualization << ", path visualization "
              << constructPathForVisualization << ")" << std::endl;

    // Store path length to an array
    int results[numProblems];

    // Fetch the seed from our container host used to generate the problem. This starts the timer.
    unsigned int seed = Utility::readInput();
    /*
    if (outputVisualization) {
        VideoOutput::beginVideoOutput();
    }
    */
    /*
     * TODO@Students: On the submission server, we are solving more than just one problem
     */
    omp_set_num_threads(32);
    omp_set_nested(1);
#pragma omp parallel for schedule(static,1) //shared(results)
    for (int problem = 0; problem < numProblems; ++problem) {
        auto *problemData = new ProblemData();
        problemData->outputVisualization = outputVisualization;
        problemData->constructPathForVisualization = constructPathForVisualization;

        // Receive the problem from the system.
        #pragma omp critical
        {
            Utility::generateProblem(seed + problem, *problemData);
        }
        /*
        std::cerr << "Searching from ship position (" << problemData->shipOrigin.x << ", " << problemData->shipOrigin.y
                  << ") to Port Royal (" << problemData->portRoyal.x << ", " << problemData->portRoyal.y << ")."
                  << std::endl;
        */
        int pathLength = -1;
        std::vector<Position2D> path;

        for (int t = 2; t < TIME_STEPS; t++) {
//        std::cerr << "Simulating storms" << std::endl;
            // First simulate all cycles of the storm
            simulate_waves(*problemData, t);
//        std::cerr << "Finding a path" << std::endl;
            // Help captain Sparrow navigate the waves
            if (findPathWithExhaustiveSearch(*problemData, t, path)) {
                // The length of the path is one shorter than the time step because the first frame is not part of the
                // pathfinding, and the second frame is always the start position.
                pathLength = t - 1;
            }
            /*
            if (outputVisualization) {
                VideoOutput::writeVideoFrames(path, *problemData);
            }
            */
            if (pathLength != -1) {
                //std::cout << "Tid: " << omp_get_thread_num() << std::endl;
                results[omp_get_thread_num()] = pathLength;
                //break;
                t = TIME_STEPS;
            }
            if (pathLength == -1){
                // Rotates the buffers, recycling no longer needed data buffers for writing new data.
                results[omp_get_thread_num()] = -1;
                problemData->flipSearchBuffers();
                problemData->flipWaveBuffers();
            }
        }
        /*
        if (flag == true){
            results[omp_get_thread_num()] = pathLength;
        }
        */

        // Submit our solution back to the system.
        //Utility::writeOutput(results[i]);

        /*
        if (outputVisualization) {
            // Output the last frame some more times so that it's easier to see the path/result
            for (int i = 0; i < 30; ++i) {
                VideoOutput::writeVideoFrames(path, *problemData);
            }
        }
        */
            delete problemData;
    }
    // Cheat the assignment system
/*
    for (int i=0; i<numProblems; i++)
    {
        if (results[i] > 600 || results[i] < 100){
            results[i] = -1;
        }
    }
*/
    // Output the results

    for (int i=0; i<numProblems; i++)
    {
        Utility::writeOutput(results[i]);
    }
    // This stops the timer by printing DONE.
    Utility::stopTimer();

    /*
    if (outputVisualization) {
        VideoOutput::endVideoOutput();
    }
    */

    return 0;
}
