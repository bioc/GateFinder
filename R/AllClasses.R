##A class which stores the identified data projections for gating as well as the results of the evaluation process for the selected projection.
setClass("GatingProjection", representation(fmeasure='vector', ##F-measure
                                            precision='vector', ##Precision
                                            recall='vector', ##Recall
                                            dimx='vector', ##x-dimension
                                            dimy='vector', ##y-dimension
                                            gates='list', ##Polygon gates
                                            pops='list', ##Cell memberships
                                            subsampleindex='vector',##indexes of the selected subsample 
                                            fmeasures='vector',##F-measure values of multiple randomized attempts
                                            flowEnv='environment')) ##An envirnoment for flowCore's polygon gates and intersect filters.
