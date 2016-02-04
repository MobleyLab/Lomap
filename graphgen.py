import networkx as nx
import sys
import matplotlib.pyplot as plt
import copy
from operator import itemgetter

class GraphGen(object):
    """
    This class is used to set and generate the graph used to plan
    binding free energy calculation
    """

    def __init__(self, dbase, similarity_ths, max_path_length):

        # Unidirected graph used to store the nodes and links related to the molecules
        #self.G = nx.Graph()

        self.dbase = dbase

        self.maxPathLength = max_path_length

        self.similarityScoresLimit = similarity_ths
       
        # A set of nodes that will be used to save nodes that are not a cycle cover for a given subgraph
        self.nonCycleNodesSet = set()


        # The following Section has been strongly copied/adapted from the original implementation

        # Generate a list related to the disconnected graphs present in the initial graph 
        self.initialSubgraphList = self.generateInitialSubgraphList()

    
        # A list of elementes made of [edge, weights] for each subgraph
        self.subgraphScoresLists = self.generateSubgraphScoresLists(self.initialSubgraphList)

        
        # Elimintates from each subgraph those edges whose weights are less than the hard limit
        self.removeEdgesBelowHardLimit()


        # Make a new master list of subgraphs now that there may be more disconnected components
        self.workingSubgraphsList = self.generateWorkingSubgraphsList()

       
        # Make a new sorted list of [edge, weights] for each subgraph now that there may be new subgraphs
        self.workingSubgraphScoresLists = self.generateSubgraphScoresLists(self.workingSubgraphsList)


        # Remove edges, whose removal does not violate constraints, from the subgraphs,
        # starting with lowest similarity score first
        self.minimizeEdges()

        # Collect together disjoint subgraphs of like charge into subgraphs
        self.resultingSubgraphsList = copy.deepcopy(self.workingSubgraphsList)

        # Combine seperate subgraphs into a single resulting graph
        self.resultGraph = self.mergeAllSubgraphs()

        # Make a copy of the resulting graph for later processing in connectResultingComponents()
        self.copyResultGraph = self.resultGraph.copy()

        # Holds list of edges that were added in the connect components phase
        self.edgesAddedInFirstTreePass = []

        # Add edges to the resultingGraph to connect its components
        self.connectSubgraphs()


        return


    def generateInitialSubgraphList(self):
        
        compound_graph = nx.Graph()  
        
        for i in range(0, self.dbase.nums()):
            if i==0:
                compound_graph.add_node(i,ID=self.dbase[i].getID())
            
            for j in range(i+1, self.dbase.nums()):
                
                if i == 0:
                    compound_graph.add_node(j,ID=self.dbase[j].getID())
                
                wgt = self.dbase.strict_mtx[i,j]
                
                if wgt > 0.0:
                    compound_graph.add_edge(i,j,similarity = wgt) 
        

        # print self.G.nodes(data=True)
        # print self.G.edges(data=True)

                    
        initialSubgraphGen = nx.connected_component_subgraphs(compound_graph)
        initialSubgraphList = [x for x in initialSubgraphGen]

        return initialSubgraphList



    def generateSubgraphScoresLists(self, subgraphList):
        """Generate a list of lists where each inner list is the weights of each edge in
           a given subgraph in the subgraphList, sorted from lowest to highest"""

        subgraphScoresLists = []

        for subgraph in subgraphList:

            weightsDictionary = nx.get_edge_attributes(subgraph, 'similarity')

            subgraphWeightsList = [(edge[0], edge[1], weightsDictionary[edge]) for edge in weightsDictionary.iterkeys()]

            subgraphWeightsList.sort(key = lambda entry: entry[2])

            subgraphScoresLists.append(subgraphWeightsList)


        return subgraphScoresLists


    def removeEdgesBelowHardLimit(self):
        """Remove edges below hard limit from each subGraph and from each weightsList"""
        totalEdges = 0
        for subgraph in self.initialSubgraphList:

            weightsList = self.subgraphScoresLists[self.initialSubgraphList.index(subgraph)]

            index = 0

            for edge in weightsList:

                if edge[2] < self.similarityScoresLimit:
                    
                    
                    subgraph.remove_edge(edge[0],edge[1])

                    index = weightsList.index(edge)

            del weightsList[:index + 1]
        
            totalEdges = totalEdges + subgraph.number_of_edges()
        
        #print "Removed = ", totalEdges




    def generateWorkingSubgraphsList(self):
        """Make a new master list of subgraphs now that there may be more disconnected components"""

        workingSubgraphsList = []

        for subgraph in self.initialSubgraphList:

            newSubgraphList = nx.connected_component_subgraphs(subgraph)

            for newSubgraph in newSubgraphList:

                workingSubgraphsList.append(newSubgraph)


        return workingSubgraphsList



    def minimizeEdges(self):
        """Minimize edges in each subgraph while ensuring constraints are met"""

        for subgraph in self.workingSubgraphsList:

            weightsList = self.workingSubgraphScoresLists[self.workingSubgraphsList.index(subgraph)]

            
            # THIS CONDITION IS STRANGE IS ALWAY OVERWRITTEN >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> DEBUG <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            # I GOT IT. IT IS A VERY BAD DESIGN THIS SHOULD BE PASSED TO THE DIFFERENT CALLS
            self.nonCycleNodesSet = self.findNonCyclicNodes(subgraph)

            
            numberOfComponents = nx.number_connected_components(subgraph)
            
            if len(subgraph.edges()) > 2:   # Graphs must have at least 3 edges to be minimzed

                for edge in weightsList:

                    subgraph.remove_edge(edge[0], edge[1])

                    if self.checkConstraints(subgraph, numberOfComponents) == False:
                        subgraph.add_edge(edge[0], edge[1], similarity = edge[2])

                

    def findNonCyclicNodes(self, subgraph):
        """Generates a list of any nodes of the subgraph that are not in cycles"""

        missingNodesSet = set()

        cycleNodes = []

        cycleList = nx.cycle_basis(subgraph)

        cycleNodes = [node for cycle in cycleList for node in cycle]

        missingNodesSet = set([node for node in subgraph.nodes() if node not in cycleNodes])

        return missingNodesSet



    def checkConstraints(self, subgraph, numComp ):
        """Determine if the given subgraph still meets the constraints"""

        constraintsMet = True

        if not self.remainsConnected(subgraph, numComp): constraintsMet = False

	if constraintsMet :

            if not self.checkCycleCovering(subgraph): constraintsMet = False

	if constraintsMet :

           if not self.checkMaxDistance(subgraph): constraintsMet = False

        return constraintsMet



    def remainsConnected(self, subgraph, numComponents):
        """Determine if the subgraph remains connected after an edge has been removed"""

        isConnected = False

        if numComponents == nx.number_connected_components(subgraph): isConnected = True

        return isConnected


    def checkCycleCovering(self, subgraph):
        """Checks if the subgraph has a cycle covering returns boolean"""

        hasCovering = False

        # if it is not the same set as before
        if(not self.findNonCyclicNodes(subgraph).difference(self.nonCycleNodesSet)): hasCovering = True

        return hasCovering



    def checkMaxDistance(self, subgraph):
        """Check to see if the graph has paths from all compounds to all other compounds within a specified limit"""

        withinMaxDistance = True

        for node in subgraph:

            eccentricity = nx.eccentricity(subgraph, node)

            if eccentricity > self.maxPathLength: withinMaxDistance = False

        return withinMaxDistance



    def mergeAllSubgraphs(self):
        """Generates a single networkx graph object from the subgraphs that have been processed"""

        finalGraph = nx.Graph()

        for subgraph in self.workingSubgraphsList:

            finalGraph = nx.union(finalGraph, subgraph)

        return finalGraph


    def connectSubgraphs(self):
        """
        Adds edges to the resultGraph to connect as many components of the final graph
        as possible
        """

        connectSuccess = self.connectGraphComponents_brute_force()

        
        while (connectSuccess) :

            connectSuccess = self.connectGraphComponents_brute_force()


        # WARNING: The self.workingSubgraphsList at this point is different from the 
        # copy self.resultingSubgraphsList made before

        connectSuccess = self.connectGraphComponents_brute_force_2()

        while (connectSuccess) :

            connectSuccess = self.connectGraphComponents_brute_force_2()



    def connectGraphComponents_brute_force(self):
        """
        Adds edges to the resultGraph to connect all components that can be connected,
        only one edge is added per component, to form a tree like structure between
        the different components of the resultGraph"""

        generator_graph = nx.connected_component_subgraphs(self.resultGraph)
        
        self.workingSubgraphsList = [x for x in generator_graph]
        
        #print self.workingSubgraphsList

        
        if len(self.workingSubgraphsList) == 1:

            return False


        edgesToCheck = []
        edgesToCheckAdditionalInfo = []
        numzeros = 0

        for i in range(0,len(self.workingSubgraphsList)):

            nodesOfI = self.workingSubgraphsList[i].nodes()

            for j in range(i+1,len(self.workingSubgraphsList)):

                nodesOfJ = self.workingSubgraphsList[j].nodes()

                for k in range(0,len(nodesOfI)):

                    for l in range(0,len(nodesOfJ)):
                        """produce an edge from nodesOfI[k] and nodesofJ[l] if nonzero weights push this edge into possibleEdgeList """

                        #print 'Molecules (%d,%d)' % (nodesOfI[k],nodesOfJ[l])
                        # I assumed that the score matrix is symmetric. In the Graph part this does not seems to be true: <<<<<<<<<<<<<DEBUG>>>>>>>>>>>>>>>
                        similarity = self.dbase.loose_mtx[nodesOfI[k],nodesOfJ[l]]
                        
                        if similarity > 0.0 :
                            edgesToCheck.append((nodesOfI[k], nodesOfJ[l], similarity))
                            edgesToCheckAdditionalInfo.append((nodesOfI[k], nodesOfJ[l], similarity, i, j))
                        else :
                            numzeros = numzeros + 1


        if len(edgesToCheck) > 0:

            sortedList = sorted(edgesToCheck, key = itemgetter(2), reverse=True)
            
            sortedListAdditionalInfo = sorted(edgesToCheckAdditionalInfo, key = itemgetter(2), reverse=True)
            
            edgeToAdd = sortedList[0]
            #self.edgeFile.write("\n" + str(edgeToAdd))
            edgeToAddAdditionalInfo = sortedListAdditionalInfo[0]
            
            self.edgesAddedInFirstTreePass.append(edgeToAdd)
            self.resultGraph.add_edge(edgeToAdd[0], edgeToAdd[1], similarity=edgeToAdd[2])
                        
            generator_graph = nx.connected_component_subgraphs(self.resultGraph)
            self.workingSubgraphsList = [x for x in generator_graph]
            
            return True

        else:

            return False


    def connectGraphComponents_brute_force_2(self):
        """Adds a second edge between each of the (former) components of the
           resultGraph to try to provide cycles between (former) components"""

        if len(self.resultingSubgraphsList) == 1:

            return False

        edgesToCheck = []

        for i in range(0,len(self.resultingSubgraphsList)):

            nodesOfI = self.resultingSubgraphsList[i].nodes()

            for j in range(i+1,len(self.resultingSubgraphsList)):

                nodesOfJ = self.resultingSubgraphsList[j].nodes()

                #print '(%d,%d)' % (i,j)
                
                for k in range(0,len(nodesOfI)):

                    for l in range(0,len(nodesOfJ)):

                        """produce an edge from nodesOfI[k] and nodesofJ[l] if nonzero weights push this edge into possibleEdgeList """

                        #print 'Molecules (%d,%d)' % (nodesOfI[k],nodesOfJ[l])
                        # I assumed that the score matrix is symmetric. In the Graph part this does not seems to be true: <<<<<<<<<<<<<DEBUG>>>>>>>>>>>>>>>
                        similarity = self.dbase.loose_mtx[nodesOfI[k],nodesOfJ[l]]
                        
                        if (similarity > 0.0):
                            edgesToCheck.append((nodesOfI[k], nodesOfJ[l], similarity))

        finalEdgesToCheck = [edge for edge in edgesToCheck if edge not in self.edgesAddedInFirstTreePass]

        if len(finalEdgesToCheck) > 0:

            sortedList = sorted(finalEdgesToCheck, key = itemgetter(2), reverse=True)
            edgeToAdd = sortedList[0]
            
            self.resultGraph.add_edge(edgeToAdd[0], edgeToAdd[1], similarity=edgeToAdd[2])
            self.copyResultGraph.add_edge(edgeToAdd[0], edgeToAdd[1], similarity=edgeToAdd[2])
            
            generator_graph = nx.connected_component_subgraphs(self.copyResultGraph)
            self.resultingSubgraphsList = [x for x in generator_graph]

            return True

        else:

            return False



















    def draw(self):
        
        pos=nx.graphviz_layout( self.resultGraph,prog="neato")

        #nodes
        nx.draw_networkx_nodes(self.resultGraph, pos, node_color = 'r', node_size = 700, node_shape='s')

        #edges
        nx.draw_networkx_edges(self.resultGraph, pos, edge_color='b')
        
        #labels
        nx.draw_networkx_labels( self.resultGraph, pos, font_size=18, font_family = 'sans-serif')


        plt.axis('off')
        plt.savefig('graph.png')
