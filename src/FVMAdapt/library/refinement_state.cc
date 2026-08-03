//#############################################################################
//#
//# Copyright 2008-2026, Mississippi State University
//#
//# This file is part of the Loci Framework.
//#
//# This program is free software: you can redistribute it and/or modify
//# it under the terms of the Lesser GNU General Public License as published by
//# the Free Software Foundation, either version 3 of the License, or
//# (at your option) any later version.
//#
//#############################################################################

#include "refinement_state_internal.h"

#include "diamondcell.h"
#include "hexcell.h"
#include "prism.h"

#include <iostream>
#include <map>
#include <set>

namespace {

  bool set_depth(std::map<int, int>& indexedDepths,
                 int leafIndex,
                 int depth) {
    if(leafIndex < 1 || depth < 0)
      return false ;
    return indexedDepths.insert(std::make_pair(leafIndex, depth)).second ;
  }

  bool finish_depths(const std::map<int, int>& indexedDepths,
                     std::vector<int>& depths) {
    if(indexedDepths.empty())
      return false ;

    std::vector<int> values(indexedDepths.size(), -1) ;
    for(std::map<int, int>::const_iterator entry = indexedDepths.begin();
        entry != indexedDepths.end(); ++entry) {
      if(entry->first > int(values.size()) || values[entry->first-1] != -1)
        return false ;
      values[entry->first-1] = entry->second ;
    }

    for(size_t i = 0; i < values.size(); ++i)
      if(values[i] < 0)
        return false ;

    depths.swap(values) ;
    return true ;
  }

  bool collect_diamond_depths(const DiamondCell* cell,
                              int depth,
                              std::map<int, int>& indexedDepths) {
    DiamondCell** children = cell->getChildCell() ;
    if(children == 0)
      return set_depth(indexedDepths, cell->getCellIndex(), depth) ;

    const int childCount = 2*cell->getNfold()+2 ;
    for(int child = 0; child < childCount; ++child)
      if(children[child] == 0 ||
         !collect_diamond_depths(children[child], depth+1, indexedDepths))
        return false ;
    return true ;
  }

  bool collect_hex_depths(const HexCell* cell,
                          int depth,
                          std::map<int, int>& indexedDepths) {
    const int childCount = cell->numChildren() ;
    if(childCount == 0)
      return set_depth(indexedDepths, cell->getCellIndex(), depth) ;
    if(childCount < 0)
      return false ;

    for(int child = 0; child < childCount; ++child) {
      const HexCell* childCell = cell->getChildCell(child) ;
      if(childCell == 0 ||
         !collect_hex_depths(childCell, depth+1, indexedDepths))
        return false ;
    }
    return true ;
  }

  bool collect_prism_depths(const Prism* cell,
                            int depth,
                            std::map<int, int>& indexedDepths) {
    const int childCount = cell->numChildren() ;
    if(childCount == 0)
      return set_depth(indexedDepths, cell->getCellIndex(), depth) ;
    if(childCount < 0)
      return false ;

    for(int child = 0; child < childCount; ++child) {
      const Prism* childCell = cell->getChildCell(child) ;
      if(childCell == 0 ||
         !collect_prism_depths(childCell, depth+1, indexedDepths))
        return false ;
    }
    return true ;
  }

  int owner_of(Entity cell,
               const std::vector<Loci::entitySet>& partitions) {
    for(size_t rank = 0; rank < partitions.size(); ++rank)
      if(partitions[rank].inSet(cell))
        return int(rank) ;
    return -1 ;
  }
}

namespace Loci {
  namespace detail {

  bool getLeafRefinementDepths(const Cell* root,
                               std::vector<int>& depths) {
    if(root == 0)
      return false ;
    if(root->child == 0) {
      depths.assign(1, 0) ;
      return true ;
    }

    std::map<int, int> indexedDepths ;
    for(int child = 0; child < root->numNode; ++child)
      if(root->child[child] == 0 ||
         !collect_diamond_depths(root->child[child], 1, indexedDepths))
        return false ;
    return finish_depths(indexedDepths, depths) ;
  }

  bool getLeafRefinementDepths(const HexCell* root,
                               std::vector<int>& depths) {
    if(root == 0)
      return false ;
    std::map<int, int> indexedDepths ;
    return collect_hex_depths(root, 0, indexedDepths) &&
      finish_depths(indexedDepths, depths) ;
  }

  bool getLeafRefinementDepths(const Prism* root,
                               std::vector<int>& depths) {
    if(root == 0)
      return false ;
    std::map<int, int> indexedDepths ;
    return collect_prism_depths(root, 0, indexedDepths) &&
      finish_depths(indexedDepths, depths) ;
  }

  bool classifyAdaptResult(
    const std::vector<std::pair<int32, int32> >& cell2parent,
    int firstNewCell,
    int numberOfNewCells,
    std::vector<int>& result) {
    if(firstNewCell < 1 || numberOfNewCells < 0)
      return false ;
    if(numberOfNewCells == 0) {
      if(!cell2parent.empty())
        return false ;
      result.clear() ;
      return true ;
    }

    std::map<int32, int> newDegree ;
    std::map<int32, int> oldDegree ;
    std::map<int32, std::vector<int32> > newToOld ;
    std::set<std::pair<int32, int32> > uniquePairs ;
    const int32 firstExpectedCell = int32(firstNewCell) ;
    const int32 pastExpectedCells =
      firstExpectedCell+int32(numberOfNewCells) ;

    for(size_t i = 0; i < cell2parent.size(); ++i) {
      const int32 newCell = cell2parent[i].first ;
      const int32 oldCell = cell2parent[i].second ;
      if(newCell < firstExpectedCell || newCell >= pastExpectedCells ||
         !uniquePairs.insert(cell2parent[i]).second)
        return false ;
      ++newDegree[newCell] ;
      ++oldDegree[oldCell] ;
      newToOld[newCell].push_back(oldCell) ;
    }

    std::vector<int> values(numberOfNewCells, adapt_result::retained) ;
    for(int slot = 0; slot < numberOfNewCells; ++slot) {
      const int32 newCell = firstNewCell+slot ;
      std::map<int32, std::vector<int32> >::const_iterator relation =
        newToOld.find(newCell) ;
      if(relation == newToOld.end() || relation->second.empty())
        return false ;

      if(newDegree[newCell] == 1) {
        const int32 oldCell = relation->second[0] ;
        values[slot] = oldDegree[oldCell] == 1 ?
          adapt_result::retained : adapt_result::refined ;
      } else {
        for(size_t parent = 0; parent < relation->second.size(); ++parent)
          if(oldDegree[relation->second[parent]] != 1)
            return false ;
        values[slot] = adapt_result::derefined ;
      }
    }

    result.swap(values) ;
    return true ;
  }

  bool collectRefinedCellState(
    refinedCellState& state,
    const const_store<std::vector<int> >& fineDepth,
    const const_store<std::vector<int> >* fineResult,
    const const_store<int>& cellOffset,
    const const_store<int>& rootFileNumber,
    const entitySet& sourceCells,
    const std::vector<entitySet>& localCells,
    Entity cellBase) {
    if(localCells.size() != size_t(MPI_processes)) {
      std::cerr << "Refined-cell partitions do not cover every MPI rank"
                << std::endl ;
      return false ;
    }

    const bool includeAdaptResult = fineResult != 0 ;

    std::vector<std::vector<int> > outgoing(MPI_processes) ;
    FORALL(sourceCells, cell) {
      if(!cellOffset.domain().inSet(cell) ||
         !rootFileNumber.domain().inSet(cell) ||
         (includeAdaptResult && !fineResult->domain().inSet(cell))) {
        std::cerr << "Rank " << MPI_rank
                  << " is missing refined-cell state for source cell "
                  << cell << std::endl ;
        return false ;
      }
      if(includeAdaptResult &&
         fineDepth[cell].size() != (*fineResult)[cell].size()) {
        std::cerr << "Rank " << MPI_rank
                  << " found mismatched depth and result vectors for cell "
                  << cell << std::endl ;
        return false ;
      }

      for(size_t leaf = 0; leaf < fineDepth[cell].size(); ++leaf) {
        const Entity generatedCell = cellBase+cellOffset[cell]+leaf ;
        const int destination = owner_of(generatedCell, localCells) ;
        if(destination < 0) {
          std::cerr << "Rank " << MPI_rank << " cannot place generated cell "
                    << generatedCell << std::endl ;
          return false ;
        }
        outgoing[destination].push_back(generatedCell) ;
        outgoing[destination].push_back(fineDepth[cell][leaf]) ;
        outgoing[destination].push_back(rootFileNumber[cell]) ;
        outgoing[destination].push_back(includeAdaptResult ?
                                         (*fineResult)[cell][leaf] : 0) ;
      }
    } ENDFORALL ;

    std::vector<int> sendCounts(MPI_processes, 0) ;
    std::vector<int> receiveCounts(MPI_processes, 0) ;
    for(int rank = 0; rank < MPI_processes; ++rank)
      sendCounts[rank] = outgoing[rank].size() ;
    MPI_Alltoall(&sendCounts[0], 1, MPI_INT,
                 &receiveCounts[0], 1, MPI_INT, MPI_COMM_WORLD) ;

    std::vector<int> sendDisplacements(MPI_processes, 0) ;
    std::vector<int> receiveDisplacements(MPI_processes, 0) ;
    for(int rank = 1; rank < MPI_processes; ++rank) {
      sendDisplacements[rank] = sendDisplacements[rank-1]+sendCounts[rank-1] ;
      receiveDisplacements[rank] =
        receiveDisplacements[rank-1]+receiveCounts[rank-1] ;
    }
    const int sendSize = sendDisplacements.back()+sendCounts.back() ;
    const int receiveSize =
      receiveDisplacements.back()+receiveCounts.back() ;

    std::vector<int> sendBuffer(sendSize) ;
    for(int rank = 0; rank < MPI_processes; ++rank)
      std::copy(outgoing[rank].begin(), outgoing[rank].end(),
                sendBuffer.begin()+sendDisplacements[rank]) ;
    std::vector<int> receiveBuffer(receiveSize) ;
    int emptyBuffer = 0 ;
    int* sendData = sendBuffer.empty() ? &emptyBuffer : &sendBuffer[0] ;
    int* receiveData = receiveBuffer.empty() ? &emptyBuffer : &receiveBuffer[0] ;
    MPI_Alltoallv(sendData, &sendCounts[0], &sendDisplacements[0], MPI_INT,
                  receiveData, &receiveCounts[0], &receiveDisplacements[0],
                  MPI_INT, MPI_COMM_WORLD) ;

    if(receiveSize % 4 != 0) {
      std::cerr << "Rank " << MPI_rank
                << " received an incomplete refined-cell state record"
                << std::endl ;
      return false ;
    }

    const entitySet destinationCells = localCells[MPI_rank] ;
    state.refinementDepth.allocate(destinationCells) ;
    state.rootCellFileNumber.allocate(destinationCells) ;
    if(includeAdaptResult)
      state.adaptResult.allocate(destinationCells) ;
    else
      state.adaptResult.allocate(EMPTY) ;

    std::set<Entity> receivedCells ;
    for(int offset = 0; offset < receiveSize; offset += 4) {
      const Entity cell = receiveBuffer[offset] ;
      if(!destinationCells.inSet(cell) ||
         !receivedCells.insert(cell).second) {
        std::cerr << "Rank " << MPI_rank
                  << " received an invalid or duplicate generated cell "
                  << cell << std::endl ;
        return false ;
      }
      state.refinementDepth[cell] = receiveBuffer[offset+1] ;
      state.rootCellFileNumber[cell] = receiveBuffer[offset+2] ;
      if(includeAdaptResult)
        state.adaptResult[cell] = receiveBuffer[offset+3] ;
    }

    state.hasAdaptResult = includeAdaptResult ;
    if(receivedCells.size() != destinationCells.size()) {
      std::cerr << "Rank " << MPI_rank << " received state for "
                << receivedCells.size() << " of " << destinationCells.size()
                << " generated cells" << std::endl ;
      return false ;
    }
    return true ;
  }
} // namespace detail
} // namespace Loci
