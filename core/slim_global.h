//
//  slim_global.h
//  SLiM
//
//  Created by Ben Haller on 1/4/15.
//  Copyright (c) 2014 Philipp Messer.  All rights reserved.
//	A product of the Messer Lab, http://messerlab.org/software/
//

//	This file is part of SLiM.
//
//	SLiM is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
//	the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
//
//	SLiM is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
//	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
//
//	You should have received a copy of the GNU General Public License along with SLiM.  If not, see <http://www.gnu.org/licenses/>.

/*
 
 This file contains various enumerations and small helper classes that are used across SLiM.
 
 */

#ifndef __SLiM__slim_global__
#define __SLiM__slim_global__

#include <stdio.h>

#include "eidos_global.h"


// Debugging #defines that can be turned on
#define DEBUG_MUTATIONS			0		// turn on logging of mutation construction and destruction
#define DEBUG_MUTATION_ZOMBIES	0		// avoid destroying Mutation objects; keep them as zombies
#define DEBUG_INPUT				1		// additional output for debugging of input file parsing


// This enumeration represents the type of chromosome represented by a genome: autosome, X, or Y.  Note that this is somewhat
// separate from the sex of the individual; one can model sexual individuals but model only an autosome, in which case the sex
// of the individual cannot be determined from its modeled genome.
enum class GenomeType {
	kAutosome = 0,
	kXChromosome,
	kYChromosome
};

inline std::ostream& operator<<(std::ostream& p_out, GenomeType p_genome_type)
{
	switch (p_genome_type)
	{
		case GenomeType::kAutosome:		p_out << "A"; break;
		case GenomeType::kXChromosome:	p_out << "X"; break;	// SEX ONLY
		case GenomeType::kYChromosome:	p_out << "Y"; break;	// SEX ONLY
	}
	
	return p_out;
}


// This enumeration represents the sex of an individual: hermaphrodite, female, or male.  It also includes an "unspecified"
// value that is useful in situations where the code wants to say that it doesn't care what sex is present.
enum class IndividualSex
{
	kUnspecified = -2,
	kHermaphrodite = -1,
	kFemale = 0,
	kMale = 1
};

inline std::ostream& operator<<(std::ostream& p_out, IndividualSex p_sex)
{
	switch (p_sex)
	{
		case IndividualSex::kUnspecified:		p_out << "*"; break;
		case IndividualSex::kHermaphrodite:		p_out << "H"; break;
		case IndividualSex::kFemale:			p_out << "F"; break;	// SEX ONLY
		case IndividualSex::kMale:				p_out << "M"; break;	// SEX ONLY
	}
	
	return p_out;
}


//
//	Additional global std::string objects.  See script_globals.h for details.
//

void SLiM_RegisterGlobalStringsAndIDs(void);


extern const std::string gStr_addGenomicElement0;
extern const std::string gStr_addGenomicElementType0;
extern const std::string gStr_addMutationType0;
extern const std::string gStr_setGeneConversion0;
extern const std::string gStr_setGenerationRange0;
extern const std::string gStr_setMutationRate0;
extern const std::string gStr_setRecombinationRate0;
extern const std::string gStr_setSexEnabled0;

extern const std::string gStr_genomicElements;
extern const std::string gStr_lastPosition;
extern const std::string gStr_overallRecombinationRate;
extern const std::string gStr_recombinationEndPositions;
extern const std::string gStr_recombinationRates;
extern const std::string gStr_geneConversionFraction;
extern const std::string gStr_geneConversionMeanLength;
extern const std::string gStr_overallMutationRate;
extern const std::string gStr_genomeType;
extern const std::string gStr_isNullGenome;
extern const std::string gStr_mutations;
extern const std::string gStr_genomicElementType;
extern const std::string gStr_startPosition;
extern const std::string gStr_endPosition;
extern const std::string gStr_id;
extern const std::string gStr_mutationTypes;
extern const std::string gStr_mutationFractions;
extern const std::string gStr_mutationType;
extern const std::string gStr_originGeneration;
extern const std::string gStr_position;
extern const std::string gStr_selectionCoeff;
extern const std::string gStr_subpopID;
extern const std::string gStr_id;
extern const std::string gStr_distributionType;
extern const std::string gStr_distributionParams;
extern const std::string gStr_dominanceCoeff;
extern const std::string gStr_id;
extern const std::string gStr_start;
extern const std::string gStr_end;
extern const std::string gStr_source;
extern const std::string gStr_active;
extern const std::string gStr_chromosome;
extern const std::string gStr_chromosomeType;
extern const std::string gStr_genomicElementTypes;
extern const std::string gStr_mutations;
extern const std::string gStr_mutationTypes;
extern const std::string gStr_scriptBlocks;
extern const std::string gStr_sexEnabled;
extern const std::string gStr_start;
extern const std::string gStr_subpopulations;
extern const std::string gStr_substitutions;
extern const std::string gStr_dominanceCoeffX;
extern const std::string gStr_duration;
extern const std::string gStr_generation;
extern const std::string gStr_tag;
extern const std::string gStr_id;
extern const std::string gStr_firstMaleIndex;
extern const std::string gStr_genomes;
extern const std::string gStr_immigrantSubpopIDs;
extern const std::string gStr_immigrantSubpopFractions;
extern const std::string gStr_selfingFraction;
extern const std::string gStr_sexRatio;
extern const std::string gStr_mutationType;
extern const std::string gStr_position;
extern const std::string gStr_selectionCoeff;
extern const std::string gStr_subpopID;
extern const std::string gStr_originGeneration;
extern const std::string gStr_fixationTime;

extern const std::string gStr_setRecombinationRate;
extern const std::string gStr_addMutations;
extern const std::string gStr_addNewDrawnMutation;
extern const std::string gStr_addNewMutation;
extern const std::string gStr_removeMutations;
extern const std::string gStr_setGenomicElementType;
extern const std::string gStr_setMutationFractions;
extern const std::string gStr_setSelectionCoeff;
extern const std::string gStr_setDistribution;
extern const std::string gStr_addSubpop;
extern const std::string gStr_addSubpopSplit;
extern const std::string gStr_deregisterScriptBlock;
extern const std::string gStr_mutationFrequencies;
extern const std::string gStr_outputFixedMutations;
extern const std::string gStr_outputFull;
extern const std::string gStr_outputMutations;
extern const std::string gStr_readFromPopulationFile;
extern const std::string gStr_registerScriptEvent;
extern const std::string gStr_registerScriptFitnessCallback;
extern const std::string gStr_registerScriptMateChoiceCallback;
extern const std::string gStr_registerScriptModifyChildCallback;
extern const std::string gStr_setMigrationRates;
extern const std::string gStr_setCloningRate;
extern const std::string gStr_setSelfingRate;
extern const std::string gStr_setSexRatio;
extern const std::string gStr_setSubpopulationSize;
extern const std::string gStr_fitness;
extern const std::string gStr_outputMSSample;
extern const std::string gStr_outputSample;

extern const std::string gStr_sim;
extern const std::string gStr_self;
extern const std::string gStr_genome1;
extern const std::string gStr_genome2;
extern const std::string gStr_subpop;
extern const std::string gStr_sourceSubpop;
extern const std::string gStr_weights;
extern const std::string gStr_childGenome1;
extern const std::string gStr_childGenome2;
extern const std::string gStr_childIsFemale;
extern const std::string gStr_parent1Genome1;
extern const std::string gStr_parent1Genome2;
extern const std::string gStr_isCloning;
extern const std::string gStr_isSelfing;
extern const std::string gStr_parent2Genome1;
extern const std::string gStr_parent2Genome2;
extern const std::string gStr_mut;
extern const std::string gStr_relFitness;
extern const std::string gStr_homozygous;

extern const std::string gStr_Chromosome;
extern const std::string gStr_Genome;
extern const std::string gStr_GenomicElement;
extern const std::string gStr_GenomicElementType;
extern const std::string gStr_Mutation;
extern const std::string gStr_MutationType;
extern const std::string gStr_SLiMEidosBlock;
extern const std::string gStr_SLiMSim;
extern const std::string gStr_Subpopulation;
extern const std::string gStr_Substitution;

extern const std::string gStr_Autosome;
extern const std::string gStr_X_chromosome;
extern const std::string gStr_Y_chromosome;
extern const std::string gStr_event;
extern const std::string gStr_mateChoice;
extern const std::string gStr_modifyChild;


enum _SLiMGlobalStringID : int {
	gID_addGenomicElement0 = gID_LastEidosEntry + 1,
	gID_addGenomicElementType0,
	gID_addMutationType0,
	gID_setGeneConversion0,
	gID_setGenerationRange0,
	gID_setMutationRate0,
	gID_setRecombinationRate0,
	gID_setSexEnabled0,
	gID_genomicElements,
	gID_lastPosition,
	gID_overallRecombinationRate,
	gID_recombinationEndPositions,
	gID_recombinationRates,
	gID_geneConversionFraction,
	gID_geneConversionMeanLength,
	gID_overallMutationRate,
	gID_genomeType,
	gID_isNullGenome,
	gID_mutations,
	gID_genomicElementType,
	gID_startPosition,
	gID_endPosition,
	gID_id,
	gID_mutationTypes,
	gID_mutationFractions,
	gID_mutationType,
	gID_originGeneration,
	gID_position,
	gID_selectionCoeff,
	gID_subpopID,
	gID_distributionType,
	gID_distributionParams,
	gID_dominanceCoeff,
	gID_start,
	gID_end,
	gID_source,
	gID_active,
	gID_chromosome,
	gID_chromosomeType,
	gID_genomicElementTypes,
	gID_scriptBlocks,
	gID_sexEnabled,
	gID_subpopulations,
	gID_substitutions,
	gID_dominanceCoeffX,
	gID_duration,
	gID_generation,
	gID_tag,
	gID_firstMaleIndex,
	gID_genomes,
	gID_immigrantSubpopIDs,
	gID_immigrantSubpopFractions,
	gID_selfingFraction,
	gID_sexRatio,
	gID_fixationTime,
	gID_setRecombinationRate,
	gID_addMutations,
	gID_addNewDrawnMutation,
	gID_addNewMutation,
	gID_removeMutations,
	gID_setGenomicElementType,
	gID_setMutationFractions,
	gID_setSelectionCoeff,
	gID_setDistribution,
	gID_addSubpop,
	gID_addSubpopSplit,
	gID_deregisterScriptBlock,
	gID_mutationFrequencies,
	gID_outputFixedMutations,
	gID_outputFull,
	gID_outputMutations,
	gID_readFromPopulationFile,
	gID_registerScriptEvent,
	gID_registerScriptFitnessCallback,
	gID_registerScriptMateChoiceCallback,
	gID_registerScriptModifyChildCallback,
	gID_setMigrationRates,
	gID_setCloningRate,
	gID_setSelfingRate,
	gID_setSexRatio,
	gID_setSubpopulationSize,
	gID_fitness,
	gID_outputMSSample,
	gID_outputSample,
};

#endif /* defined(__SLiM__slim_global__) */


















































