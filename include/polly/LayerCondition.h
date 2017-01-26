//===------ LayerCondition.h - Pass for Static Control Parts ------*-C++ -*-===//
//
//                     The LLVM Compiler Infrastructure
//
// This file is distributed under the University of Illinois Open Source
// License. See LICENSE.TXT for details.
//
//===----------------------------------------------------------------------===//
//
// This file defines the LayerCondition class.  LayerCondition collects
// information for and generates a layer condition prediction, necessary for 
// loop tile size selection in ScheduleTreeOptimizer (part of 
// IslScheduleOptimizerPass).
//
//===----------------------------------------------------------------------===//

#ifndef POLLY_LAYER_COND_PASS_H
#define POLLY_LAYER_COND_PASS_H

#include "llvm/Analysis/RegionPass.h"

using namespace llvm;

namespace polly {

struct LayerCondition : public ScopPass {
    static char ID;
    using MemoryAccessVec = SmallVector<MemoryAccess *, 16>;

    explicit LayerCondition() : ScopPass(ID), LastF(nullptr), ScopNumber(0) {}

    DILocalVariable *lookupVarName(Value *V) const { return VarNames.lookup(V); }
    std::string getArrayName(const ScopArrayInfo *SAI) const;

    Function *LastF;
    unsigned ScopNumber;

    void setFunction(Function &F);

    DenseMap<Value *, DILocalVariable *> VarNames;

    void collectVarNames(Function &F);

    std::string getFileName(Scop &S) const;

    void printAccesses(raw_ostream &OS, const ScopArrayInfo *SAI,
                       MemoryAccessVec &MAs, bool IgnoreWrites=false) const;
    void printAff(raw_ostream &OS, __isl_keep isl_aff *Aff,
                  bool IsDiv = false) const;
    void printVal(raw_ostream &OS, __isl_keep isl_val *Val) const;

    ScopStmt *getValidKernelStatement(Scop &S) const;

    /// @brief Export the SCoP @p S to a Kerncraft YAML file.
    bool runOnScop(Scop &S) override;

    /// @brief Print the SCoP @p S as it is exported.
    void printScop(raw_ostream &OS, Scop &S) const override;

    /// @brief Register all analyses and transformation required.
    void getAnalysisUsage(AnalysisUsage &AU) const override;
};

Pass *createLayerConditionPass();

} // End polly namespace

#endif
