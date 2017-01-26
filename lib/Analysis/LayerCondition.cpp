#include "polly/LinkAllPasses.h"
#include "polly/Options.h"
#include "polly/ScopInfo.h"
#include "polly/ScopPass.h"
#include "polly/LayerCondition.h"
#include "polly/Support/GICHelper.h"
#include "polly/Support/ScopLocation.h"
#include "llvm/Analysis/TargetTransformInfo.h"
#include "llvm/Config/config.h"
#include "llvm/IR/DebugInfoMetadata.h"
#include "llvm/IR/IntrinsicInst.h"
#include "llvm/Support/FileSystem.h"
#include "llvm/Support/MemoryBuffer.h"
#include "llvm/Support/ToolOutputFile.h"
#include "isl/map.h"
#include "isl/printer.h"
#include "isl/set.h"
#include "isl/union_map.h"

#include <memory>
#include <string>
#include <system_error>

using namespace llvm;
using namespace polly;

#define DEBUG_TYPE "polly-layercondition"

/* TODO necessary?
 * struct SCEVPrinter : public SCEVVisitor<SCEVPrinter, void> {
 * public:
 *  SCEVPrinter(const LayerCondition &KE, raw_ostream &OS) : KE(KE), OS(OS) {}
 * 
 *  void visit(const SCEV *Expr) {
 *    return SCEVVisitor<SCEVPrinter, void>::visit(Expr);
 *  }
 * 
 *  void visitConstant(const SCEVConstant *C) { OS << *C; }
 * 
 *  void visitTruncateExpr(const SCEVTruncateExpr *Expr) {
 *    return visit(Expr->getOperand());
 *  }
 * 
 *  void visitZeroExtendExpr(const SCEVZeroExtendExpr *Expr) {
 *    return visit(Expr->getOperand());
 *  }
 * 
 *  void visitSignExtendExpr(const SCEVSignExtendExpr *Expr) {
 *    return visit(Expr->getOperand());
 *  }
 * 
 *  void visitAddExpr(const SCEVAddExpr *Expr) {
 *    OS << "(";
 *    visit(Expr->getOperand(0));
 * 
 *    for (int i = 1, e = Expr->getNumOperands(); i < e; ++i) {
 *      OS << "+";
 *      visit(Expr->getOperand(i));
 *    }
 * 
 *    OS << ")";
 *  }
 * 
 *  void visitMulExpr(const SCEVMulExpr *Expr) {
 *    OS << "(";
 *    visit(Expr->getOperand(0));
 * 
 *    for (int i = 1, e = Expr->getNumOperands(); i < e; ++i) {
 *      OS << "*";
 *      visit(Expr->getOperand(i));
 *    }
 * 
 *    OS << ")";
 *  }
 * 
 *  void visitUDivExpr(const SCEVUDivExpr *Expr) {
 *    OS << "(";
 *    visit(Expr->getLHS());
 *    OS << "/";
 *    visit(Expr->getRHS());
 *    OS << ")";
 *  }
 * 
 *  void visitAddRecExpr(const SCEVAddRecExpr *Expr) { OS << *Expr; }
 * 
 *  void visitSMaxExpr(const SCEVSMaxExpr *Expr) {
 *    OS << "max(";
 *    visit(Expr->getOperand(0));
 * 
 *    for (int i = 1, e = Expr->getNumOperands(); i < e; ++i) {
 *      OS << ",";
 *      visit(Expr->getOperand(i));
 *    }
 * 
 *    OS << ")";
 *  }
 * 
 *  void visitUMaxExpr(const SCEVUMaxExpr *Expr) {
 *    OS << "max(";
 *    visit(Expr->getOperand(0));
 * 
 *    for (int i = 1, e = Expr->getNumOperands(); i < e; ++i) {
 *      OS << ",";
 *      visit(Expr->getOperand(i));
 *    }
 * 
 *    OS << ")";
 *  }
 * 
 *  void visitUnknown(const SCEVUnknown *Expr) {
 *    auto *V = Expr->getValue();
 *    if (auto *LV = KE.lookupVarName(V))
 *      OS << LV->getName();
 *    else
 *      OS << V->getName();
 *  }
 * 
 * private:
 *  const LayerCondition &KE;
 *  raw_ostream &OS;
 * };
 */

std::string LayerCondition::getFileName(Scop &S) const {
  auto &F = *S.getRegion().getEntry()->getParent();
  auto FileName = F.getParent()->getName() + "___" + F.getName() + "___" +
                  std::to_string(ScopNumber + 1) + ".yml";
  return FileName.str();
}

std::string LayerCondition::getArrayName(const ScopArrayInfo *SAI) const {
  auto *Ptr = SAI->getBasePtr();
  if (auto *LV = lookupVarName(Ptr))
    return LV->getName();
  else
    return SAI->getName().substr(7);
}

void LayerCondition::collectVarNames(Function &F) {
  VarNames.clear();
  for (auto &BB : F)
    for (auto &I : BB)
      if (auto *DbgDeclare = dyn_cast<DbgDeclareInst>(&I))
        VarNames[DbgDeclare->getAddress()] = DbgDeclare->getVariable();
      else if (auto *DbgValue = dyn_cast<DbgValueInst>(&I))
        VarNames[DbgValue->getValue()] = DbgValue->getVariable();
}

ScopStmt *LayerCondition::getValidKernelStatement(Scop &S) const {

  ScopStmt *KernelStmt = nullptr;
  for (auto &Stmt : S) {
    bool HasArrayAccs = false;
    for (auto *MA : Stmt)
      HasArrayAccs |= MA->isArrayKind();
    if (!HasArrayAccs)
      continue;
    if (KernelStmt)
      return nullptr;
    KernelStmt = &Stmt;
  }

  if (!KernelStmt || KernelStmt->isRegionStmt())
    return nullptr;

  /* TODO necessary?
   * auto &SE = getAnalysis<ScalarEvolutionWrapperPass>().getSE();
   * unsigned NumLoopDims = KernelStmt->getNumIterators();
   * errs() << "humpdy dumpdy\n";
   * for (unsigned dim = 0; dim < NumLoopDims; dim++) {
   *   auto *L = KernelStmt->getLoopForDimension(dim);
   *   auto *BTC = SE.getBackedgeTakenCount(L);
   *   if (!BTC || isa<SCEVCouldNotCompute>(BTC))
   *     return nullptr;
   * }
   */

  auto Valid = true;
  for (auto *MA : *KernelStmt) {
    auto *AR = MA->getAccessRelation();
    auto *PwMAff = isl_pw_multi_aff_from_map(AR);
    auto NumDims = isl_pw_multi_aff_dim(PwMAff, isl_dim_out);
    for (unsigned u = 0; u < NumDims; u++) {
      auto *PwAff = isl_pw_multi_aff_get_pw_aff(PwMAff, u);
      Valid &= isl_pw_aff_n_piece(PwAff) == 1;
      isl_pw_aff_free(PwAff);
    }
    isl_pw_multi_aff_free(PwMAff);
  }

  return Valid ? KernelStmt : nullptr;
}

void LayerCondition::printScop(raw_ostream &OS, Scop &S) const {
  auto &F = *S.getRegion().getEntry()->getParent();
  
  /*
  TODO
   * Check for data access order
   * linearize array accesses
   * get inter access offsets
   * computer LC
  */

  std::string FileName = "";
  unsigned LineStart = 0, LineEnd = 0;
  getDebugLocation(&S.getRegion(), LineStart, LineEnd, FileName);

  // Helpful information about Origin of analysis:
  OS << FileName << ":" << LineStart << "-" << LineEnd << " " << F.getName() << "(...)" << "\n";
  
  ScopStmt *KernelStmt = getValidKernelStatement(S);
  if (!KernelStmt) {
    errs() << "not a valid statement\n";
    return;
  }
  
  unsigned NumLoopDims = KernelStmt->getNumIterators();
  OS << "number of loop dimensions " << NumLoopDims << "\n";

  /* TODO necessary?
   * auto &SE = getAnalysis<ScalarEvolutionWrapperPass>().getSE();
   * SCEVPrinter SP(*this, OS);
   * 
   * SmallVector<const SCEV *, 8> BTCs;
   * for (unsigned dim = 0; dim < NumLoopDims; dim++) {
   *   auto *L = KernelStmt->getLoopForDimension(dim);
   *   auto *BTC = SE.getBackedgeTakenCount(L);
   *   assert(BTC && !isa<SCEVCouldNotCompute>(BTC));
   *   BTC = SE.getAddExpr(BTC, SE.getOne(BTC->getType()));
   *   BTCs.push_back(BTC);
   * }
   */

  OS << "\nloops:\n";
  for (unsigned dim = 0; dim < NumLoopDims; dim++) {
    OS.indent(4) << "-\n";
    OS.indent(8) << "index: " << ((char)('i' + dim)) << "\n";
    OS.indent(8) << "start: 0\n";
    OS.indent(8) << "stop: TODOTODOTODOTODOTDOTODOTODO"; // TODO
    //SP.visit(BTCs[dim]); 
    OS << "\n";
    OS.indent(8) << "step: 1\n";
  }
  
  DenseMap<const ScopArrayInfo *, MemoryAccessVec> AccsPerPtr;
  for (auto *MA : *KernelStmt) {
    if (!MA->isArrayKind())
      continue;
    auto &MAs = AccsPerPtr[MA->getScopArrayInfo()];
    MAs.push_back(MA);
  }
  
  errs() << "\ndata accesses:\n";
  for (auto &AccsPerPtrItem : AccsPerPtr) {
    auto *SAI = AccsPerPtrItem.first;
    auto &MAs = AccsPerPtrItem.second;
    printAccesses(OS, SAI, MAs);
  }
  
}

static isl_stat getSinglePiece(isl_set *Dom, isl_aff *Aff, void *User) {
  auto *UserPtr = static_cast<isl_aff **>(User);
  assert(*UserPtr == nullptr);
  *UserPtr = Aff;
  isl_set_free(Dom);
  return isl_stat_ok;
}

void LayerCondition::printVal(raw_ostream &OS,
                                 __isl_keep isl_val *Val) const {
  auto AI = APIntFromVal(isl_val_copy(Val));
  AI.print(OS, true);
}

void LayerCondition::printAff(raw_ostream &OS, __isl_keep isl_aff *Aff,
                                 bool IsDiv) const {
  auto First = true;

  auto *Val = isl_aff_get_constant_val(Aff);
  if (!isl_val_is_zero(Val)) {
    First = false;
    printVal(OS, Val);
  }
  isl_val_free(Val);

  int DimIn = isl_aff_dim(Aff, isl_dim_in);
  for (int dim = 0; dim < DimIn; dim++) {
    Val = isl_aff_get_coefficient_val(Aff, isl_dim_in, dim);
    if (isl_val_is_zero(Val)) {
      isl_val_free(Val);
      continue;
    }

    if (!First)
      OS << "+";
    First = false;

    if (isl_val_is_one(Val)) {
      isl_val_free(Val);
      OS << ((char)('i' + dim));
      continue;
    }

    printVal(OS, Val);
    isl_val_free(Val);
    OS << "*" << ((char)('i' + dim));
  }

  if (IsDiv)
    return;

  int DimDiv = isl_aff_dim(Aff, isl_dim_div);
  for (int dim = 0; dim < DimDiv; dim++) {
    auto *DivAff = isl_aff_get_div(Aff, dim);
    if (!First)
      OS << "+";

    Val = isl_aff_get_denominator_val(DivAff);
    DivAff = isl_aff_scale_val(DivAff, isl_val_copy(Val));

    OS << "floor(";
    printAff(OS, DivAff, true);
    isl_aff_free(DivAff);
    OS << ")";

    if (!isl_val_is_one(Val)) {
      OS << "/";
      printVal(OS, Val);
    }
    isl_val_free(Val);

    First = false;
  }
}

void LayerCondition::printAccesses(raw_ostream &OS, const ScopArrayInfo *SAI,
                                  MemoryAccessVec &MAs, bool IgnoreWrites) const {
  auto First = true;
  auto NumDims = SAI->getNumberOfDimensions();
  for (auto *MA : MAs) {
    if(MA->isWrite() && IgnoreWrites)
      // Ignore writes (for non write-allocate architectures)
      continue;

    if (First) {
      OS.indent(4) << getArrayName(SAI) << ":\n";
      First = false;
    }

    auto *AR = MA->getAccessRelation();
    auto *IP = isl_printer_to_str(isl_map_get_ctx(AR));
    isl_printer_print_map(IP, AR);
    OS << isl_printer_get_str(IP) << "\n";
    isl_printer_flush(IP); 
    //MA->getSubscript(0)->print(OS);
    OS << "\n";
    isl_printer_free(IP);
    auto *PwMAff = isl_pw_multi_aff_from_map(isl_map_copy(AR));
    OS.indent(8) << "- [";
    for (unsigned u = 0; u < NumDims; u++) {
      auto *PwAff = isl_pw_multi_aff_get_pw_aff(PwMAff, u);
      assert(isl_pw_aff_n_piece(PwAff) == 1);

      isl_aff *Aff = nullptr;
      isl_pw_aff_foreach_piece(PwAff, getSinglePiece, &Aff);
      isl_pw_aff_free(PwAff);

      printAff(OS, Aff);
      isl_aff_free(Aff);
      OS << ((u + 1 < NumDims) ? ", " : "]\n");
    }
    isl_pw_multi_aff_free(PwMAff);
    isl_map_free(AR);
  }
}

void LayerCondition::setFunction(Function &F) {
  ScopNumber++;

  if (&F == LastF)
    return;

  LastF = &F;
  ScopNumber = 0;
  collectVarNames(F);
}

bool LayerCondition::runOnScop(Scop &S) {
  printScop(errs(), S);
  return true;
}

void LayerCondition::getAnalysisUsage(AnalysisUsage &AU) const {
  ScopPass::getAnalysisUsage(AU);
  AU.addRequired<TargetTransformInfoWrapperPass>();
}

Pass *polly::createLayerConditionPass() { return new LayerCondition(); }

char LayerCondition::ID = 0;

INITIALIZE_PASS_BEGIN(LayerCondition, "polly-layercondition-pass",
                      "Polly - Analyzes SCoPs for layer conditions.", true, true);
INITIALIZE_PASS_DEPENDENCY(ScopInfoRegionPass);
INITIALIZE_PASS_DEPENDENCY(TargetTransformInfoWrapperPass);
INITIALIZE_PASS_END(LayerCondition, "polly-layercondition-pass",
                    "Polly - Analyzes SCoPs for layer conditions.", true, true);