//===- IslTest.cpp ----------------------------------------------------===//
//
//                     The LLVM Compiler Infrastructure
//
// This file is distributed under the University of Illinois Open Source
// License. See LICENSE.TXT for details.
//
//===----------------------------------------------------------------------===//

#include "polly/Support/GICHelper.h"
#include "gtest/gtest.h"
#include "isl/val.h"

using namespace llvm;
using namespace polly;

namespace {

TEST(Isl, APIntToIslVal) {
  isl_ctx *IslCtx = isl_ctx_alloc();

  {
    APInt APZero(1, 0, true);
    auto *IslZero = isl_valFromAPInt(IslCtx, APZero, true);
    EXPECT_EQ(isl_bool_true, isl_val_is_zero(IslZero));
    isl_val_free(IslZero);
  }

  {
    APInt APNOne(1, -1, true);
    auto *IslNOne = isl_valFromAPInt(IslCtx, APNOne, true);
    EXPECT_EQ(isl_bool_true, isl_val_is_negone(IslNOne));
    isl_val_free(IslNOne);
  }

  {
    APInt APZero(1, 0, false);
    auto *IslZero = isl_valFromAPInt(IslCtx, APZero, false);
    EXPECT_EQ(isl_bool_true, isl_val_is_zero(IslZero));
    isl_val_free(IslZero);
  }

  {
    APInt APOne(1, 1, false);
    auto *IslOne = isl_valFromAPInt(IslCtx, APOne, false);
    EXPECT_EQ(isl_bool_true, isl_val_is_one(IslOne));
    isl_val_free(IslOne);
  }

  {
    APInt APNTwo(2, -2, true);
    auto *IslNTwo = isl_valFromAPInt(IslCtx, APNTwo, true);
    auto *IslNTwoCmp = isl_val_int_from_si(IslCtx, -2);
    EXPECT_EQ(isl_bool_true, isl_val_eq(IslNTwo, IslNTwoCmp));
    isl_val_free(IslNTwo);
    isl_val_free(IslNTwoCmp);
  }

  {
    APInt APNOne(32, -1, true);
    auto *IslNOne = isl_valFromAPInt(IslCtx, APNOne, true);
    EXPECT_EQ(isl_bool_true, isl_val_is_negone(IslNOne));
    isl_val_free(IslNOne);
  }

  {
    APInt APZero(32, 0, false);
    auto *IslZero = isl_valFromAPInt(IslCtx, APZero, false);
    EXPECT_EQ(isl_bool_true, isl_val_is_zero(IslZero));
    isl_val_free(IslZero);
  }

  {
    APInt APOne(32, 1, false);
    auto *IslOne = isl_valFromAPInt(IslCtx, APOne, false);
    EXPECT_EQ(isl_bool_true, isl_val_is_one(IslOne));
    isl_val_free(IslOne);
  }

  {
    APInt APTwo(32, 2, false);
    auto *IslTwo = isl_valFromAPInt(IslCtx, APTwo, false);
    EXPECT_EQ(0, isl_val_cmp_si(IslTwo, 2));
    isl_val_free(IslTwo);
  }

  {
    APInt APNOne(32, (1ull << 32) - 1, false);
    auto *IslNOne = isl_valFromAPInt(IslCtx, APNOne, false);
    auto *IslRef = isl_val_int_from_ui(IslCtx, (1ull << 32) - 1);
    EXPECT_EQ(isl_bool_true, isl_val_eq(IslNOne, IslRef));
    isl_val_free(IslNOne);
    isl_val_free(IslRef);
  }

  {
    APInt APLarge(130, 2, false);
    APLarge = APLarge.shl(70);
    auto *IslLarge = isl_valFromAPInt(IslCtx, APLarge, false);
    auto *IslRef = isl_val_int_from_ui(IslCtx, 71);
    IslRef = isl_val_2exp(IslRef);
    EXPECT_EQ(isl_bool_true, isl_val_eq(IslLarge, IslRef));
    isl_val_free(IslLarge);
    isl_val_free(IslRef);
  }

  isl_ctx_free(IslCtx);
}

TEST(Isl, IslValToAPInt) {
  isl_ctx *IslCtx = isl_ctx_alloc();

  {
    auto *IslNOne = isl_val_int_from_si(IslCtx, -1);
    auto APNOne = APIntFromVal(IslNOne);
    // Compare with the two's complement of -1 in a 1-bit integer.
    EXPECT_EQ(1, APNOne);
    EXPECT_EQ(1u, APNOne.getBitWidth());
  }

  {
    auto *IslNTwo = isl_val_int_from_si(IslCtx, -2);
    auto APNTwo = APIntFromVal(IslNTwo);
    // Compare with the two's complement of -2 in a 2-bit integer.
    EXPECT_EQ(2, APNTwo);
    EXPECT_EQ(2u, APNTwo.getBitWidth());
  }

  {
    auto *IslNThree = isl_val_int_from_si(IslCtx, -3);
    auto APNThree = APIntFromVal(IslNThree);
    // Compare with the two's complement of -3 in a 3-bit integer.
    EXPECT_EQ(5, APNThree);
    EXPECT_EQ(3u, APNThree.getBitWidth());
  }

  {
    auto *IslNFour = isl_val_int_from_si(IslCtx, -4);
    auto APNFour = APIntFromVal(IslNFour);
    // Compare with the two's complement of -4 in a 3-bit integer.
    EXPECT_EQ(4, APNFour);
    EXPECT_EQ(3u, APNFour.getBitWidth());
  }

  {
    auto *IslZero = isl_val_int_from_ui(IslCtx, 0);
    auto APZero = APIntFromVal(IslZero);
    EXPECT_EQ(0, APZero);
    EXPECT_EQ(1u, APZero.getBitWidth());
  }

  {
    auto *IslOne = isl_val_int_from_ui(IslCtx, 1);
    auto APOne = APIntFromVal(IslOne);
    EXPECT_EQ(1, APOne);
    EXPECT_EQ(2u, APOne.getBitWidth());
  }

  {
    auto *IslTwo = isl_val_int_from_ui(IslCtx, 2);
    auto APTwo = APIntFromVal(IslTwo);
    EXPECT_EQ(2, APTwo);
    EXPECT_EQ(3u, APTwo.getBitWidth());
  }

  {
    auto *IslThree = isl_val_int_from_ui(IslCtx, 3);
    auto APThree = APIntFromVal(IslThree);
    EXPECT_EQ(3, APThree);
    EXPECT_EQ(3u, APThree.getBitWidth());
  }

  {
    auto *IslFour = isl_val_int_from_ui(IslCtx, 4);
    auto APFour = APIntFromVal(IslFour);
    EXPECT_EQ(4, APFour);
    EXPECT_EQ(4u, APFour.getBitWidth());
  }

  {
    auto *IslNOne = isl_val_int_from_ui(IslCtx, (1ull << 32) - 1);
    auto APNOne = APIntFromVal(IslNOne);
    EXPECT_EQ((1ull << 32) - 1, APNOne);
    EXPECT_EQ(33u, APNOne.getBitWidth());
  }

  {
    auto *IslLargeNum = isl_val_int_from_ui(IslCtx, 60);
    IslLargeNum = isl_val_2exp(IslLargeNum);
    IslLargeNum = isl_val_sub_ui(IslLargeNum, 1);
    auto APLargeNum = APIntFromVal(IslLargeNum);
    EXPECT_EQ((1ull << 60) - 1, APLargeNum);
    EXPECT_EQ(61u, APLargeNum.getBitWidth());
  }

  {
    auto *IslExp = isl_val_int_from_ui(IslCtx, 500);
    auto *IslLargePow2 = isl_val_2exp(IslExp);
    auto APLargePow2 = APIntFromVal(IslLargePow2);
    EXPECT_TRUE(APLargePow2.isPowerOf2());
    EXPECT_EQ(502u, APLargePow2.getBitWidth());
    EXPECT_EQ(502u, APLargePow2.getMinSignedBits());
  }

  {
    auto *IslExp = isl_val_int_from_ui(IslCtx, 500);
    auto *IslLargeNPow2 = isl_val_neg(isl_val_2exp(IslExp));
    auto APLargeNPow2 = APIntFromVal(IslLargeNPow2);
    EXPECT_EQ(501u, APLargeNPow2.getBitWidth());
    EXPECT_EQ(501u, APLargeNPow2.getMinSignedBits());
    EXPECT_EQ(500, (-APLargeNPow2).exactLogBase2());
  }

  {
    auto *IslExp = isl_val_int_from_ui(IslCtx, 512);
    auto *IslLargePow2 = isl_val_2exp(IslExp);
    auto APLargePow2 = APIntFromVal(IslLargePow2);
    EXPECT_TRUE(APLargePow2.isPowerOf2());
    EXPECT_EQ(514u, APLargePow2.getBitWidth());
    EXPECT_EQ(514u, APLargePow2.getMinSignedBits());
  }

  {
    auto *IslExp = isl_val_int_from_ui(IslCtx, 512);
    auto *IslLargeNPow2 = isl_val_neg(isl_val_2exp(IslExp));
    auto APLargeNPow2 = APIntFromVal(IslLargeNPow2);
    EXPECT_EQ(513u, APLargeNPow2.getBitWidth());
    EXPECT_EQ(513u, APLargeNPow2.getMinSignedBits());
    EXPECT_EQ(512, (-APLargeNPow2).exactLogBase2());
  }

  isl_ctx_free(IslCtx);
}

TEST(Isl, Foreach) {
  std::unique_ptr<isl_ctx, decltype(&isl_ctx_free)> Ctx(isl_ctx_alloc(),
                                                        &isl_ctx_free);

  auto MapSpace = give(isl_space_alloc(Ctx.get(), 0, 1, 1));
  auto TestBMap = give(isl_basic_map_universe(MapSpace.copy()));
  TestBMap = give(isl_basic_map_fix_si(TestBMap.take(), isl_dim_in, 0, 0));
  TestBMap = give(isl_basic_map_fix_si(TestBMap.take(), isl_dim_out, 0, 0));
  auto TestMap = give(isl_map_from_basic_map(TestBMap.copy()));
  auto TestUMap = give(isl_union_map_from_map(TestMap.copy()));

  auto SetSpace = give(isl_space_set_alloc(Ctx.get(), 0, 1));
  auto TestBSet =
      give(isl_basic_set_from_point(isl_point_zero(SetSpace.copy())));
  auto TestSet = give(isl_set_from_basic_set(TestBSet.copy()));
  auto TestUSet = give(isl_union_set_from_set(TestSet.copy()));

  {
    auto NumBMaps = 0;
    foreachElt(TestMap, [&](IslPtr<isl_basic_map> BMap) {
      EXPECT_EQ(isl_bool_true,
                isl_basic_map_is_equal(BMap.keep(), TestBMap.keep()));
      NumBMaps++;
    });
    EXPECT_EQ(1, NumBMaps);
  }

  {
    auto NumBSets = 0;
    foreachElt(TestSet, [&](IslPtr<isl_basic_set> BSet) {
      EXPECT_EQ(isl_bool_true,
                isl_basic_set_is_equal(BSet.keep(), TestBSet.keep()));
      NumBSets++;
    });
    EXPECT_EQ(1, NumBSets);
  }

  {
    auto NumMaps = 0;
    foreachElt(TestUMap, [&](IslPtr<isl_map> Map) {
      EXPECT_EQ(isl_bool_true, isl_map_is_equal(Map.keep(), TestMap.keep()));
      NumMaps++;
    });
    EXPECT_EQ(1, NumMaps);
  }

  {
    auto NumSets = 0;
    foreachElt(TestUSet, [&](IslPtr<isl_set> Set) {
      EXPECT_EQ(isl_bool_true, isl_set_is_equal(Set.keep(), TestSet.keep()));
      NumSets++;
    });
    EXPECT_EQ(1, NumSets);
  }

  {
    auto UPwAff = give(isl_union_pw_aff_val_on_domain(TestUSet.copy(),
                                                      isl_val_zero(Ctx.get())));
    auto NumPwAffs = 0;
    foreachElt(UPwAff, [&](IslPtr<isl_pw_aff> PwAff) {
      EXPECT_EQ(isl_bool_true, isl_pw_aff_is_cst(PwAff.keep()));
      NumPwAffs++;
    });
    EXPECT_EQ(1, NumPwAffs);
  }

  {
    auto NumBMaps = 0;
    EXPECT_EQ(isl_stat_error,
              foreachEltWithBreak(
                  TestMap, [&](IslPtr<isl_basic_map> BMap) -> isl_stat {
                    EXPECT_EQ(
                        isl_bool_true,
                        isl_basic_map_is_equal(BMap.keep(), TestBMap.keep()));
                    NumBMaps++;
                    return isl_stat_error;
                  }));
    EXPECT_EQ(1, NumBMaps);
  }

  {
    auto NumMaps = 0;
    EXPECT_EQ(
        isl_stat_error,
        foreachEltWithBreak(TestUMap, [&](IslPtr<isl_map> Map) -> isl_stat {
          EXPECT_EQ(isl_bool_true,
                    isl_map_is_equal(Map.keep(), TestMap.keep()));
          NumMaps++;
          return isl_stat_error;
        }));
    EXPECT_EQ(1, NumMaps);
  }

  {
    auto TestPwAff =
        give(isl_pw_aff_val_on_domain(TestSet.copy(), isl_val_zero(Ctx.get())));
    auto NumPieces = 0;
    foreachPieceWithBreak(
        TestPwAff,
        [&](IslPtr<isl_set> Domain, IslPtr<isl_aff> Aff) -> isl_stat {
          EXPECT_EQ(isl_bool_true,
                    isl_set_is_equal(Domain.keep(), TestSet.keep()));
          EXPECT_EQ(isl_bool_true, isl_aff_is_cst(Aff.keep()));
          NumPieces++;
          return isl_stat_error;
        });
    EXPECT_EQ(1, NumPieces);
  }
}

} // anonymous namespace
