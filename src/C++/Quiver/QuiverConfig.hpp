// Copyright (c) 2011-2013, Pacific Biosciences of California, Inc.
//
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted (subject to the limitations in the
// disclaimer below) provided that the following conditions are met:
//
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//
//  * Redistributions in binary form must reproduce the above
//    copyright notice, this list of conditions and the following
//    disclaimer in the documentation and/or other materials provided
//    with the distribution.
//
//  * Neither the name of Pacific Biosciences nor the names of its
//    contributors may be used to endorse or promote products derived
//    from this software without specific prior written permission.
//
// NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE
// GRANTED BY THIS LICENSE. THIS SOFTWARE IS PROVIDED BY PACIFIC
// BIOSCIENCES AND ITS CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
// WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
// OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
// USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
// OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
// SUCH DAMAGE.

// Author: David Alexander

#pragma once

#include <list>
#include <string>
#include <utility>

#include "Utils.hpp"

namespace ConsensusCore
{
    enum Move
    {
        INVALID_MOVE = 0x0,
        INCORPORATE  = 0x1,
        EXTRA        = 0x2,
        DELETE       = 0x4,
        MERGE        = 0x8,
        BASIC_MOVES  = (INCORPORATE | EXTRA | DELETE),
        ALL_MOVES    = (BASIC_MOVES | MERGE)
    };

    /// \brief The banding optimizations to be used by a recursor
    ///
    /// The recursors may implement a dynamic adjustment of the ScoreDiff:
    /// dynamicScoreDiff = max(ScoreDiff, min(10*ScoreDiff, dynamicScoreDiff - adjustment))
    /// where
    ///   adjustment = DynamicAdjustFactor * (columnBestScoreChange - DynamicAdjustOffset)
    ///   columnBestScoreChange = the difference between best scores obtained in the current and
    ///                           previous columns.
    ///
    struct BandingOptions
    {
        int DiagonalCross;
        float ScoreDiff;
        float DynamicAdjustFactor;
        float DynamicAdjustOffset;

        BandingOptions(int diagonalCross, float scoreDiff)
            : DiagonalCross(diagonalCross),
              ScoreDiff(scoreDiff),
              DynamicAdjustFactor(0.0f),
              DynamicAdjustOffset(0.0f)
        {}

        BandingOptions(int diagonalCross, float scoreDiff,
                       float dynamicAdjustFactor, float dynamicAdjustOffset)
            : DiagonalCross(diagonalCross),
              ScoreDiff(scoreDiff),
              DynamicAdjustFactor(dynamicAdjustFactor),
              DynamicAdjustOffset(dynamicAdjustOffset)
        {}
    };


    /// \brief A parameter vector for analysis using the QV model
    struct QvModelParams
    {
        float Match;
        float Mismatch;
        float MismatchS;
        float Branch;
        float BranchS;
        float DeletionN;
        float DeletionWithTag;
        float DeletionWithTagS;
        float Nce;
        float NceS;
        float Merge[4];
        float MergeS[4];

        //
        // Constructor for single merge rate and merge rate slope
        //
        QvModelParams(float Match,
                      float Mismatch,
                      float MismatchS,
                      float Branch,
                      float BranchS,
                      float DeletionN,
                      float DeletionWithTag,
                      float DeletionWithTagS,
                      float Nce,
                      float NceS,
                      float Merge,
                      float MergeS)
            : Match(Match)
            , Mismatch(Mismatch)
            , MismatchS(MismatchS)
            , Branch(Branch)
            , BranchS(BranchS)
            , DeletionN(DeletionN)
            , DeletionWithTag(DeletionWithTag)
            , DeletionWithTagS(DeletionWithTagS)
            , Nce(Nce)
            , NceS(NceS)
        {
            for (int base = 0; base < 4; base++)
            {
                this->Merge[base]  = Merge;
                this->MergeS[base] = MergeS;
            }
        }

        //
        // Constructor for per-channel merge rate and merge rate slope
        //
        QvModelParams(float Match,
                      float Mismatch,
                      float MismatchS,
                      float Branch,
                      float BranchS,
                      float DeletionN,
                      float DeletionWithTag,
                      float DeletionWithTagS,
                      float Nce,
                      float NceS,
                      float Merge_A,
                      float Merge_C,
                      float Merge_G,
                      float Merge_T,
                      float MergeS_A,
                      float MergeS_C,
                      float MergeS_G,
                      float MergeS_T)
            : Match(Match)
            , Mismatch(Mismatch)
            , MismatchS(MismatchS)
            , Branch(Branch)
            , BranchS(BranchS)
            , DeletionN(DeletionN)
            , DeletionWithTag(DeletionWithTag)
            , DeletionWithTagS(DeletionWithTagS)
            , Nce(Nce)
            , NceS(NceS)
        {
            this->Merge[0] = Merge_A;
            this->Merge[1] = Merge_C;
            this->Merge[2] = Merge_G;
            this->Merge[3] = Merge_T;
            this->MergeS[0] = MergeS_A;
            this->MergeS[1] = MergeS_C;
            this->MergeS[2] = MergeS_G;
            this->MergeS[3] = MergeS_T;
        }


        // Access to the array-stored params

        float Merge_A() const { return this->Merge[0]; }
        float Merge_C() const { return this->Merge[1]; }
        float Merge_G() const { return this->Merge[2]; }
        float Merge_T() const { return this->Merge[3]; }

        float MergeS_A() const { return this->MergeS[0]; }
        float MergeS_C() const { return this->MergeS[1]; }
        float MergeS_G() const { return this->MergeS[2]; }
        float MergeS_T() const { return this->MergeS[3]; }
    };


    struct QuiverConfig
    {
        const QvModelParams QvParams;
        const int MovesAvailable;
        const BandingOptions Banding;
        const float FastScoreThreshold;

        QuiverConfig(const QvModelParams& qvParams,
                     int movesAvailable,
                     const BandingOptions& bandingOptions,
                     float fastScoreThreshold);

        QuiverConfig(const QuiverConfig& qvConfig);
    };



    class QuiverConfigTable
    {
    private:
        typedef std::pair<const std::string, const QuiverConfig> QuiverConfigTableEntry;

        std::list<QuiverConfigTableEntry> table;

    public:
        QuiverConfigTable();

        typedef std::list<QuiverConfigTableEntry>::const_iterator const_iterator;

        bool insert(const std::string& name, const QuiverConfig& config);

        const QuiverConfig& at(const std::string& name) const throw(InvalidInputError);

#ifndef SWIG
        const_iterator begin() const;
        const_iterator end() const;
#endif
    };
}
