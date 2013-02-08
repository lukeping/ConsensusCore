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
    struct BandingOptions
    {
        int DiagonalCross;
        float ScoreDiff;

        BandingOptions(int diagonalCross, float scoreDiff)
            : DiagonalCross(diagonalCross),
              ScoreDiff(scoreDiff)
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


        static const QvModelParams Untrained();
        static const QvModelParams C2_130();
        static const QvModelParams Default();
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
    };
}
