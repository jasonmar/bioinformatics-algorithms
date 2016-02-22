/*
 * Bioinformatics Algorithms in Scala
 * Copyright (C) 2016  Jason Mar
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package bioinf.comparinggenomes

import bioinf.UnitSpec
import bioinf.ComparingGenomes._

class ComparingGenomesSpec extends UnitSpec{
  "ComparingGenomes" should "score global alignmeent correctly" in {
    scoreAlignment("AAA","AAA",fn_BLOSUM62) should be (12)
  }

  it should "calculate BLOSUM62 score " in {
    fn_BLOSUM62('Y','Y') should be (7)
  }

  it should "calculate global alignment" in {
    val v = "PLEASANTLY"
    val w = "MEANLY"
    val sigma = 5
    val res = globalAlignment(v,w,sigma,fn_BLOSUM62)
    scoreGlobalAlignment(res) should be (8)
  }

  it should "calculate PAM250 score " in {
    fn_PAM250('Y','Y') should be (10)
  }

  it should "calculate local alignment" in {
    val v = "MEANLY"
    val w = "PENALTY"
    val sigma = 5
    val res = localAlignment(v,w,sigma,fn_PAM250)
    scoreLocalAlignment(res) should be (15)
  }

  it should "calculate affine gap alignment" in {
    val v = "PRTEINS"
    val w = "PRTWPSEIN"
    val sigma = 11
    val eta = 1
    val res = alignmentWithAffineGap(v, w, sigma, eta, fn_BLOSUM62)
    val score = scoreAffineGapAlignment(res)
    res.score should be (score)
    score should be (8)

  }

  it should "calculate longer affine gap alignment" in {
    val v = "KRYINALREEAYHCNNIHLFARCDDQRDNNYTQCTGYMGGVYYKWQFLIIQLYLCHSKVYAMSQMVVTPLRVTMYIV"
    val w = "KRALRNNIHLFARCDDPRDNNYTACTGYMGDVYYKWQFMIIHLYLCHSFQVYAMSQMVVEPLRVTMYEV"
    val sigma = 11
    val eta = 1
    val res = alignmentWithAffineGap(v, w, sigma, eta, fn_BLOSUM62)
    val score = scoreAffineGapAlignment(res)
    res.score should be (score)
    score should be (288)
    val print = res.print.split(System.lineSeparator())
    print(1) should be ("KRYINALREEAYHCNNIHLFARCDDQRDNNYTQCTGYMGGVYYKWQFLIIQLYLCHS-KVYAMSQMVVTPLRVTMYIV")
    print(2) should be ("KR---ALR------NNIHLFARCDDPRDNNYTACTGYMGDVYYKWQFMIIHLYLCHSFQVYAMSQMVVEPLRVTMYEV")

  }

}
