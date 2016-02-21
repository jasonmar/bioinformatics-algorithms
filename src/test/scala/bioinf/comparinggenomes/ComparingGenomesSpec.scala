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
    scoreGlobalAlignment("AAA","AAA") should be (12)
  }

  it should "calculate BLOSUM62 score " in {
    BLOSUM62score('Y','Y') should be (7)
  }

  it should "calculate global alignment" in {
    val v = "PLEASANTLY"
    val w = "MEANLY"
    val sigma = 5
    val res = globalAlignment(v,w,sigma)
    scoreGlobalAlignment(res) should be (8)
  }
}
